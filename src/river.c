#include "river.h"
#include "adaptive.h"

static void cell_interpolated_face_values (FttCell * cell,
					   const GfsRiverParams * par)
{
  FttComponent c;
  guint v;

  for (c = 0; c < FTT_DIMENSION; c++)
    for (v = 0; v < GFS_RIVER_NVAR; v++) {
      gdouble val = GFS_VALUE (cell, par->v[v]);
      gdouble g = (* par->gradient) (cell, c, par->v[v]->i)/2.;
      
      GFS_VALUE (cell, par->fv[2*c][v])     = val + g;
      GFS_VALUE (cell, par->fv[2*c + 1][v]) = val - g;
    }
}

static void boundary_face_values (FttCell * cell,
				  const GfsRiverParams * p)
{
  FttCell * neighbor = ftt_cell_neighbor (cell, p->d);
  if (neighbor != NULL) {
    FttDirection d = FTT_OPPOSITE_DIRECTION (p->d);
    guint v = 0;
    g_assert (GFS_CELL_IS_BOUNDARY (neighbor));
    for (v = 0; v < GFS_RIVER_NVAR; v++)
      /* fixme: what about g? */
      GFS_VALUE (neighbor, p->fv[d][v]) = GFS_VALUE (neighbor, p->v[v]);
  }
}

static void flux (const gdouble * u, gdouble g, gdouble * f)
{
  f[0] = u[0]*u[1];                                     /* h*u */
  f[1] = u[0]*u[1]*u[1] + g*(u[0]*u[0] - u[3]*u[3])/2.; /* h*u*u + g*(h*h - zb*zb)/2 */
  f[2] = u[0]*u[1]*u[2];                                /* h*u*v */
}

static gdouble min (gdouble a, gdouble b)
{
  return a < b ? a : b;
}

static gdouble max (gdouble a, gdouble b)
{
  return a > b ? a : b;
}

/*
 * uL: left state vector [h,u,v,zb].
 * uR: right state vector.
 * g: acceleration of gravity.
 * f: flux vector.
 *
 * Fills @f by solving an approximate Riemann problem using the HLLC
 * scheme. See e.g. Liang, Borthwick, Stelling, IJNMF, 2004.
 */
static void riemann (const gdouble * uL, const gdouble * uR,
		     gdouble g,
		     gdouble * f)
{
#if 0
  gdouble fL[3], fR[3];
  flux (uL, g, fL);
  flux (uR, g, fR);
  guint i;
  for (i = 0; i < 3; i++)
    f[i] = (fL[i] + fR[i])/2.;
  return;
#endif
  gdouble cL = sqrt (g*uL[0]), cR = sqrt (g*uR[0]);
  gdouble ustar = (uL[1] + uR[1])/2. + cL - cR;
  gdouble cstar = (cL + cR)/2. + (uL[1] - uR[1])/4.;
  g_assert (cstar >= 0.);
  gdouble SL = uL[0] == 0. ? uR[1] - 2.*cR : min (uL[1] - cL, ustar - cstar);
  gdouble SR = uR[0] == 0. ? uL[1] + 2.*cL : max (uR[1] + cR, ustar + cstar);

  if (0. <= SL)
    flux (uL, g, f);
  else if (0. >= SR)
    flux (uR, g, f);
  else {
    gdouble fL[3], fR[3];
    flux (uL, g, fL);
    flux (uR, g, fR);
    guint i;
    for (i = 0; i < 2; i++)
      f[i] = (SR*fL[i] - SL*fR[i] + SL*SR*(uR[i] - uL[i]))/(SR - SL);

    gdouble SM = ((SL*uR[0]*(uR[1] - SR) - SR*uL[0]*(uL[1] - SL))/
		  (uR[0]*(uR[1] - SR) - uL[0]*(uL[1] - SL)));
    if (SL <= 0. && 0. <= SM)
      f[2] = uL[2]*f[0];
    else if (SM <= 0. && 0. <= SR)
      f[2] = uR[2]*f[0];
    else
      g_assert_not_reached ();
  }
}

static void face_fluxes (FttCellFace * face, GfsRiverParams * p)
{
  gdouble zbL = 0.;
  gdouble uL[4], uR[4], f[3];
  
  uL[0] = GFS_VALUE (face->cell, p->fv[face->d][0]) - zbL; /* h = eta - zb */
  g_assert (uL[0] > 0.);
  uL[1] = GFS_VALUE (face->cell, p->fv[face->d][1])/uL[0]; /* u = uh/h */
  uL[2] = GFS_VALUE (face->cell, p->fv[face->d][2])/uL[0]; /* v = vh/h */
  uL[3] = zbL;

  g_assert (ftt_face_type (face) == FTT_FINE_FINE);
  gdouble zbR = 0.;
  FttDirection d = FTT_OPPOSITE_DIRECTION (face->d);
  uR[0] = GFS_VALUE (face->neighbor, p->fv[d][0]) - zbR; /* h = eta - zb */
  g_assert (uR[0] > 0.);
  uR[1] = GFS_VALUE (face->neighbor, p->fv[d][1])/uR[0]; /* u = uh/h */
  uR[2] = GFS_VALUE (face->neighbor, p->fv[d][2])/uR[0]; /* v = vh/h */
  uR[3] = zbR;

  riemann (uL, uR, p->g, f);

  gdouble dt = gfs_domain_face_fraction (p->v[0]->domain, face)*p->dt/ftt_cell_size (face->cell);
  if (!FTT_FACE_DIRECT (face))
    dt = - dt;
  GFS_VALUE (face->cell, p->flux[0]) -= dt*f[0];
  GFS_VALUE (face->cell, p->flux[1]) -= dt*f[1];
  GFS_VALUE (face->cell, p->flux[2]) -= dt*f[2];

  g_assert (ftt_face_type (face) == FTT_FINE_FINE);
  GFS_VALUE (face->neighbor, p->flux[0]) += dt*f[0];
  GFS_VALUE (face->neighbor, p->flux[1]) += dt*f[1];
  GFS_VALUE (face->neighbor, p->flux[2]) += dt*f[2];
}

/* GfsRiver: Object */

static void reset_fluxes (FttCell * cell, GfsRiverParams * p)
{
  guint v;
  for (v = 0; v < GFS_RIVER_NVAR; v++)
    GFS_VALUE (cell, p->flux[v]) = 0.;
}

static void river_run (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsRiverParams * p = &GFS_RIVER (sim)->p;

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  p->g = sim->physical_params.g;
  p->gradient = sim->advection_params.gradient;
  //  p->gradient = gfs_center_minmod_gradient;

  gfs_simulation_set_timestep (sim);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    p->dt = sim->advection_params.dt;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) cell_interpolated_face_values, p);
    guint v;
#if 1
    for (p->d = 0; p->d < FTT_NEIGHBORS; p->d++)
      gfs_domain_cell_traverse_boundary (domain, p->d, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
					 (FttCellTraverseFunc) boundary_face_values, p);
#else    
    /* fixme: works only for periodic BC */
    FttDirection d;
    for (d = 0; d < FTT_NEIGHBORS; d++)
      for (v = 0; v < GFS_RIVER_NVAR; v++)
	gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, p->fv[d][v]);
#endif
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) reset_fluxes, p);
    gfs_domain_face_traverse (domain, FTT_X,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) face_fluxes, p);
    for (v = 0; v < GFS_RIVER_NVAR; v++) {
      GfsAdvectionParams par;
      par.v = p->v[v];
      par.fv = p->flux[v];
      gfs_domain_traverse_merged (domain, (GfsMergedTraverseFunc) gfs_advection_update, &par);
      gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, par.v);
    }

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gfs_simulation_adapt (sim);

    sim->time.t = sim->tnext;
    sim->time.i++;

    gfs_simulation_set_timestep (sim);

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);
}

static void minimum_cfl (FttCell * cell, GfsRiverParams * p)
{
  gdouble size = ftt_cell_size (cell);
  gdouble eta = GFS_VALUE (cell, p->v[0]);
  gdouble zb = 0.;
  gdouble h = eta - zb;
  if (h > 0.) {
    gdouble uh = fabs (GFS_VALUE (cell, p->v[1]));
    gdouble c = sqrt (p->g*h);
    gdouble cfl = size/(uh/h + c);
    if (cfl < p->cfl)
      p->cfl = cfl;

    gdouble vh = fabs (GFS_VALUE (cell, p->v[2]));
    cfl = size/(vh/h + c);
    if (cfl < p->cfl)
      p->cfl = cfl;
  }
}

static gdouble river_cfl (GfsSimulation * sim)
{
  GfsRiverParams * p = &GFS_RIVER (sim)->p;
  p->cfl = G_MAXDOUBLE;
  gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) minimum_cfl, p);
  gfs_all_reduce (domain, p->cfl, MPI_DOUBLE, MPI_MIN);
  return p->cfl;
}

static void river_class_init (GfsSimulationClass * klass)
{
  klass->run = river_run;
  klass->cfl = river_cfl;
}

static void river_init (GfsRiver * c)
{
  GfsDomain * domain = GFS_DOMAIN (c);
  GfsRiverParams * p = &c->p;
  FttDirection d;

  p->v[0] = gfs_variable_from_name (domain->variables, "P");
  p->v[1] = gfs_variable_from_name (domain->variables, "U");
  p->v[2] = gfs_variable_from_name (domain->variables, "V");

  for (d = 0; d < FTT_NEIGHBORS; d++) {
    gchar * name = g_strdup_printf ("fvP%d", d);
    p->fv[d][0] = gfs_domain_add_variable (domain, name, NULL);
    g_free (name);
    name = g_strdup_printf ("fvU%d", d);
    p->fv[d][1] = gfs_domain_add_variable (domain, name, NULL);
    g_free (name);
    name = g_strdup_printf ("fvV%d", d);
    p->fv[d][2] = gfs_domain_add_variable (domain, name, NULL);
    g_free (name);
  }
  
  p->flux[0] = gfs_domain_add_variable (domain, "fluxP", NULL);
  p->flux[1] = gfs_domain_add_variable (domain, "fluxU", NULL);
  p->flux[2] = gfs_domain_add_variable (domain, "fluxV", NULL);
}

GfsSimulationClass * gfs_river_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_river_info = {
      "GfsRiver",
      sizeof (GfsRiver),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) river_class_init,
      (GtsObjectInitFunc) river_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()), 
				  &gfs_river_info);
  }

  return klass;
}

/* GfsBcSubcritical: Object */

static void subcritical (FttCellFace * f, GfsBc * b)
{
  gdouble hb = gfs_function_face_value (GFS_BC_VALUE (b)->val, f);
  GfsRiver * river = GFS_RIVER (b->v->domain);
  gdouble hi = GFS_VALUE (f->neighbor, river->p.v[0]);

  GFS_VALUE (f->cell, b->v) = GFS_VALUE (f->neighbor, b->v) -
    2.*hi*(sqrt (river->p.g*hi) - sqrt (river->p.g*hb));
}

static void bc_subcritical_read (GtsObject ** o, GtsFile * fp)
{
  GfsBc * bc = GFS_BC (*o);

  if (GTS_OBJECT_CLASS (gfs_bc_subcritical_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_bc_subcritical_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (!GFS_IS_RIVER (bc->v->domain)) {
    gts_file_error (fp, "GfsBcSubcritical only makes sense for GfsRiver simulations");
    return;
  }

  gfs_function_set_units (GFS_BC_VALUE (bc)->val, 1.);
}

static void gfs_bc_subcritical_init (GfsBc * object)
{
  object->bc =  (FttFaceTraverseFunc) subcritical;
}

static void gfs_bc_subcritical_class_init (GtsObjectClass * klass)
{
  klass->read = bc_subcritical_read;
}

GfsBcClass * gfs_bc_subcritical_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_subcritical_info = {
      "GfsBcSubcritical",
      sizeof (GfsBcValue),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_subcritical_class_init,
      (GtsObjectInitFunc) gfs_bc_subcritical_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_value_class ()),
				  &gfs_bc_subcritical_info);
  }

  return klass;
}
