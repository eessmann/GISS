#include "river.h"
#include "adaptive.h"

#define DRY 1e-6

/* GfsRiver: Object */

static void flux (const gdouble * u, gdouble g, gdouble * f)
{
  f[0] = u[0]*u[1];                                     /* h*u */
  f[1] = u[0]*u[1]*u[1] + g*(u[0]*u[0] - u[3]*u[3])/2.; /* h*u*u + g*(h*h - zb*zb)/2 */
  f[2] = u[0]*u[1]*u[2];                                /* h*u*v */
}

#define matrix_vector_product(C, M, V) ((C)[0] = (M)[0][0]*(V)[0] + (M)[0][1]*(V)[1] + (M)[0][2]*(V)[2],\
                                        (C)[1] = (M)[1][0]*(V)[0] + (M)[1][1]*(V)[1] + (M)[1][2]*(V)[2],\
                                        (C)[2] = (M)[2][0]*(V)[0] + (M)[2][1]*(V)[1] + (M)[2][2]*(V)[2])

static void riemann_roe (const gdouble * uL, const gdouble * uR,
			 gdouble g,
			 gdouble * f)
{
  GtsVector V, alpha, beta, lambda;
  gdouble Hl, Hr, H, u, v, sqrtgH;
  gdouble cpsi, cU, cV;
  GtsMatrix * M, * IM;

  Hl = uL[0];
  Hr = uR[0];

  if (Hr <= 0. || Hl <= 0.) {
    g_warning ("Hr: %g Hl: %g\n", Hr, Hl);
    g_assert_not_reached ();
  }

  Hl = 1./sqrt (Hl);
  Hr = 1./sqrt (Hr);

  cpsi = (Hl*uL[0] + Hr*uR[0])/(Hl + Hr);
  cU = (Hl*uL[0]*uL[1] + Hr*uR[0]*uR[1])/(Hl + Hr);
  cV = (Hl*uL[0]*uL[2] + Hr*uR[0]*uR[2])/(Hl + Hr);

  H = (1./Hl + 1./Hr)/(Hl + Hr);
  g_assert (H > 0.0);
  u = cU/H;
  v = cV/H;
  sqrtgH = sqrt (g*H);

  lambda[0] = u - sqrtgH;
  lambda[1] = u;
  lambda[2] = u + sqrtgH;

  M = gts_matrix_new (1.,                    0.,          1., 0.,
		      u - sqrtgH,            0.,  u + sqrtgH, 0.,
		      v,                     1.,          v,  0.,
		      0., 0., 0., 0.);
#if 0
  IM = gts_matrix3_inverse (M);
  g_assert (IM != NULL);
#else
  IM = gts_matrix_new ((sqrtgH + u)/(2.*sqrtgH), -1./(2.*sqrtgH), 0., 0.,
		       - v, 0., 1., 0.,
		       (sqrtgH - u)/(2.*sqrtgH),  1./(2.*sqrtgH), 0.,  0.,
		       0., 0., 0., 0.);
#endif

  V[0] = uR[0] - uL[0];
  V[1] = uR[0]*uR[1] - uL[0]*uL[1];
  V[2] = uR[0]*uR[2] - uL[0]*uL[2];
  matrix_vector_product (alpha, IM, V);
  gts_matrix_destroy (IM);

  beta[0] = alpha[0]*MIN (lambda[0], 0.0);
  beta[1] = alpha[1]*MIN (lambda[1], 0.0);
  beta[2] = alpha[2]*MIN (lambda[2], 0.0);
  matrix_vector_product (V, M, beta);

  flux (uL, g, f);
  f[0] += V[0];
  f[1] += V[1];
  f[2] += V[2];

#if 0
  /* Flux right -> left */

  u = - u;

  lambda[0] = u - sqrtgH;
  lambda[1] = u;
  lambda[2] = u + sqrtgH;

  gts_matrix_assign (M,
		     1.,                    0.,               1., 0.,
		     u - sqrtgH,            0.,       u + sqrtgH, 0.,
		     v,                     1.,                v, 0.,
		     0., 0., 0., 0.);
  beta[0] = - alpha[2]*MIN (lambda[0], 0.0);
  beta[1] =   alpha[1]*MIN (lambda[1], 0.0);
  beta[2] = - alpha[0]*MIN (lambda[2], 0.0);
  matrix_vector_product (V, M, beta);

  gdouble f2[3];  
  flux (uR, g, f2);
  f2[0] = - f2[0] + V[0];
  f2[1] = - f2[1] + V[1];
  f2[2] = - f2[2] + V[2];

  g_assert (f2[0] + f[0] < 1e-10);
#endif

  gts_matrix_destroy (M);
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
static void riemann_hllc (const gdouble * uL, const gdouble * uR,
			  gdouble g,
			  gdouble * f)
{
  gdouble cL = sqrt (g*uL[0]), cR = sqrt (g*uR[0]);
  gdouble ustar = (uL[1] + uR[1])/2. + cL - cR;
  gdouble cstar = (cL + cR)/2. + (uL[1] - uR[1])/4.;
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
    f[0] = (SR*fL[0] - SL*fR[0] + SL*SR*(uR[0] - uL[0]))/(SR - SL);
    f[1] = (SR*fL[1] - SL*fR[1] + SL*SR*(uR[0]*uR[1] - uL[0]*uL[1]))/(SR - SL);
    gdouble SM = ((SL*uR[0]*(uR[1] - SR) - SR*uL[0]*(uL[1] - SL))/
		  (uR[0]*(uR[1] - SR) - uL[0]*(uL[1] - SL)));
    if (SL <= 0. && 0. <= SM)
      f[2] = uL[2]*f[0];
    else if (SM <= 0. && 0. <= SR)
      f[2] = uR[2]*f[0];
    else {
      fprintf (stderr, "L: %g %g %g R: %g %g %g\n",
	       uL[0], uL[1], uL[2],
	       uR[0], uR[1], uR[2]);
      fprintf (stderr, "SL: %g SR: %g SM: %g\n", SL, SR, SM);
      g_assert_not_reached ();
    }
  }
}

#define U 1
#define V 2

typedef struct {
  FttComponent u;
  gdouble du;
  FttComponent v;
  gdouble dv;
} Sym;

static void face_fluxes (FttCellFace * face, GfsRiver * r)
{
  static Sym sym[4] = {
    {U,  1., V,  1.},
    {U, -1., V, -1.},
    {V,  1., U, -1.},
    {V, -1., U,  1.}
  };
  Sym * s = &sym[face->d];
  gdouble etaL = (GFS_VALUE (face->cell, r->v1[0]) + 
		  s->du*GFS_VALUE (face->cell, r->dv[face->d/2][0]));
  gdouble zbL = (GFS_VALUE (face->cell, r->v[3]) + 
		 s->du*GFS_VALUE (face->cell, r->dv[face->d/2][3]));
  gdouble zbR = (GFS_VALUE (face->neighbor, r->v[3]) -
		 s->du*GFS_VALUE (face->neighbor, r->dv[face->d/2][3]));
  gdouble zbLR = MAX (zbL, zbR);
  gdouble uL[4], uR[4], f[3];

  if (etaL > DRY) {
    uL[1] = s->du*(GFS_VALUE (face->cell, r->v1[s->u]) +
		   s->du*GFS_VALUE (face->cell, r->dv[face->d/2][s->u]))/etaL; /* u = uh/h */
    uL[2] = s->dv*(GFS_VALUE (face->cell, r->v1[s->v]) +
		   s->du*GFS_VALUE (face->cell, r->dv[face->d/2][s->v]))/etaL; /* v = vh/h */
  }
  else
    uL[1] = uL[2] = 0.;
  uL[0] = MAX (0., etaL + zbL - zbLR);
  uL[3] = 0.;

  gdouble etaR = (GFS_VALUE (face->neighbor, r->v1[0]) -
		  s->du*GFS_VALUE (face->neighbor, r->dv[face->d/2][0]));
  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE: case FTT_FINE_COARSE:
    /* fixme: this is only first-order accurate for fine/coarse */
    if (etaR > DRY) {
      uR[1] = s->du*(GFS_VALUE (face->neighbor, r->v1[s->u]) -
		     s->du*GFS_VALUE (face->neighbor, r->dv[face->d/2][s->u]))/etaR; /* u = uh/h */
      uR[2] = s->dv*(GFS_VALUE (face->neighbor, r->v1[s->v]) -
		     s->du*GFS_VALUE (face->neighbor, r->dv[face->d/2][s->v]))/etaR; /* v = vh/h */
    }
    else
      uR[1] = uR[2] = 0.;
    uR[0] = MAX (0., etaR + zbR - zbLR);
    uR[3] = 0.;
    break;

  default:
    g_assert_not_reached ();
  }

  riemann_hllc (uL, uR, r->g, f);

  gdouble dt = gfs_domain_face_fraction (r->v[0]->domain, face)*r->dt/ftt_cell_size (face->cell);
  f[0] *= dt;
  f[2] = s->dv*dt*f[2];
  GFS_VALUE (face->cell, r->flux[0])    -= f[0];
  GFS_VALUE (face->cell, r->flux[s->u]) -= s->du*dt*(f[1] - r->g/2.*(uL[0]*uL[0] - etaL*etaL));
  GFS_VALUE (face->cell, r->flux[s->v]) -= f[2];

  f[1] = s->du*dt*(f[1] - r->g/2.*(uR[0]*uR[0] - etaR*etaR));
  if (ftt_face_type (face) == FTT_FINE_COARSE) {
    f[0] /= FTT_CELLS;
    f[1] /= FTT_CELLS;
    f[2] /= FTT_CELLS;
  }
  GFS_VALUE (face->neighbor, r->flux[0])    += f[0];
  GFS_VALUE (face->neighbor, r->flux[s->u]) += f[1];
  GFS_VALUE (face->neighbor, r->flux[s->v]) += f[2];
}

static void reset_fluxes (FttCell * cell, const GfsRiver * r)
{
  guint v;
  for (v = 0; v < GFS_RIVER_NVAR; v++)
    GFS_VALUE (cell, r->flux[v]) = 0.;
}

static void sources (FttCell * cell, GfsRiver * r)
{
  gdouble delta = ftt_cell_size (cell);

  gdouble etaL = GFS_VALUE (cell, r->v1[0]) - GFS_VALUE (cell, r->dv[0][0]);
  gdouble zbL = GFS_VALUE (cell, r->v[3]) - GFS_VALUE (cell, r->dv[0][3]);
  gdouble etaR = GFS_VALUE (cell, r->v1[0]) + GFS_VALUE (cell, r->dv[0][0]);
  gdouble zbR = GFS_VALUE (cell, r->v[3]) + GFS_VALUE (cell, r->dv[0][3]);

  GFS_VALUE (cell, r->vc[1]) += r->dt*r->g/2.*(etaL + etaR)*(zbL - zbR)/delta;

  etaL = GFS_VALUE (cell, r->v1[0]) - GFS_VALUE (cell, r->dv[1][0]);
  zbL = GFS_VALUE (cell, r->v[3]) - GFS_VALUE (cell, r->dv[1][3]);
  etaR = GFS_VALUE (cell, r->v1[0]) + GFS_VALUE (cell, r->dv[1][0]);
  zbR = GFS_VALUE (cell, r->v[3]) + GFS_VALUE (cell, r->dv[1][3]);

  GFS_VALUE (cell, r->vc[2]) += r->dt*r->g/2.*(etaL + etaR)*(zbL - zbR)/delta;
}

static void advance (GfsRiver * r, GfsVariable ** v, gdouble dt)
{
  GfsDomain * domain = GFS_DOMAIN (r);
  guint i;

  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) reset_fluxes, r);
  r->dt = dt;
  gfs_domain_face_traverse (domain, FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) face_fluxes, r);
  r->vc = v;
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) sources, r);
  for (i = 0; i < GFS_RIVER_NVAR; i++) {
    GfsAdvectionParams par;
    par.v = v[i];
    par.fv = r->flux[i];
    par.average = FALSE;
    gfs_domain_traverse_merged (domain, (GfsMergedTraverseFunc) gfs_advection_update, &par);
    gfs_domain_copy_bc (domain, FTT_TRAVERSE_LEAFS, -1, r->v[i], par.v);
  }
}

static void copy (FttCell * cell, const GfsRiver * r)
{
  guint v;
  for (v = 0; v < GFS_RIVER_NVAR; v++)
    GFS_VALUE (cell, r->v1[v]) = GFS_VALUE (cell, r->v[v]);
}

static void cell_H (FttCell * cell,
		    const GfsRiver * r)
{
  GFS_VALUE (cell, r->H) = GFS_VALUE (cell, r->zb) + GFS_VALUE (cell, r->v[0]);
}

static void cell_gradients (FttCell * cell,
			    const GfsRiver * r)
{
  FttComponent c;
  guint v;

  for (c = 0; c < FTT_DIMENSION; c++) {
    for (v = 0; v < GFS_RIVER_NVAR; v++)
      GFS_VALUE (cell, r->dv[c][v]) = (* r->gradient) (cell, c, r->v[v]->i)/2.;
    /* recontruct Zb + eta rather than Zb: see Theorem 3.1 of Audusse et al, 2004 */
    GFS_VALUE (cell, r->dv[c][3]) = 
      (* r->gradient) (cell, c, r->H->i)/2. - GFS_VALUE (cell, r->dv[c][0]);
  }
}

typedef struct { 
  FttCellTraverseFunc func;
  FttDirection d;
  gpointer data;
} FaceTraverseData;

static void face_traverse (FttCell * cell, FaceTraverseData * p)
{
  FttCell * neighbor = ftt_cell_neighbor (cell, p->d);
  if (neighbor)
    (* p->func) (neighbor, p->data);
}

static void domain_traverse_all_leaves (GfsDomain * domain,
					FttCellTraverseFunc func,
					gpointer data)
{
  FaceTraverseData p;

  gfs_domain_traverse_leaves (domain, func, data);
  /* now traverses boundary cells */
  p.func = func;
  p.data = data;
  for (p.d = 0; p.d < FTT_NEIGHBORS; p.d++)
    gfs_domain_cell_traverse_boundary (domain, p.d, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				       (FttCellTraverseFunc) face_traverse, &p);
}

static void river_run (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsRiver * r = GFS_RIVER (sim);

  r->g = sim->physical_params.g/sim->physical_params.L;
  r->gradient = sim->advection_params.gradient;

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  gfs_simulation_set_timestep (sim);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    /* update H */
    domain_traverse_all_leaves (domain, (FttCellTraverseFunc) cell_H, r);
    
    /* events */
    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    /* gradients */
    gfs_domain_timer_start (domain, "gradients");
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) cell_gradients, r);
    FttComponent c;
    guint v;
    for (c = 0; c < FTT_DIMENSION; c++)
      for (v = 0; v < GFS_RIVER_NVAR + 1; v++)
	gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, r->dv[c][v]);
    gfs_domain_timer_stop (domain, "gradients");

    /* predictor */
    gfs_domain_timer_start (domain, "predictor");
    domain_traverse_all_leaves (domain, (FttCellTraverseFunc) copy, r);
    advance (r, r->v1, sim->advection_params.dt/2.);
    gfs_domain_timer_stop (domain, "predictor");

    /* corrector */
    gfs_domain_timer_start (domain, "corrector");
    advance (r, r->v, sim->advection_params.dt);
    gfs_domain_timer_stop (domain, "corrector");

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

static void minimum_cfl (FttCell * cell, GfsRiver * r)
{
  gdouble size = ftt_cell_size (cell);
  gdouble h = GFS_VALUE (cell, r->v[0]);
  if (h > DRY) {
    gdouble uh = fabs (GFS_VALUE (cell, r->v[1]));
    gdouble c = sqrt (r->g*h);
    gdouble cfl = size/(uh/h + c);
    if (cfl < r->cfl)
      r->cfl = cfl;

    gdouble vh = fabs (GFS_VALUE (cell, r->v[2]));
    cfl = size/(vh/h + c);
    if (cfl < r->cfl)
      r->cfl = cfl;
  }
}

static gdouble river_cfl (GfsSimulation * sim)
{
  GfsRiver * r = GFS_RIVER (sim);
  r->cfl = G_MAXDOUBLE;
  gfs_domain_traverse_leaves (GFS_DOMAIN (sim), (FttCellTraverseFunc) minimum_cfl, r);
  gfs_all_reduce (domain, r->cfl, MPI_DOUBLE, MPI_MIN);
  return r->cfl;
}

static void river_class_init (GfsSimulationClass * klass)
{
  klass->run = river_run;
  klass->cfl = river_cfl;
}

static void river_init (GfsRiver * r)
{
  GfsDomain * domain = GFS_DOMAIN (r);

  r->v[0] = gfs_variable_from_name (domain->variables, "P");
  r->v[0]->units = 1.;
  r->v[1] = gfs_variable_from_name (domain->variables, "U");
  r->v[1]->units = 2.;
  r->v[2] = gfs_variable_from_name (domain->variables, "V");
  r->v[2]->units = 2.;

  r->zb = gfs_domain_add_variable (domain, "Zb", "Bed elevation above datum");
  r->zb->units = 1.;
  r->v[3] = r->zb;

  r->H = gfs_domain_add_variable (domain, "H", "Elevation above datum");
  r->H->units = 1.;

  r->flux[0] = gfs_domain_add_variable (domain, "fP", NULL);
  r->flux[1] = gfs_domain_add_variable (domain, "fU", NULL);
  r->flux[2] = gfs_domain_add_variable (domain, "fV", NULL);

  r->v1[0] = gfs_domain_add_variable (domain, "P1", NULL);
  r->v1[1] = gfs_domain_add_variable (domain, "U1", NULL);
  r->v1[1]->component = FTT_X;
  r->v1[2] = gfs_domain_add_variable (domain, "V1", NULL);
  r->v1[2]->component = FTT_Y;

  r->dv[0][0] = gfs_domain_add_variable (domain, "Px", NULL);
  r->dv[1][0] = gfs_domain_add_variable (domain, "Py", NULL);
  r->dv[0][1] = gfs_domain_add_variable (domain, "Ux", NULL);
  r->dv[1][1] = gfs_domain_add_variable (domain, "Uy", NULL);
  r->dv[0][2] = gfs_domain_add_variable (domain, "Vx", NULL);
  r->dv[1][2] = gfs_domain_add_variable (domain, "Vy", NULL);
  r->dv[0][3] = gfs_domain_add_variable (domain, "Zbx", NULL);
  r->dv[1][3] = gfs_domain_add_variable (domain, "Zby", NULL);

  GFS_SIMULATION (r)->advection_params.gradient = gfs_center_minmod_gradient;
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
  gdouble hi = GFS_VALUE (f->neighbor, river->v[0]);

  GFS_VALUE (f->cell, b->v) = GFS_VALUE (f->neighbor, b->v) + 
    (FTT_FACE_DIRECT (f) ? -1. : 1.)*2.*hi*(sqrt (river->g*hi) - sqrt (river->g*hb));
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
