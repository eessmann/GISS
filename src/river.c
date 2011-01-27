/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2009 National Institute of Water and Atmospheric Research
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.  
 */

#include <stdlib.h>
#include "river.h"
#include "adaptive.h"
#include "source.h"
#include "init.h"

/* GfsRiver: Object */

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

#define CFL_CLAMP(u, umax) (fabs (u) <= (umax) ? (u) : (u) > 0. ? (umax) : - (umax))

static void face_fluxes (FttCellFace * face, GfsRiver * r)
{
  if (GFS_VALUE (face->cell, r->v1[0]) <= r->dry &&
      GFS_VALUE (face->neighbor, r->v1[0]) <= r->dry)
    return;

  static Sym sym[4] = {
    {U,  1., V,  1.},
    {U, -1., V, -1.},
    {V,  1., U, -1.},
    {V, -1., U,  1.}
  };
  Sym * s = &sym[face->d];
  gdouble etaL = (GFS_VALUE (face->cell, r->v1[0]) < r->dry ? 0. :
		  GFS_VALUE (face->cell, r->v1[0]) + 
		  s->du*GFS_VALUE (face->cell, r->dv[face->d/2][0]));
  gdouble zbL = (GFS_VALUE (face->cell, r->v[3]) + 
		 s->du*GFS_VALUE (face->cell, r->dv[face->d/2][3]));
  gdouble zbR = (GFS_VALUE (face->neighbor, r->v[3]) -
		 s->du*GFS_VALUE (face->neighbor, r->dv[face->d/2][3]));
  gdouble zbLR = MAX (zbL, zbR);
  gdouble uL[4], uR[4], f[3];

  if (etaL > r->dry) {
    uL[1] = s->du*(GFS_VALUE (face->cell, r->v1[s->u]) +
		   s->du*GFS_VALUE (face->cell, r->dv[face->d/2][s->u]))/etaL; /* u = uh/h */
    uL[2] = s->dv*(GFS_VALUE (face->cell, r->v1[s->v]) +
		   s->du*GFS_VALUE (face->cell, r->dv[face->d/2][s->v]))/etaL; /* v = vh/h */
  }
  else
    uL[1] = uL[2] = 0.;
  uL[0] = MAX (0., etaL + zbL - zbLR);
  uL[3] = 0.;

  gdouble etaR = (GFS_VALUE (face->neighbor, r->v1[0]) < r->dry ? 0. :
		  GFS_VALUE (face->neighbor, r->v1[0]) -
		  s->du*GFS_VALUE (face->neighbor, r->dv[face->d/2][0]));
  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE: case FTT_FINE_COARSE:
    /* fixme: this is only first-order accurate for fine/coarse */
    if (etaR > r->dry) {
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

  gdouble h = ftt_cell_size (face->cell);
  gdouble umax = GFS_SIMULATION (r)->advection_params.cfl*h/r->dt;
  uL[1] = CFL_CLAMP (uL[1], umax);
  uR[1] = CFL_CLAMP (uR[1], umax);
  uL[2] = CFL_CLAMP (uL[2], umax);
  uR[2] = CFL_CLAMP (uR[2], umax);

  riemann_hllc (uL, uR, r->g, f);

  gdouble dt = gfs_domain_face_fraction (GFS_DOMAIN (r), face)*r->dt/h;
  f[0] *= dt;
  f[2] = s->dv*dt*f[2];
  GFS_VALUE (face->cell, r->flux[0])    -= f[0];
  /* see equation 2.16 of Audusse et al, 2004 */
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
  /* metric coefficients */
  gdouble fm[FTT_NEIGHBORS], cm;

  /* Geometric source terms (see doc/figures/lonlat.tm) */
  if (GFS_DOMAIN (r)->cell_metric) {
    GfsDomain * domain = GFS_DOMAIN (r);
    FttCellFace face = { cell };
    for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++)
      fm[face.d] = (* domain->face_metric) (domain, &face);
    gdouble dh_dl = fm[FTT_RIGHT] - fm[FTT_LEFT];
    gdouble dh_dt = fm[FTT_TOP]   - fm[FTT_BOTTOM];
    cm = (* domain->cell_metric) (domain, cell)*ftt_cell_size (cell);
    gdouble dldh = cm*GFS_SIMULATION (r)->physical_params.L;
    gdouble 
      phiu = GFS_VALUE (cell, r->v1[1]), 
      phiv = GFS_VALUE (cell, r->v1[2]);
    gdouble fG = phiv*dh_dl - phiu*dh_dt;
    gdouble g = GFS_SIMULATION (r)->physical_params.g;

    gdouble etaL = GFS_VALUE (cell, r->v1[0]) - GFS_VALUE (cell, r->dv[0][0]);
    gdouble etaR = GFS_VALUE (cell, r->v1[0]) + GFS_VALUE (cell, r->dv[0][0]);
    GFS_VALUE (cell, r->v[1]) += r->dt*(g*(etaL*etaL + etaR*etaR)/4.*dh_dl + fG*phiv)/dldh;

    etaL = GFS_VALUE (cell, r->v1[0]) - GFS_VALUE (cell, r->dv[1][0]);
    etaR = GFS_VALUE (cell, r->v1[0]) + GFS_VALUE (cell, r->dv[1][0]);
    GFS_VALUE (cell, r->v[2]) += r->dt*(g*(etaL*etaL + etaR*etaR)/4.*dh_dt - fG*phiu)/dldh;
  }
  else { /* metric unity */
    FttDirection d;
    for (d = 0; d < FTT_NEIGHBORS; d++)
      fm[d] = 1.;
    cm = ftt_cell_size (cell);
  }

  /* Second-order correction for slope source term ("Sci" of Audusse
     et al, 2004, SIAM, 25(6):2050-2065, equation 3.8) */
  gdouble etaL = GFS_VALUE (cell, r->v1[0]) - GFS_VALUE (cell, r->dv[0][0]);
  gdouble zbL = GFS_VALUE (cell, r->v[3]) - GFS_VALUE (cell, r->dv[0][3]);
  gdouble etaR = GFS_VALUE (cell, r->v1[0]) + GFS_VALUE (cell, r->dv[0][0]);
  gdouble zbR = GFS_VALUE (cell, r->v[3]) + GFS_VALUE (cell, r->dv[0][3]);

  GFS_VALUE (cell, r->v[1]) += 
    r->dt*r->g/4.*(fm[FTT_RIGHT] + fm[FTT_LEFT])*(etaL + etaR)*(zbL - zbR)/cm;

  etaL = GFS_VALUE (cell, r->v1[0]) - GFS_VALUE (cell, r->dv[1][0]);
  zbL = GFS_VALUE (cell, r->v[3]) - GFS_VALUE (cell, r->dv[1][3]);
  etaR = GFS_VALUE (cell, r->v1[0]) + GFS_VALUE (cell, r->dv[1][0]);
  zbR = GFS_VALUE (cell, r->v[3]) + GFS_VALUE (cell, r->dv[1][3]);

  GFS_VALUE (cell, r->v[2]) += 
    r->dt*r->g/4.*(fm[FTT_TOP] + fm[FTT_BOTTOM])*(etaL + etaR)*(zbL - zbR)/cm;
}

static void advance (GfsRiver * r, gdouble dt)
{
  GfsDomain * domain = GFS_DOMAIN (r);
  guint i;

  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) reset_fluxes, r);
  r->dt = dt;
  gfs_domain_face_traverse (domain, FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) face_fluxes, r);
  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) sources, r);
  for (i = 0; i < GFS_RIVER_NVAR; i++) {
    GfsAdvectionParams par;
    par.v = r->v[i];
    par.fv = r->flux[i];
    par.average = FALSE;
    gfs_domain_traverse_merged (domain, (GfsMergedTraverseFunc) gfs_advection_update, &par);
    gfs_domain_variable_centered_sources (domain, par.v, par.v, dt);
  }
  gfs_source_coriolis_implicit (domain, dt);
  for (i = 0; i < GFS_RIVER_NVAR; i++)
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, r->v[i]);
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

  if (GFS_VALUE (cell, r->v[0]) < r->dry) {
    for (c = 0; c < FTT_DIMENSION; c++) {
      for (v = 0; v < GFS_RIVER_NVAR; v++)
	GFS_VALUE (cell, r->dv[c][v]) = 0.;
      GFS_VALUE (cell, r->dv[c][3]) = 0.;
    }
  }
  else { /* wet */
    for (c = 0; c < FTT_DIMENSION; c++) {
      for (v = 0; v < GFS_RIVER_NVAR; v++)
	GFS_VALUE (cell, r->dv[c][v]) = (* r->gradient) (cell, c, r->v[v]->i)/2.;
      /* recontruct Zb + eta rather than Zb: see Theorem 3.1 of Audusse et al, 2004 */
      GFS_VALUE (cell, r->dv[c][3]) =
	(* r->gradient) (cell, c, r->H->i)/2.
	- GFS_VALUE (cell, r->dv[c][0]);
    }
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

  r->v[3] = r->zb = gfs_variable_from_name (domain->variables, "Zb");

  r->g = sim->physical_params.g/sim->physical_params.L;
  r->gradient = sim->advection_params.gradient;

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  gfs_simulation_set_timestep (sim);

  domain_traverse_all_leaves (domain, (FttCellTraverseFunc) cell_H, r);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    /* events */
    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    /* update H */
    domain_traverse_all_leaves (domain, (FttCellTraverseFunc) cell_H, r);
    
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
    domain_traverse_all_leaves (domain, (FttCellTraverseFunc) copy, r);
    if (r->time_order == 2) {
      gfs_domain_timer_start (domain, "predictor");
      for (v = 0; v < GFS_RIVER_NVAR; v++)
	gfs_variables_swap (r->v[v], r->v1[v]);
      advance (r, sim->advection_params.dt/2.);
      for (v = 0; v < GFS_RIVER_NVAR; v++)
	gfs_variables_swap (r->v[v], r->v1[v]);
      gfs_domain_timer_stop (domain, "predictor");
    }
    /* corrector */
    gfs_domain_timer_start (domain, "corrector");
    advance (r, sim->advection_params.dt);
    gfs_domain_timer_stop (domain, "corrector");

    /* update H */
    domain_traverse_all_leaves (domain, (FttCellTraverseFunc) cell_H, r);

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

static gdouble maximum_face_metric (FttCell * cell, GfsDomain * domain, FttComponent c)
{
  if (domain->face_metric) {
    FttCellFace f;
    f.cell = cell; f.d = 2*c;
    gdouble fm1 = (* domain->face_metric) (domain, &f);
    f.d = 2*c + 1;
    gdouble fm2 = (* domain->face_metric) (domain, &f);
    return MAX (fm1, fm2);
  }
  else
    return 1.;
}

static void minimum_cfl (FttCell * cell, GfsRiver * r)
{
  gdouble h = GFS_VALUE (cell, r->v[0]);
  if (h > r->dry) {
    GfsDomain * domain = GFS_DOMAIN (r);
    gdouble vol = ftt_cell_size (cell);
    if (domain->cell_metric)
      vol *= (* domain->cell_metric) (domain, cell);
    gdouble cg = sqrt (r->g*h);
    FttComponent c;
    for (c = FTT_X; c <= FTT_Y; c++) {
      gdouble uh = fabs (GFS_VALUE (cell, r->v[c + 1]));
      gdouble fm = maximum_face_metric (cell, domain, c);
      gdouble cfl = vol/(fm*(uh/h + cg));
      if (cfl < r->cfl)
	r->cfl = cfl;
    }
  }
}

static gdouble river_cfl (GfsSimulation * sim)
{
  GfsRiver * r = GFS_RIVER (sim);
  r->cfl = G_MAXDOUBLE;
  gfs_domain_traverse_leaves (GFS_DOMAIN (sim), (FttCellTraverseFunc) minimum_cfl, r);
  gfs_all_reduce (GFS_DOMAIN (sim), r->cfl, MPI_DOUBLE, MPI_MIN);
  return r->cfl;
}

static void river_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_river_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsRiver * river = GFS_RIVER (*o);
  if (fp->type == '{') {
    double dry;
    GtsFileVariable var[] = {
      {GTS_UINT,   "time_order", TRUE, &river->time_order},
      {GTS_DOUBLE, "dry",        TRUE, &dry},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
    if (var[1].set)
      river->dry = dry/GFS_SIMULATION (river)->physical_params.L;
  }
}

static void river_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_river_class ())->parent_class->write) (o, fp);

  GfsRiver * river = GFS_RIVER (o);
  fprintf (fp, " {\n"
	   "  time_order = %d\n"
	   "  dry = %g\n"
	   "}",
	   river->time_order,
	   river->dry*GFS_SIMULATION (river)->physical_params.L);
}

static void river_class_init (GfsSimulationClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = river_read;
  GTS_OBJECT_CLASS (klass)->write = river_write;
  klass->run = river_run;
  klass->cfl = river_cfl;
}

static gdouble cell_velocity (FttCell * cell, FttCellFace * face, GfsRiver * r)
{
  g_return_val_if_fail (cell != NULL, 0.);

  gdouble D = GFS_VALUE (cell, r->v[0]);
  gdouble L = GFS_SIMULATION (r)->physical_params.L;
  return D > r->dry ? L*gfs_vector_norm (cell, gfs_domain_velocity (GFS_DOMAIN (r)))/D : 0.;
}

static gdouble cell_velocity2 (FttCell * cell, FttCellFace * face, GfsRiver * r)
{
  g_return_val_if_fail (cell != NULL, 0.);

  gdouble D = GFS_VALUE (cell, r->v[0]);
  gdouble L = GFS_SIMULATION (r)->physical_params.L;
  return D > r->dry ? 
    L*L*gfs_vector_norm2 (cell, gfs_domain_velocity (GFS_DOMAIN (r)))/(D*D) : 0.;
}

static void momentum_coarse_fine (FttCell * parent, GfsVariable * v)
{
  /* Only initializes momentum when @parent is "deep enough". This
     assumes that shallow parents have just been submerged. For these
     shallow cells, childrens' initial momentum defaults to zero. This
     prevents creating spurious large velocities. */
  GfsRiver * r = GFS_RIVER (v->domain);
  if (GFS_VALUE (parent, r->v[0]) > 2.*r->dry)
    gfs_cell_coarse_fine (parent, v);
}

static void river_init (GfsRiver * r)
{
  GfsDomain * domain = GFS_DOMAIN (r);

  gts_object_destroy (GTS_OBJECT (gfs_variable_from_name (domain->variables, "Pmac")));

  r->v[0] = gfs_variable_from_name (domain->variables, "P");
  r->v[0]->units = 1.;
  g_free (r->v[0]->description);
  r->v[0]->description = g_strdup ("Fluid depth");

  r->v[1] = gfs_variable_from_name (domain->variables, "U");
  r->v[1]->units = 2.;
  g_free (r->v[1]->description);
  r->v[1]->description = g_strdup ("x-component of the fluid flux");
  r->v[1]->coarse_fine = momentum_coarse_fine;

  r->v[2] = gfs_variable_from_name (domain->variables, "V");
  r->v[2]->units = 2.;
  g_free (r->v[2]->description);
  r->v[2]->description = g_strdup ("y-component of the fluid flux");
  r->v[2]->coarse_fine = momentum_coarse_fine;

  r->zb = gfs_domain_add_variable (domain, "Zb", "Bed elevation above datum");
  r->zb->units = 1.;

  r->H = gfs_domain_add_variable (domain, "H", "Elevation above datum (Zb + P)");
  r->H->units = 1.;

  r->flux[0] = gfs_domain_add_variable (domain, NULL, NULL);
  r->flux[1] = gfs_domain_add_variable (domain, NULL, NULL);
  r->flux[2] = gfs_domain_add_variable (domain, NULL, NULL);

  r->v1[0] = gfs_domain_add_variable (domain, NULL, NULL);
  r->v1[1] = gfs_domain_add_variable (domain, NULL, NULL);
  r->v1[2] = gfs_domain_add_variable (domain, NULL, NULL);
  gfs_variable_set_vector (&r->v1[1], 2);

  r->dv[0][0] = gfs_domain_add_variable (domain, "Px", "x-component of the depth gradient");
  r->dv[1][0] = gfs_domain_add_variable (domain, "Py", "y-component of the depth gradient");
  r->dv[0][1] = gfs_domain_add_variable (domain, "Ux", "x-component of the flux gradient");
  r->dv[1][1] = gfs_domain_add_variable (domain, "Uy", "y-component of the flux gradient");
  r->dv[0][2] = gfs_domain_add_variable (domain, "Vx", "x-component of the flux gradient");
  r->dv[1][2] = gfs_domain_add_variable (domain, "Vy", "y-component of the flux gradient");
  r->dv[0][3] = gfs_domain_add_variable (domain, "Zbx", "x-component of the bed slope");
  r->dv[1][3] = gfs_domain_add_variable (domain, "Zby", "y-component of the bed slope");

  GFS_SIMULATION (r)->advection_params.gradient = gfs_center_minmod_gradient;
  GFS_SIMULATION (r)->advection_params.cfl = 0.5;

  GfsDerivedVariable * v = gfs_derived_variable_from_name (domain->derived_variables, "Velocity");
  v->func = cell_velocity;
  v = gfs_derived_variable_from_name (domain->derived_variables, "Velocity2");
  v->func = cell_velocity2;

  gfs_domain_remove_derived_variable (domain, "Vorticity");
  gfs_domain_remove_derived_variable (domain, "Divergence");
  gfs_domain_remove_derived_variable (domain, "Lambda2");
  gfs_domain_remove_derived_variable (domain, "Curvature");
  gfs_domain_remove_derived_variable (domain, "D2");

  r->time_order = 2;
  r->dry = 1e-6;
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

  g_assert (hi >= 0.);
  GFS_VALUE (f->cell, b->v) = GFS_VALUE (f->neighbor, b->v) + 
    (FTT_FACE_DIRECT (f) ? -1. : 1.)*2.*hi*(sqrt (river->g*hi) - sqrt (river->g*MAX (hb, 0.)));
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

/* GfsDischargeElevation: Object */

static void discharge_elevation_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_DISCHARGE_ELEVATION (o)->Q));
  gts_object_destroy (GTS_OBJECT (GFS_DISCHARGE_ELEVATION (o)->profile));

  (* GTS_OBJECT_CLASS (gfs_discharge_elevation_class ())->parent_class->destroy) (o);
}

static void discharge_elevation_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_discharge_elevation_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!GFS_IS_RIVER (domain)) {
    gts_file_error (fp, "GfsDischargeElevation only makes sense for GfsRiver simulations");
    return;
  }
  GfsDischargeElevation * bd = GFS_DISCHARGE_ELEVATION (*o);
  gfs_function_read (bd->Q, domain, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '\n')
    gfs_function_read (bd->profile, domain, fp);
  else
    gfs_object_simulation_set (bd->profile, domain);

  bd->P = GFS_RIVER (domain)->v[0];
  g_free (GFS_CONSTANT (bd)->derived->description);
  GFS_CONSTANT (bd)->derived->description = g_strdup ("Elevation for a given discharge");
}

static void discharge_elevation_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_discharge_elevation_class ())->parent_class->write) (o, fp);
  gfs_function_write (GFS_DISCHARGE_ELEVATION (o)->Q, fp);
  if (gfs_function_get_constant_value (GFS_DISCHARGE_ELEVATION (o)->profile) != 0.)
    gfs_function_write (GFS_DISCHARGE_ELEVATION (o)->profile, fp);
}

static void boundary_flux (FttCellFace * f, GfsDischargeElevation * b)
{
  gdouble profile = gfs_function_face_value (b->profile, f);
  if (profile != GFS_NODATA) {
    GfsRiver * river = GFS_RIVER (gfs_object_simulation (b));
    GFS_VALUE (f->cell, river->flux[0]) = 0.;
    gdouble v1 = GFS_VALUE (f->cell, river->v1[0]);
    GFS_VALUE (f->cell, river->v1[0]) = MAX (0.,
					     profile + GFS_CONSTANT (b)->val - 
					     gfs_face_interpolated_value (f, river->zb->i));
    gdouble dt = river->dt;
    river->dt = 1.;
    face_fluxes (f, river);
    river->dt = dt;
    GFS_VALUE (f->cell, river->v1[0]) = v1;
    double h = ftt_cell_size (f->cell);
    b->flow -= GFS_VALUE (f->cell, river->flux[0])*h*h;
  }
}

static void traverse_dirichlet_boundaries (GfsBox * box, GfsDischargeElevation * bd)
{
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, bd->P);
      if (GFS_IS_BC_DIRICHLET (bc))
	ftt_face_traverse_boundary (b->root, b->d,
				    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				    (FttFaceTraverseFunc) boundary_flux, bd);
    }
}

static gboolean discharge_elevation_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS 
			  (gfs_discharge_elevation_class ())->parent_class)->event)
      (event, sim)) {
    GfsConstant * c = GFS_CONSTANT (event);
    GfsDischargeElevation * bd = GFS_DISCHARGE_ELEVATION (event);
    GfsRiver * r = GFS_RIVER (sim);
    guint v;
      
    for (v = 0; v < GFS_RIVER_NVAR; v++)
      gfs_variables_swap (r->v[v], r->v1[v]);
      
    gfs_catch_floating_point_exceptions ();
    gdouble Q = gfs_function_value (bd->Q, NULL);
    gfs_restore_fpe_for_function (bd->Q);
    gdouble hmax, hmin = 0.;
    gdouble L = sim->physical_params.L;
      
    hmax = c->val*2./L;
    bd->flow = 0.;
    c->val = hmax;
    gts_container_foreach (GTS_CONTAINER (sim), (GtsFunc) traverse_dirichlet_boundaries, bd);
    gfs_all_reduce (GFS_DOMAIN (sim), bd->flow, MPI_DOUBLE, MPI_SUM);
    if (Q > bd->flow)
      hmax = 1.;
      
    guint n = 0, nitermin = 4, nitermax = 100;
    c->val = hmax/2.;
    do {
      bd->flow = 0.;
      gts_container_foreach (GTS_CONTAINER (sim), (GtsFunc) traverse_dirichlet_boundaries, bd);
      gfs_all_reduce (GFS_DOMAIN (sim), bd->flow, MPI_DOUBLE, MPI_SUM);
      if (bd->flow > Q) {
	hmax = c->val;
	c->val = (c->val + hmin)/2.;
      }
      else {
	hmin = c->val;
	c->val = (c->val + hmax)/2.;
      }
      n++;
    } while (n < nitermax && (n < nitermin || fabs (Q - bd->flow)/Q > bd->tolerance));
    if (n == nitermax)
      g_warning ("discharge_elevation_event() did not converge after %d iterations: %g", n, 
		 fabs (Q - bd->flow)/Q);
    c->val *= L;
    /* g_message ("### flow: %g H: %g nitermax: %d\n", bd->flow*pow(L,3), c->val, n); */

    for (v = 0; v < GFS_RIVER_NVAR; v++)
      gfs_variables_swap (r->v[v], r->v1[v]);

    return TRUE;
  }
  return FALSE;
}

static void gfs_discharge_elevation_class_init (GtsObjectClass * klass)
{
  klass->destroy  = discharge_elevation_destroy;
  klass->read     = discharge_elevation_read;
  klass->write    = discharge_elevation_write;
  GFS_EVENT_CLASS (klass)->event = discharge_elevation_event;
}

static void gfs_discharge_elevation_init (GfsDischargeElevation * b)
{
  b->tolerance = 1e-2;
  b->Q = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_set_units (b->Q, 3.);
  b->profile = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_set_units (b->profile, 1.);
}

GfsEventClass * gfs_discharge_elevation_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_discharge_elevation_info = {
      "GfsDischargeElevation",
      sizeof (GfsDischargeElevation),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_discharge_elevation_class_init,
      (GtsObjectInitFunc) gfs_discharge_elevation_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_constant_class ()),
				  &gfs_discharge_elevation_info);
  }

  return klass;
}
