/* Gerris - The GNU Flow Solver
 * Copyright (C) 2008-2012 National Institute of Water and Atmospheric Research
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
/*! \file
 * \brief GfsRiver model.
 */

#include <stdlib.h>
#include "river.h"
#include "adaptive.h"
#include "source.h"
#include "solid.h"
#include "init.h"

/* generalisation of the limited gradients (in fluid.c) to mixed cells */

static gdouble generic_limiter (gdouble r, gdouble beta)
{
  gdouble v1 = MIN (r, beta), v2 = MIN (beta*r, 1.);
  v1 = MAX (0., v1);
  return MAX (v1, v2);
}

static gdouble minmod_limiter (gdouble r)
{
  return generic_limiter (r, 1.);
}

static gdouble superbee_limiter (gdouble r)
{
  return generic_limiter (r, 2.);
}

static gdouble sweby_limiter (gdouble r)
{
  return generic_limiter (r, 1.5);
}

static gdouble center_limited_gradient_full (FttCell * cell,
					     FttComponent c,
					     guint v,
					     gdouble (* limiter) (gdouble))
{
  FttDirection d = 2*c;
  FttCellFace f1;
  gdouble v0;

  f1 = gfs_cell_face (cell, FTT_OPPOSITE_DIRECTION (d));
  v0 = GFS_VALUEI (cell, v);
  if (f1.neighbor) {
    FttCellFace f2 = gfs_cell_face (cell, d);
    if (f2.neighbor) {
      /* two neighbors */
      gdouble x1 = 1., v1, x2 = 1., v2;
      v1 = gfs_neighbor_value (&f1, v, &x1);
      v2 = gfs_neighbor_value (&f2, v, &x2);

      gdouble g;
      if (v0 == v1)
	g = 0.;
      else
	g = (* limiter) ((v2 - v0)*x1/((v0 - v1)*x2))*(v0 - v1)/x1;
      return g;
    }
  }
  /* only one or no neighbors */
  return 0.;
}

static gdouble center_limited_gradient (FttCell * cell,
					FttComponent c,
					guint v,
					gdouble (* limiter) (gdouble))
{
  FttDirection d = 2*c;
  FttCellFace f1, f2;
  f1 = gfs_cell_face (cell, FTT_OPPOSITE_DIRECTION (d));
  f2 = gfs_cell_face (cell, d);
  if (!GFS_IS_MIXED (cell) && 
      (!f1.neighbor || !GFS_IS_MIXED (f1.neighbor)) &&
      (!f2.neighbor || !GFS_IS_MIXED (f2.neighbor)))
    return center_limited_gradient_full (cell, c, v, limiter);
  gdouble h = ftt_cell_size (cell);
  FttVector cm;
  gfs_cell_cm (cell, &cm);  
  gdouble v0 = GFS_VALUEI (cell, v), g = 0.;
  
  if (f1.neighbor && f2.neighbor) {
    /* two neighbors */
    gdouble x1, x2;
    gdouble v1 = gfs_neighbor_value (&f1, v, &x1);
    gdouble v2 = gfs_neighbor_value (&f2, v, &x2);
    if (v0 != v1) {
      FttVector cm1, cm2;
      gfs_cell_cm (f1.neighbor, &cm1);
      gfs_cell_cm (f2.neighbor, &cm2);       
      /* fixme: this is not correct at coarse/fine boundaries */
      x1 = ((&cm.x)[c] - (&cm1.x)[c])/h;
      x2 = ((&cm2.x)[c] - (&cm.x)[c])/h;
      g = (* limiter) ((v2 - v0)*x1/((v0 - v1)*x2))*(v0 - v1)/x1;
    }
  }
 
  /* mixed cells gradient following Causon et al. (2000) */
  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * s = GFS_STATE (cell)->solid;
    FttVector ca = s->ca;
    FttVector n;
    gdouble nn;

    gfs_solid_normal (cell, &n);
    nn = sqrt (n.x*n.x + n.y*n.y);
    n.x /= nn;
    n.y /= nn;
    
    /* solid is on the right side of the cell */
    if (s->s[2*c] < s->s[2*c + 1]) {
      if (f1.neighbor) {
	gdouble vr;
	/* fixme: this relies on specific indices for U and V. Not recommended. */
	if (v == 2) vr = v0 - 2.*(n.x*GFS_VALUEI (cell, 2) + n.y*GFS_VALUEI (cell, 3))*n.x;
	else if (v == 3) vr = v0 - 2.*(n.x*GFS_VALUEI (cell, 2) + n.y*GFS_VALUEI (cell, 3))*n.y;
	else return s->s[2*c]*g/s->s[2*c + 1];
	gdouble x1, v1 = gfs_neighbor_value (&f1, v, &x1);
	FttVector cm1;
	gfs_cell_cm (f1.neighbor, &cm1);
	/* fixme: this is not correct at coarse/fine boundaries */
	x1 = ((&cm.x)[c] - (&cm1.x)[c])/h;
	gdouble x2 = 2.*((&ca.x)[c] - (&cm.x)[c])/h;
	gdouble gs = ((v0 - v1)*x2 == 0. || x1 == 0.) ? 0. :
	  (* limiter) ((vr - v0)*x1/((v0 - v1)*x2))*(v0 - v1)/x1;
	return (s->s[2*c]*g + (s->s[2*c + 1] - s->s[2*c])*gs)/s->s[2*c + 1];
      }
      else 
	return 0;
    }    
    /* solid is on the left side of the cell */
    else if (s->s[2*c] > s->s[2*c + 1]) {
      if (f2.neighbor) {
	gdouble vr;
	/* fixme: this relies on specific indices for U and V. Not recommended. */
	if (v == 2) vr = v0 - 2.*(n.x*GFS_VALUEI (cell, 2) + n.y*GFS_VALUEI (cell, 3))*n.x; 
	else if (v == 3) vr = v0 - 2.*(n.x*GFS_VALUEI (cell, 2) + n.y*GFS_VALUEI (cell, 3))*n.y; 
	else return s->s[2*c + 1]*g/s->s[2*c];
	gdouble x2, v2 = gfs_neighbor_value (&f2, v, &x2);
	FttVector cm2;
	gfs_cell_cm (f2.neighbor, &cm2);
	/* fixme: this is not correct at coarse/fine boundaries */
	gdouble x1 = 2.*((&cm.x)[c] - (&ca.x)[c])/h;
	x2 = ((&cm2.x)[c] - (&cm.x)[c])/h;
 	gdouble gs = ((v0 - vr)*x2 == 0. || x1 == 0.) ? 0. :
	  (* limiter) ((v2 - v0)*x1/((v0 - vr)*x2))*(v0 - vr)/x1;
	return (s->s[2*c + 1]*g + (s->s[2*c] - s->s[2*c + 1])*gs)/s->s[2*c];
      }
      else 
	return 0;
    }
  }
  return g;
}

static gdouble center_minmod_gradient (FttCell * cell,
				       FttComponent c,
				       guint v)
{
  return center_limited_gradient (cell, c, v, minmod_limiter);
}

static gdouble center_superbee_gradient (FttCell * cell,
					 FttComponent c,
					 guint v)
{
  return center_limited_gradient (cell, c, v, superbee_limiter);
}

static gdouble center_sweby_gradient (FttCell * cell,
				      FttComponent c,
				      guint v)
{
  return center_limited_gradient (cell, c, v, sweby_limiter);
}

/**
 * Solves the Saint-Venant equations.
 * \beginobject{GfsRiver}
 */

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

/*
 * uL: left state vector [h,u,v,zb].
 * uR: right state vector.
 * g: acceleration of gravity.
 * f: flux vector.
 *
 * Fills @f by solving an approximate Riemann problem using the kinetic
 * scheme. See Audusse & Bristeau, JCP, 206, 2005.
 */

#define SQRT3 1.73205080756888

static void riemann_kinetic (const gdouble * uL, const gdouble * uR,
			     gdouble g,
			     gdouble * f)
{
  gdouble ci, Mp, Mm, cig;

  ci = sqrt (g*uL[0]/2.);
  Mp = MAX (uL[1] + ci*SQRT3, 0.);
  Mm = MAX (uL[1] - ci*SQRT3, 0.);
  cig = ci/(6.*g*SQRT3);
  f[0] = cig*3.*(Mp*Mp - Mm*Mm);
  f[1] = cig*2.*(Mp*Mp*Mp - Mm*Mm*Mm);

  ci = sqrt (g*uR[0]/2.);
  Mp = MIN (uR[1] + ci*SQRT3, 0.);
  Mm = MIN (uR[1] - ci*SQRT3, 0.);
  cig = ci/(6.*g*SQRT3);
  f[0] += cig*3.*(Mp*Mp - Mm*Mm);
  f[1] += cig*2.*(Mp*Mp*Mp - Mm*Mm*Mm);

  f[2] = (f[0] > 0. ? uL[2] : uR[2])*f[0];
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
  gdouble eta = GFS_VALUE (face->cell, r->v1[0]), etan = GFS_VALUE (face->neighbor, r->v1[0]);

  if (eta <= r->dry && etan <= r->dry)
    return;

  static Sym sym[4] = {
    {U,  1., V,  1.},
    {U, -1., V, -1.},
    {V,  1., U, -1.},
    {V, -1., U,  1.}
  };
  Sym * s = &sym[face->d];

  gdouble a = 1., b = 1.;
  if (GFS_IS_MIXED (face->cell)) {
    FttVector ca, cm;
    gfs_face_ca (face, &ca);
    gfs_cell_cm (face->cell, &cm);
    FttComponent c = face->d/2;
    a = fabs (2.*((&ca.x)[c] - (&cm.x)[c])/ftt_cell_size (face->cell));
  }
  if (GFS_IS_MIXED (face->neighbor)) {
    FttVector ca, cm;
    gfs_face_ca (face, &ca); /* fixme?: this is not symmetric with the above for face->cell */
    gfs_cell_cm (face->neighbor, &cm);
    FttComponent c = face->d/2;
    b = fabs (2.*((&ca.x)[c] - (&cm.x)[c])/ftt_cell_size (face->neighbor));
  } 

  gdouble etaL = (eta <= r->dry ? 0. : 
		  eta + a*s->du*GFS_VALUE (face->cell, r->dv[face->d/2][0]));
  gdouble zbL = (GFS_VALUE (face->cell, r->v[3])
		 + a*s->du*GFS_VALUE (face->cell, r->dv[face->d/2][3]));
  gdouble zbR = (GFS_VALUE (face->neighbor, r->v[3])
		 - b*s->du*GFS_VALUE (face->neighbor, r->dv[face->d/2][3])); 
  gdouble zbLR = MAX (zbL, zbR);
  gdouble uL[4], uR[4], f[3];

  if (etaL > r->dry) {
    uL[1] = s->du*(GFS_VALUE (face->cell, r->v1[s->u]) +
		   a*s->du*GFS_VALUE (face->cell, r->dv[face->d/2][s->u]))/etaL; /* u = uh/h */
    uL[2] = s->dv*(GFS_VALUE (face->cell, r->v1[s->v]) +
		   a*s->du*GFS_VALUE (face->cell, r->dv[face->d/2][s->v]))/etaL; /* v = vh/h */
  }
  else
    uL[1] = uL[2] = 0.;
  uL[0] = MAX (0., etaL + zbL - zbLR);
  uL[3] = 0.;

  gdouble etaR = (etan <= r->dry ? 0. :
		  etan - b*s->du*GFS_VALUE (face->neighbor, r->dv[face->d/2][0]));
  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE: case FTT_FINE_COARSE:
    /* fixme: this is only first-order accurate for fine/coarse */
    if (etaR > r->dry) {
      /* u = uh/h */
      uR[1] = s->du*(GFS_VALUE (face->neighbor, r->v1[s->u]) -
		     b*s->du*GFS_VALUE (face->neighbor, r->dv[face->d/2][s->u]))/etaR; 
      /* v = vh/h */
      uR[2] = s->dv*(GFS_VALUE (face->neighbor, r->v1[s->v]) -
		     b*s->du*GFS_VALUE (face->neighbor, r->dv[face->d/2][s->v]))/etaR;
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

  (* r->scheme) (uL, uR, r->g, f);

  gdouble dt = gfs_domain_face_fraction (GFS_DOMAIN (r), face)*r->dt/h;
  f[0] *= dt;
  f[2] = s->dv*dt*f[2];
  GFS_VALUE (face->cell, r->flux[0]) -= f[0];
  /* see Audusse and Bristeau, 2005, JCP, 206:311-333 */
  /* Slope source term "S_{i,j,p}" and second-order correction for
     slope source term "Sc_{i,j,p}" of An et al. 2012, Advances in
     Water Resources. 39:60-70, equations (11) and (12) */
  if (eta <= r->dry) eta = 0.;
  gdouble zb = GFS_VALUE (face->cell, r->v[3]);
  GFS_VALUE (face->cell, r->flux[s->u]) -= 
    s->du*dt*(f[1] - r->g/2.*(uL[0]*uL[0] - etaL*etaL - (etaL + eta)*(zbL - zb)));
  GFS_VALUE (face->cell, r->flux[s->v]) -= f[2];

  if (etan <= r->dry) etan = 0.;
  zb = GFS_VALUE (face->neighbor, r->v[3]);
  f[1] = s->du*dt*(f[1] - r->g/2.*(uR[0]*uR[0] - etaR*etaR - (etaR + etan)*(zbR - zb)));
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

static void solid_boundary_fluxes (FttCell * cell, GfsRiver * r)
{
  gdouble h = ftt_cell_size (cell);
  GfsSolidVector * s = GFS_STATE (cell)->solid;
  FttVector cm;

  gdouble hh = MAX (GFS_VALUE (cell, r->v1[0]), 0.);
  gfs_cell_cm (cell, &cm);
  gdouble hs = hh + 2.*((s->ca.x - cm.x)*GFS_VALUE (cell, r->dv[0][0]) + 
			(s->ca.y - cm.y)*GFS_VALUE (cell, r->dv[1][0]))/h;
  gdouble zbs = 2.*((s->ca.x - cm.x)*GFS_VALUE (cell, r->dv[0][3]) +
		    (s->ca.y - cm.y)*GFS_VALUE (cell, r->dv[1][3]))/h;
  gdouble hszbs = 
    r->dt/h*r->g/2.*(
		     /* the normal component of the velocity is zero at the solid boundary. */  
		     hs*hs +
		     /* Second-order correction for slope source term ("Sc_{i,j,s}" of 
			An et al. 2012, Advances in Water Resources. 39:60-70, equation (27) */
		     (hs + hh)*zbs);
  FttVector n;
  gfs_solid_normal (cell, &n);
  GFS_VALUE (cell, r->flux[1]) -= n.x*hszbs;
  GFS_VALUE (cell, r->flux[2]) -= n.y*hszbs;
}

/* Metric source terms (see doc/figures/lonlat.tm) */
static void metric_sources (FttCell * cell, GfsRiver * r)
{
  if (GFS_VALUE (cell, r->v1[0]) > r->dry) {
    /* fixme: this will probably not work when combined with solids */
    GfsDomain * domain = GFS_DOMAIN (r);
    gdouble fm[FTT_NEIGHBORS], cm;
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

    gdouble eta = GFS_VALUE (cell, r->v1[0]);
    GFS_VALUE (cell, r->v[1]) += r->dt*(g*eta*eta/2.*dh_dl + fG*phiv)/dldh;
    GFS_VALUE (cell, r->v[2]) += r->dt*(g*eta*eta/2.*dh_dt - fG*phiu)/dldh;
  }
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
  gfs_domain_traverse_mixed (domain,
  			     FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
  			     (FttCellTraverseFunc) solid_boundary_fluxes, domain);
  if (domain->cell_metric)
    gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) metric_sources, r);
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

  if (GFS_VALUE (cell, r->v[0]) <= r->dry) {
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
  if (r->gradient == gfs_center_minmod_gradient)
    r->gradient = center_minmod_gradient;
  else if (r->gradient == gfs_center_superbee_gradient)
    r->gradient = center_superbee_gradient;
  else if (r->gradient == gfs_center_sweby_gradient)
    r->gradient = center_sweby_gradient;

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
    gchar * scheme = NULL;
    GtsFileVariable var[] = {
      {GTS_UINT,   "time_order", TRUE, &river->time_order},
      {GTS_DOUBLE, "dry",        TRUE, &dry},
      {GTS_STRING, "scheme",     TRUE, &scheme},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
    if (var[1].set)
      river->dry = dry/GFS_SIMULATION (river)->physical_params.L;
    if (scheme) {
      if (!strcmp (scheme, "hllc"))
	river->scheme = riemann_hllc;
      else if (!strcmp (scheme, "kinetic"))
	river->scheme = riemann_kinetic;
      else
	gts_file_error (fp, "unknown scheme '%s'", scheme);
      g_free (scheme);
    }
  }
}

static void river_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_river_class ())->parent_class->write) (o, fp);

  GfsRiver * river = GFS_RIVER (o);
  fprintf (fp, " {\n"
	   "  time_order = %d\n"
	   "  dry = %g\n"
	   "  scheme = %s\n"
	   "}",
	   river->time_order,
	   river->dry*GFS_SIMULATION (river)->physical_params.L,
	   river->scheme == riemann_hllc ? "hllc" : "kinetic");
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
  r->scheme = riemann_hllc;
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

/** \endobject{GfsRiver} */

/**
 * 
 * \beginobject{GfsBcSubcritical}
 */

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

/** \endobject{GfsBcSubcritical} */

/**
 *
 * \beginobject{GfsDischargeElevation}
 */

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
  GFS_EVENT (b)->start = 0.; /* this is not an "Init" event */
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

/** \endobject{GfsDischargeElevation} */

/**
 * "Pipe" between two locations.
 * \beginobject{GfsSourcePipe}
 */

static gboolean read_position (GtsFile * fp, FttVector * p)
{
  gchar * start[FTT_DIMENSION];
  if (!gfs_read_vector (fp, start))
    return FALSE;
  FttComponent c;
  p->z = 0.;
  for (c = 0; c < FTT_DIMENSION; c++) {
    (&p->x)[c] = atof (start[c]);
    g_free (start[c]);
  }
  return TRUE;
}

static void source_pipe_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_pipe_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsSimulation * sim = gfs_object_simulation (*o);
  if (!GFS_IS_RIVER (sim)) {
    gts_file_error (fp, "%s only makes sense for GfsRiver simulations",
		    (*o)->klass->info.name);
    return;
  }
  GfsVariable * v = GFS_RIVER (sim)->v[0];
  if (v->sources == NULL)
    v->sources = gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
  gts_container_add (v->sources, GTS_CONTAINEE (*o));

  GfsSourcePipe * p = GFS_SOURCE_PIPE (*o);
  if (!read_position (fp, &p->start))
    return;
  if (!read_position (fp, &p->end))
    return;

  p->diameter = gfs_read_constant (fp, GFS_DOMAIN (gfs_object_simulation (p)));
  if (p->diameter == G_MAXDOUBLE)
    return;
}

static void source_pipe_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_pipe_class ())->parent_class->write) (o, fp);
  GfsSourcePipe * p = GFS_SOURCE_PIPE (o);
  fprintf (fp, " (%f,%f) (%f,%f) %f",
	   p->start.x, p->start.y,
	   p->end.x, p->end.y,
	   p->diameter);
}

#define DQ (1e-4/L3)

static double flow_rate_Q (double z1, double h1, double z2, double h2,
			   double l, double g, GfsSourcePipe * p,
			   double a1, double a2, double Q)
{
  double Q1 = (*p->flow_rate) (z1, h1 - Q/a1, z2, h2 + Q/a2, l, g, p);
  if (Q1 > 0.) Q1 = MIN (Q1, a1*h1);
  if (Q1 < 0.) Q1 = MAX (Q1, - a2*h2);
  return Q1;
}

static gboolean source_pipe_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* gfs_event_class ()->event) (event, sim)) {
    GfsSourcePipe * p = GFS_SOURCE_PIPE (event);
    GfsDomain * domain = GFS_DOMAIN (sim);
    FttVector start = p->start, end = p->end;
    gfs_simulation_map (sim, &start);
    gfs_simulation_map (sim, &end);
    /* fixme: this won't work in parallel if the ends of the pipe are on different PEs */
    p->scell = gfs_domain_locate (domain, start, -1, NULL);
    p->ecell = gfs_domain_locate (domain, end, -1, NULL);
    p->Q = 0.;
    if (p->scell && p->ecell && p->scell != p->ecell) {
      gdouble L = sim->physical_params.L, g = sim->physical_params.g;
      GfsVariable * h = GFS_RIVER (sim)->v[0], * zb = GFS_RIVER (sim)->zb;
      gdouble h1 = L*GFS_VALUE (p->scell, h), z1 = L*GFS_VALUE (p->scell, zb);
      gdouble h2 = L*GFS_VALUE (p->ecell, h), z2 = L*GFS_VALUE (p->ecell, zb);      
      /* fixme: the length below does not take into account metric
	 properly (e.g. won't work for MetricLonLat) */
      gdouble l = L*sqrt ((start.x - end.x)*(start.x - end.x) +
			  (start.y - end.y)*(start.y - end.y));
      gdouble L2 = L*L, L3 = L*L*L;
      gdouble a1 = L2*gfs_cell_volume (p->scell, GFS_DOMAIN (sim))/sim->advection_params.dt;
      gdouble a2 = L2*gfs_cell_volume (p->ecell, GFS_DOMAIN (sim))/sim->advection_params.dt;

      /* secant-bisection root-finding: solves flow_rate(h, l, Q) - Q = 0 for the flow rate Q */
      p->Q = (*p->flow_rate) (z1, h1, z2, h2, l, g, p)/L3;
      gdouble Q1 = p->Q*2.;
      gdouble v1 = flow_rate_Q (z1, h1, z2, h2, l, g, p, a1, a2, Q1*L3)/L3 - Q1;
      gdouble Q2 = 0.;
      gdouble v2 = p->Q;
      if (fabs (v1) > DQ && fabs (v2) > DQ) {
	if (v1 > v2) {
	  gdouble v = v1;
	  v1 = v2; v2 = v;
	  v = Q1;
	  Q1 = Q2; Q2 = v;
	}
	if (v1*v2 >= 0.)
	  g_warning ("source_pipe_event: v1: %g v2: %g", v1*L3, v2*L3);
	else {
	  guint nitermax = 1000;
	  gdouble Qb;
	  p->Q = (v1*Q2 - v2*Q1)/(v1 - v2);
	  do {
	    Qb = p->Q;
	    gdouble v = flow_rate_Q (z1, h1, z2, h2, l, g, p, a1, a2, p->Q*L3)/L3 - p->Q;
	    if (v < 0.) {
	      v1 = v; Q1 = p->Q;
	    }
	    else {
	      v2 = v; Q2 = p->Q;
	    }
	    if (v2 > v1)
	      p->Q = (v1*Q2 - v2*Q1)/(v1 - v2);
	    nitermax--;
	  } while (fabs (p->Q - Qb) > DQ && nitermax);
	  if (nitermax == 0)
	    g_warning ("source_pipe_event: failed to converge! %g %g", 
		       p->Q*L3, fabs (p->Q - Qb)*L3);
	}
      }
    }
    return TRUE;
  }
  return FALSE;
}

static void source_pipe_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = source_pipe_read;
  GTS_OBJECT_CLASS (klass)->write = source_pipe_write;
  GFS_EVENT_CLASS (klass)->event = source_pipe_event;
}

static gdouble source_pipe_value (GfsSourceGeneric * s, 
				  FttCell * cell, 
				  GfsVariable * v)
{
  GfsSourcePipe * p = GFS_SOURCE_PIPE (s);
  if (cell == p->scell)
    return - p->Q/gfs_cell_volume (cell, v->domain);
  if (cell == p->ecell)
    return   p->Q/gfs_cell_volume (cell, v->domain);
  return 0.;
}

/* This is a simplistic flow rate model for a circular pipe. 
   The pipe is assumed to be always fully submerged. */
static double pipe_flow_rate (double z1, double h1, /* terrain elevation and flow depth at inlet */
			      double z2, double h2, /* terrain elevation and flow depth at outlet */
			      double l,             /* pipe length */
			      double g,             /* acceleration of gravity */
			      GfsSourcePipe * p)
{
  gdouble r = p->diameter/2.;
  gdouble A = M_PI*r*r; /* area */
  gdouble P = 2.*M_PI*r; /* perimeter */
  gdouble Rh = A/P; /* hydraulic radius */
  gdouble S = fabs (z1 + h1 - z2 - h2)/l; /* slope */
  gdouble n = 0.03; /* Gauckler-Manning coefficient */
  /* Gauckler-Manning-Strickler formula for the (signed) flow rate */
  return (z1 + h1 > z2 + h2 ? 1. : -1.)*A/n*pow (Rh, 2./3.)*sqrt (S);
}

static void source_pipe_init (GfsSourceGeneric * s)
{
  s->mac_value = s->centered_value = source_pipe_value;
  GFS_SOURCE_PIPE (s)->flow_rate = pipe_flow_rate;
}

GfsSourceGenericClass * gfs_source_pipe_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsSourcePipe",
      sizeof (GfsSourcePipe),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_pipe_class_init,
      (GtsObjectInitFunc) source_pipe_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsSourcePipe} */
