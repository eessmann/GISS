/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001 National Institute of Water and Atmospheric Research
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

#include <math.h>
#include <stdlib.h>

#include "timestep.h"
#include "source.h"
#include "solid.h"

/**
 * gfs_multilevel_params_write:
 * @par: the multilevel parameters.
 * @fp: a file pointer.
 *
 * Writes in @fp a text representation of the multilevel parameters
 * @par.  
 */
void gfs_multilevel_params_write (GfsMultilevelParams * par, FILE * fp)
{
  g_return_if_fail (par != NULL);
  g_return_if_fail (fp != NULL);

  fprintf (fp,
           "{\n"
	   "  tolerance = %g\n"
	   "  nrelax    = %u\n"
	   "  minlevel  = %u\n"
	   "  nitermax  = %u\n"
	   "}",
	   par->tolerance,
	   par->nrelax,
	   par->minlevel,
	   par->nitermax);
}

void gfs_multilevel_params_init (GfsMultilevelParams * par)
{
  g_return_if_fail (par != NULL);

  par->tolerance = 1e-3;
  par->nrelax    = 4;
  par->minlevel  = 0;
  par->nitermax  = 100;

  par->dimension = FTT_DIMENSION;
}

void gfs_multilevel_params_read (GfsMultilevelParams * par, GtsFile * fp)
{
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "tolerance", TRUE},
    {GTS_UINT,   "nrelax",    TRUE},
    {GTS_UINT,   "minlevel",  TRUE},
    {GTS_UINT,   "nitermax",  TRUE},
    {GTS_NONE}
  };

  g_return_if_fail (par != NULL);
  g_return_if_fail (fp != NULL);

  var[0].data = &par->tolerance;
  var[1].data = &par->nrelax;
  var[2].data = &par->minlevel;
  var[3].data = &par->nitermax;

  gfs_multilevel_params_init (par);
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (par->tolerance <= 0.) {
    gts_file_variable_error (fp, var, "tolerance",
			     "tolerance `%g' must be strictly positive",
			     par->tolerance);
    return;
  }
  if (par->nrelax == 0)
    gts_file_variable_error (fp, var, "nrelax", "nrelax must be non zero");
}

static void reset_gradients (FttCell * cell)
{
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++)
    GFS_STATE (cell)->g[c] = 0.;
}

static void correct_normal_velocity (FttCellFace * face,
				     gpointer * data)
{
  GfsGradient g;
  gdouble dp;
  FttFaceType type;
  GfsStateVector * s;
  GfsVariable * p = data[0];
  gdouble * dt = data[1];
  FttComponent c;

  if (GFS_FACE_FRACTION (face) == 0.)
    return;

  s = GFS_STATE (face->cell);
  type = ftt_face_type (face);
  c = face->d/2;

  //  gfs_face_gradient_flux_centered (face, &g, GFS_P, -1);
  gfs_face_weighted_gradient (face, &g, p->i, -1);
  dp = (g.b - g.a*s->p)/ftt_cell_size (face->cell);
  if (!FTT_FACE_DIRECT (face))
    dp = - dp;

  if (s->solid && s->solid->s[face->d] > 0.)
    dp /= s->solid->s[face->d];

  GFS_FACE_NORMAL_VELOCITY_LEFT (face) -= dp*(*dt);
  s->g[c] += dp;

  switch (type) {
  case FTT_FINE_FINE:
    GFS_FACE_NORMAL_VELOCITY_RIGHT (face) -= dp*(*dt);
    GFS_STATE (face->neighbor)->g[c] += dp;
    break;
  case FTT_FINE_COARSE: {
    /* fixme: does this work (FTT_CELLS/2?) for 2D3? */
    dp *= GFS_FACE_FRACTION_LEFT (face)/(GFS_FACE_FRACTION_RIGHT (face)*FTT_CELLS/2);
    GFS_FACE_NORMAL_VELOCITY_RIGHT (face) -= dp*(*dt);
    GFS_STATE (face->neighbor)->g[c] += dp;
    break;
  }
  default:
    g_assert_not_reached ();
  }
}

static void scale_gradients (FttCell * cell)
{
  FttComponent c;
  FttCellNeighbors n;

  ftt_cell_neighbors (cell, &n);
  for (c = 0; c < FTT_DIMENSION; c++) {
    FttCell * c1 = n.c[2*c], * c2 = n.c[2*c + 1];

    if (c1 && c2 && !GFS_CELL_IS_BOUNDARY (c1) && !GFS_CELL_IS_BOUNDARY (c2))
      GFS_STATE (cell)->g[c] /= 2.;
  }
}

/**
 * gfs_correct_normal_velocities:
 * @domain: a #GfsDomain.
 * @dimension: the number of dimensions (2 or 3).
 * @p: the pressure field.
 * @dt: the timestep.
 *
 * Corrects the normal velocity field of @domain using @p and and @dt.
 *
 * Also fills the g[] field with the centered gradient of @p.
 */
void gfs_correct_normal_velocities (GfsDomain * domain,
				    guint dimension,
				    GfsVariable * p,
				    gdouble dt)
{
  gpointer data[2];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (p != NULL);
  
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) reset_gradients, NULL);
  data[0] = p;
  data[1] = &dt;
  gfs_domain_face_traverse (domain, dimension == 2 ? FTT_XY : FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) correct_normal_velocity, data);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) scale_gradients, NULL);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, gfs_gx);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, gfs_gy);
#if (!FTT_2D)
  if (dimension == 3)
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, gfs_gz);
#endif /* 2D3 or 3D */
}

static void scale_divergence (FttCell * cell, gdouble * a)
{
  GFS_STATE (cell)->div /= *a;
  GFS_STATE (cell)->g[0] = 0.;
}

/**
 * gfs_mac_projection:
 * @domain: a #GfsDomain.
 * @par: the projection control parameters.
 * @apar: the advection parameters.
 *
 * Corrects the face-centered velocity field (MAC field) on the leaf
 * level of @domain using an exact (MAC) projection. The resulting
 * face-centered velocity field is (almost) exactly divergence
 * free. The (potential) pressure field is also obtained as a
 * by-product as well as its gradient at the center of the leaf cells
 * of the domain (the gradient is stored in the %GFS_G variables and is
 * obtained by simple averaging from the face values to the center).
 *
 * The @residual field of the @par projection parameters is set to the
 * norm of the residual after the projection. The @niter field of the
 * @par projection parameters is set to the number of iterations
 * performed to solve the Poisson equation. The other projection
 * parameters are not modified.
 */
void gfs_mac_projection (GfsDomain * domain,
			 GfsMultilevelParams * par,
			 GfsAdvectionParams * apar)
{
  guint minlevel, maxlevel;
  gdouble dt;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (apar != NULL);

  gfs_domain_timer_start (domain, "mac_projection");
  
  apar->v = gfs_variable_from_name (domain->variables, "U");
  dt = apar->dt;
  apar->dt /= 2.;

  /* Initialize face coefficients */
  gfs_poisson_coefficients (domain, apar->c, apar->rho);

  /* compute MAC divergence */
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_normal_divergence, 
			    NULL);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
  			    (FttCellTraverseFunc) scale_divergence, &apar->dt);

#if 0
  {
    FILE * fp = fopen ("/tmp/mac", "wt");
    GfsNorm norm;

    gfs_write_mac_velocity (domain, 0.9, FTT_TRAVERSE_LEAFS, -1, NULL, fp);
    fclose (fp);
    norm = gfs_domain_norm_variable (domain, gfs_div, FTT_TRAVERSE_LEAFS, -1);
    fprintf (stderr, "mac div before: %g %g %g\n",
	     norm.first, norm.second, norm.infty);
  }
#endif

  /* solve for pressure */
  minlevel = domain->rootlevel;
  if (par->minlevel > minlevel)
    minlevel = par->minlevel;
  maxlevel = gfs_domain_depth (domain);
  gfs_residual (domain, par->dimension, FTT_TRAVERSE_LEAFS, -1, gfs_p, gfs_div, gfs_res);
  par->residual_before = par->residual = 
    gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, apar->dt);
  par->niter = 0;
  while (par->residual.infty > par->tolerance && 
	 par->niter < par->nitermax) {
    gfs_poisson_cycle (domain, par->dimension, minlevel, maxlevel, par->nrelax, gfs_p, gfs_div);
    par->residual = gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1,
					      apar->dt);
    par->niter++;
  }
  
  gfs_correct_normal_velocities (domain, FTT_DIMENSION, gfs_p, apar->dt);

#if 0
  {
    FILE * fp = fopen ("/tmp/macafter", "wt");
    GfsNorm norm;

    gfs_write_mac_velocity (domain, 0.9, FTT_TRAVERSE_LEAFS, -1, NULL, fp);
    fclose (fp);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_normal_divergence, 
			      NULL);
    norm = gfs_domain_norm_variable (domain, gfs_div, FTT_TRAVERSE_LEAFS, -1);
    fprintf (stderr, "mac div after: %g %g %g\n",
	     norm.first, norm.second, norm.infty);
  }
#endif

  apar->dt = dt;

  gfs_domain_timer_stop (domain, "mac_projection");
}

static void correct_2D (FttCell * cell, gdouble * dt)
{
  GfsStateVector * s = GFS_STATE (cell);

  s->u -= s->g[0]*(*dt);
  s->v -= s->g[1]*(*dt);
}

#if (!FTT_2D)
static void correct_3D (FttCell * cell, gdouble * dt)
{
  GfsStateVector * s = GFS_STATE (cell);

  s->u -= s->g[0]*(*dt);
  s->v -= s->g[1]*(*dt);
  s->w -= s->g[2]*(*dt);
}
#endif /* 2D3 or 3D */

/**
 * gfs_correct_centered_velocities:
 * @domain: a #GfsDomain.
 * @dimension: the number of dimensions (2 or 3).
 * @dt: the timestep.
 *
 * Corrects the velocity field of @domain using the pressure gradient
 * stored in g[].
 */
void gfs_correct_centered_velocities (GfsDomain * domain,
				      guint dimension,
				      gdouble dt)
{
  GfsVariable * v;
  FttComponent c;

  g_return_if_fail (domain != NULL);

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc)
#if FTT_2D
			    correct_2D,
#else  /* 2D3 or 3D */
			    (dimension == 2 ? correct_2D : correct_3D),
#endif /* 2D3 or 3D */
			    &dt);
  v = gfs_variable_from_name (domain->variables, "U");
  for (c = 0; c < dimension; c++, v = v->next)
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, v);
}

/**
 * gfs_approximate_projection:
 * @domain: a #GfsDomain.
 * @par: the projection control parameters.
 * @apar: the advection parameters.
 *
 * Corrects the centered velocity field on the leaf level of @domain
 * using an approximate projection. The resulting centered velocity
 * field is approximately divergence free. The (potential) pressure
 * field is also obtained as a by-product.
 *
 * The @residual field of the @par projection parameters is set to the
 * norm of the residual (on the MAC grid) after the projection. The
 * @niter field of the @par projection parameters is set to the number
 * of iterations performed to solve the Poisson equation. The other
 * projection parameters are not modified.
 *
 * The Poisson equation for the pressure is first solved on a MAC grid
 * where the MAC velocities are obtained from the centered velocities
 * by simple averaging. The resulting pressure gradients (defined on
 * the faces) are then averaged down on the center of the cells to
 * correct the centered velocity.  
 */
void gfs_approximate_projection (GfsDomain * domain,
				 GfsMultilevelParams * par,
				 GfsAdvectionParams * apar)
{
  guint minlevel, maxlevel;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (apar != NULL);

  gfs_domain_timer_start (domain, "approximate_projection");

  /* Initialize face coefficients */
  gfs_poisson_coefficients (domain, apar->c, apar->rho);

  /* compute MAC velocities from centered velocities */
  gfs_domain_face_traverse (domain, FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
      (FttFaceTraverseFunc) gfs_face_reset_normal_velocity,
			    NULL);
  gfs_domain_face_traverse (domain, FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
      (FttFaceTraverseFunc) gfs_face_interpolated_normal_velocity, NULL);

#if 0
  {
    FILE * fp = fopen ("/tmp/macapprox", "wt");
    GfsNorm norm;

    gfs_write_mac_velocity (domain, 0.9, FTT_TRAVERSE_LEAFS, -1, NULL, fp);
    fclose (fp);
    norm = gfs_domain_norm_variable (domain, gfs_div, FTT_TRAVERSE_LEAFS, -1);
    fprintf (stderr, "mac div before: %g %g %g\n",
	     norm.first, norm.second, norm.infty);
  }
#endif

  /* compute MAC divergence */
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_normal_divergence, 
			    NULL);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
  			    (FttCellTraverseFunc) scale_divergence, &apar->dt);
  
  /* solve for pressure */
  minlevel = domain->rootlevel;
  if (par->minlevel > minlevel)
    minlevel = par->minlevel;
  maxlevel = gfs_domain_depth (domain);
  gfs_residual (domain, par->dimension, FTT_TRAVERSE_LEAFS, -1, gfs_p, gfs_div, gfs_res);
  par->residual_before = par->residual = 
    gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, apar->dt);
  par->niter = 0;
  while (par->residual.infty > par->tolerance && 
	 par->niter < par->nitermax) {
#if 0
    fprintf (stderr, "%d bias: %g first: %g second: %g infty: %g\n",
	     par->niter, 
	     par->residual.bias, 
	     par->residual.first, 
	     par->residual.second, 
	     par->residual.infty);
#endif
    gfs_poisson_cycle (domain, par->dimension, minlevel, maxlevel, par->nrelax, gfs_p, gfs_div);
    par->residual = gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, apar->dt);
    par->niter++;
  }

  gfs_correct_normal_velocities (domain, FTT_DIMENSION, gfs_p, apar->dt);
  gfs_correct_centered_velocities (domain, FTT_DIMENSION, apar->dt);

  gfs_domain_timer_stop (domain, "approximate_projection");
}

/**
 * gfs_predicted_face_velocities:
 * @domain: a #GfsDomain.
 * @d: the number of dimensions (2 or 3).
 * @par: the advection parameters.
 *
 * Fills the face (MAC) normal velocities of each leaf cell of @domain
 * with the predicted values at time t + dt/2 using a godunov type
 * advection scheme.  
 */
void gfs_predicted_face_velocities (GfsDomain * domain,
				    guint d,
				    GfsAdvectionParams * par)
{
  FttComponent c;
  FttCellTraverseFunc face_values;
  GfsUpwinding upwinding;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);

  gfs_domain_timer_start (domain, "predicted_face_velocities");

  gfs_domain_face_traverse (domain, d == 2 ? FTT_XY : FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
      (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, 
			    NULL);
  par->use_centered_velocity = TRUE;
  if (par->scheme == GFS_NONE) {
    face_values = (FttCellTraverseFunc) gfs_cell_non_advected_face_values;
    upwinding = GFS_NO_UPWINDING;
  }
  else {
    face_values = (FttCellTraverseFunc) gfs_cell_advected_face_values;
    upwinding = GFS_CENTERED_UPWINDING;
  }
  par->v = gfs_variable_from_name (domain->variables, "U");
  for (c = 0; c < d; c++, par->v = par->v->next) {
    gfs_domain_cell_traverse (domain, 
    			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
    			      face_values, par);
    gfs_domain_face_bc (domain, c, par->v);
    gfs_domain_face_traverse (domain, c,
    			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) gfs_face_advected_normal_velocity, 
			      &upwinding);
  }

  gfs_domain_timer_stop (domain, "predicted_face_velocities");
}

/**
 * gfs_diffusion:
 * @domain: a #GfsDomain.
 * @par: the multilevel parameters.
 * @v: a #GfsVariable.
 *
 * Solves a diffusion equation for variable @v using a Crank-Nicholson
 * scheme with multilevel relaxations.
 *
 * Diffusion coefficients must have been set using
 * gfs_diffusion_coefficients() and a right-hand side defined using
 * calls to gfs_diffusion_rhs().
 */
void gfs_diffusion (GfsDomain * domain,
		    GfsMultilevelParams * par,
		    GfsVariable * v)
{
  guint minlevel, maxlevel;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (v != NULL);

  minlevel = domain->rootlevel;
  if (par->minlevel > minlevel)
    minlevel = par->minlevel;
  maxlevel = gfs_domain_depth (domain);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_diffusion_residual, v);
  par->residual_before = par->residual = 
    gfs_domain_norm_variable (domain, gfs_res, FTT_TRAVERSE_LEAFS, -1);
  par->niter = 0;
  while (par->residual.infty > par->tolerance && 
	 par->niter < par->nitermax) {
    gfs_diffusion_cycle (domain, minlevel, maxlevel, par->nrelax, v);
    par->residual = gfs_domain_norm_variable (domain, gfs_res, 
					      FTT_TRAVERSE_LEAFS, -1);
    par->niter++;
  }
}

static GfsSourceDiffusion * source_diffusion (GfsVariable * v)
{
  if (v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
    
    while (i) {
      GtsObject * o = i->data;
      
      if (GFS_IS_SOURCE_DIFFUSION (o) && !GFS_IS_SOURCE_DIFFUSION_EXPLICIT (o))
        return GFS_SOURCE_DIFFUSION (o);
      i = i->next;
    }
  }
  return NULL;
}

static void variable_sources (GfsDomain * domain,
			      GfsAdvectionParams * par,
			      GfsVariable * sv)
{
  if (par->scheme == GFS_GODUNOV) {
    GfsVariable * v = par->v;

    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, par->fv);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_advected_face_values, par);
    gfs_domain_face_bc (domain, FTT_XYZ, par->v);
    gfs_domain_face_traverse (domain, FTT_XYZ,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) par->flux, par);
    par->v = sv;
    gfs_domain_traverse_merged (domain, (GfsMergedTraverseFunc) gfs_advection_update, par);
    par->v = v;
  }
  /* fixme: time should be set to t + dt/2 here for evaluation of
     source terms in the call below */
  gfs_domain_variable_centered_sources (domain, par->v, sv, par->dt);
}

static void variable_diffusion (GfsDomain * domain,
				GfsSourceDiffusion * d,
				GfsAdvectionParams * par,
				GfsMultilevelParams * dpar,
				GfsVariable * c,
				gdouble rho)
{
  if (c != NULL)
    gfs_viscosity_coefficients (domain, d, par->dt, c, rho);
  else
    gfs_diffusion_coefficients (domain, d, par->dt);
  gfs_domain_surface_bc (domain, par->v);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_diffusion_rhs, par->v);
  /* fixme: time shoud be set to t + dt here in case boundary values are
     time-dependent in the call below */
  gfs_domain_surface_bc (domain, par->v);
  gfs_diffusion (domain, dpar, par->v);
}

/**
 * gfs_centered_velocity_advection_diffusion:
 * @domain: a #GfsDomain.
 * @dimension: the number of dimensions (2 or 3).
 * @apar: the advection parameters.
 * @dpar: the multilevel solver parameters for the diffusion equation.
 *
 * Advects the (centered) velocity field using the current
 * face-centered (MAC) velocity field and @par->flux to compute the
 * velocity flux through the faces of each cell.
 *
 * For each component of the velocity, before calling the @par->flux
 * function the face values are first defined (at time t + dt/2) and
 * can then be used within the @par->flux function.
 *
 * "Small" cut cells are treated using a cell-merging approach to
 * avoid any restrictive CFL stability condition.  
 */
void gfs_centered_velocity_advection_diffusion (GfsDomain * domain,
						guint dimension,
						GfsAdvectionParams * apar,
						GfsMultilevelParams * dpar)
{
  FttComponent c;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (apar != NULL);
  g_return_if_fail (dpar != NULL);

  gfs_domain_timer_start (domain, "centered_velocity_advection_diffusion");

  apar->use_centered_velocity = FALSE;
  apar->v = gfs_variable_from_name (domain->variables, "U");
  for (c = 0; c < dimension; c++, apar->v = apar->v->next) {
    GfsSourceDiffusion * d = source_diffusion (apar->v);

    if (d) {
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) gfs_cell_reset, gfs_div);
      variable_sources (domain, apar, gfs_div);
      variable_diffusion (domain, d, apar, dpar, apar->c, apar->rho);
    }
    else {
      variable_sources (domain, apar, apar->v);
      gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, apar->v);
    }
  }

  gfs_domain_timer_stop (domain, "centered_velocity_advection_diffusion");
}

static void save_previous (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsVariable * vh = data[1];

  GFS_VARIABLE (cell, vh->i) = GFS_VARIABLE (cell, v->i);
}

static void average_previous (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsVariable * vh = data[1];

  GFS_VARIABLE (cell, vh->i) = (GFS_VARIABLE (cell, vh->i) +
				GFS_VARIABLE (cell, v->i))/2.;
}

/**
 * gfs_tracer_advection_diffusion:
 * @domain: a #GfsDomain.
 * @par: the advection parameters.
 * @dpar: the multilevel solver parameters for the diffusion equation.
 * @half: a #GfsVariable or %NULL.
 *
 * Advects the @v field of @par using the current face-centered (MAC)
 * velocity field.
 *
 * If @half is not %NULL, the half-timestep value of @par->v is
 * stored in the corresponding variable.  
 */
void gfs_tracer_advection_diffusion (GfsDomain * domain,
				     GfsAdvectionParams * par,
				     GfsMultilevelParams * dpar,
				     GfsVariable * half)
{
  gpointer data[2];
  GfsSourceDiffusion * d;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (dpar != NULL);

  gfs_domain_timer_start (domain, "tracer_advection_diffusion");

  if (half) {
    data[0] = par->v;
    data[1] = half;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) save_previous, data);
  }

  if ((d = source_diffusion (par->v))) {
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, gfs_div);
    variable_sources (domain, par, gfs_div);
    variable_diffusion (domain, d, par, dpar, NULL, 0.);
  }
  else {
    variable_sources (domain, par, par->v);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, par->v);
  }

  if (half) {
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) average_previous, data);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, half);
  }

  gfs_domain_timer_stop (domain, "tracer_advection_diffusion");
}

/* GfsSurfaceGenericBc: Object */

static void gfs_surface_generic_bc_destroy (GtsObject * o)
{
  if (GFS_SURFACE_GENERIC_BC (o)->v)
    GFS_SURFACE_GENERIC_BC (o)->v->surface_bc = NULL;

  (* GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ())->parent_class->destroy) (o);
}

static void gfs_surface_generic_bc_read (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  GfsSurfaceGenericBc * bc = GFS_SURFACE_GENERIC_BC (*o);
  GtsObjectClass * klass;

  if (GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a class name");
    return;
  }
  if (!(klass = gfs_object_class_from_name (fp->token->str))) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_surface_generic_bc_class ())) {
    gts_file_error (fp, "class `%s' is not a GfsSurfaceGenericClass", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable name");
    return;
  }
  bc->v = gfs_variable_from_name (domain->variables, fp->token->str);
  if (!bc->v) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  if (bc->v->surface_bc) {
    gts_file_error (fp, "variable `%s' already has a surface boundary condition", 
		    fp->token->str);
    return;
  }
  bc->v->surface_bc = bc;
  gts_file_next_token (fp);
}

static void gfs_surface_generic_bc_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ())->parent_class->write) (o, fp);
  fprintf (fp, "%s %s", o->klass->info.name, GFS_SURFACE_GENERIC_BC (o)->v->name);
}

static void gfs_surface_generic_bc_class_init (GfsSurfaceGenericBcClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gfs_surface_generic_bc_destroy;
  GTS_OBJECT_CLASS (klass)->read = gfs_surface_generic_bc_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_surface_generic_bc_write;
}

GfsSurfaceGenericBcClass * gfs_surface_generic_bc_class (void)
{
  static GfsSurfaceGenericBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_surface_generic_bc_info = {
      "GfsSurfaceGenericBc",
      sizeof (GfsSurfaceGenericBc),
      sizeof (GfsSurfaceGenericBcClass),
      (GtsObjectClassInitFunc) gfs_surface_generic_bc_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &gfs_surface_generic_bc_info);
  }

  return klass;
}

/* GfsSurfaceBc: Object */

static void gfs_surface_bc_destroy (GtsObject * object)
{
  gts_object_destroy (GTS_OBJECT (GFS_SURFACE_BC (object)->type));
  gts_object_destroy (GTS_OBJECT (GFS_SURFACE_BC (object)->val));

  (* GTS_OBJECT_CLASS (gfs_surface_bc_class ())->parent_class->destroy) (object);
}

static void gfs_surface_bc_read (GtsObject ** o, GtsFile * fp)
{
  GfsSurfaceBc * bc = GFS_SURFACE_BC (*o);

  if (GTS_OBJECT_CLASS (gfs_surface_bc_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_surface_bc_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  if (!strcmp (fp->token->str, "Neumann")) {
    bc->type->val = 0.;
    gts_file_next_token (fp);
  }
  else if (!strcmp (fp->token->str, "Dirichlet")) {
    bc->type->val = 1.;
    gts_file_next_token (fp);
  }
  else {
    gfs_object_simulation (bc->type) = gfs_object_simulation (bc);
    gfs_function_read (bc->type, fp);
    if (fp->type == GTS_ERROR)
      return;
  }
  gfs_object_simulation (bc->val) = gfs_object_simulation (bc);
  gfs_function_read (bc->val, fp);
}

static void gfs_surface_bc_write (GtsObject * o, FILE * fp)
{
  GfsSurfaceBc * bc = GFS_SURFACE_BC (o);

  if (GTS_OBJECT_CLASS (gfs_surface_bc_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_surface_bc_class ())->parent_class->write) (o, fp);
  if (!bc->type->f)
    fprintf (fp, " %s", bc->type->val ? "Dirichlet" : "Neumann");
  else
    gfs_function_write (bc->type, fp);
  gfs_function_write (bc->val, fp);
}

static void gfs_surface_bc_bc (FttCell * cell, GfsSurfaceGenericBc * b)
{
  GfsSurfaceBc * bc = GFS_SURFACE_BC (b);
  gdouble t = gfs_object_simulation (bc)->time.t;
  gdouble val = gfs_function_value (bc->val, &GFS_STATE (cell)->solid->ca, t);

  if (gfs_function_value (bc->type, &GFS_STATE (cell)->solid->ca, t) > 0.) {
    cell->flags |= GFS_FLAG_DIRICHLET;
    GFS_STATE (cell)->solid->fv = val;
  }
  else {
    cell->flags &= ~GFS_FLAG_DIRICHLET;
    GFS_STATE (cell)->solid->fv = val; /* fixme: scaling is probably wrong */
  }
}

static void gfs_surface_bc_class_init (GfsSurfaceGenericBcClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_surface_bc_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_surface_bc_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_surface_bc_destroy;
  klass->bc = gfs_surface_bc_bc;
}

static void gfs_surface_bc_init (GfsSurfaceBc * object)
{
  object->type = gfs_function_new (gfs_function_class (), 0.);
  object->val  = gfs_function_new (gfs_function_class (), 0.);
}

GfsSurfaceGenericBcClass * gfs_surface_bc_class (void)
{
  static GfsSurfaceGenericBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_surface_bc_info = {
      "GfsSurfaceBc",
      sizeof (GfsSurfaceBc),
      sizeof (GfsSurfaceGenericBcClass),
      (GtsObjectClassInitFunc) gfs_surface_bc_class_init,
      (GtsObjectInitFunc) gfs_surface_bc_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ()),
				  &gfs_surface_bc_info);
  }

  return klass;
}
