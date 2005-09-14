/* Gerris - The GNU Flow Solver
 * Copyright (C) 2004 Stéphane Popinet
 * National Institute of Water and Atmospheric Research
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

#include "ocean.h"
#include "timestep.h"
#include "adaptive.h"
#include "source.h"
#include "vof.h"
#include "graphic.h"

static void reset_gradients (FttCell * cell, gpointer * data)
{
  GfsVariable ** g = data[0];
  guint * dimension = data[1];    
  FttComponent c;

  for (c = 0; c < *dimension; c++)
    GFS_VARIABLE (cell, g[c]->i) = 0.;
}

static void correct_normal_velocity_weighted (FttCellFace * face,
					      gpointer * data)
{
  GfsGradient g;
  gdouble dp;
  FttFaceType type;
  GfsStateVector * s;
  GfsVariable * p = data[0];
  GfsVariable ** gv = data[1];
  gdouble * dt = data[2];
  FttComponent c;

  if (GFS_FACE_FRACTION (face) == 0.)
    return;

  s = GFS_STATE (face->cell);
  type = ftt_face_type (face);
  c = face->d/2;

  gfs_face_weighted_gradient (face, &g, p->i, -1);
  dp = (g.b - g.a*GFS_VARIABLE (face->cell, p->i))/ftt_cell_size (face->cell);
  if (!FTT_FACE_DIRECT (face))
    dp = - dp;

  if (s->solid && s->solid->s[face->d] > 0.)
    dp /= s->solid->s[face->d];

  GFS_FACE_NORMAL_VELOCITY_LEFT (face) -= dp*(*dt);
  GFS_VARIABLE (face->cell, gv[c]->i) += dp*GFS_FACE_FRACTION_LEFT (face);

  switch (type) {
  case FTT_FINE_FINE:
    GFS_FACE_NORMAL_VELOCITY_RIGHT (face) -= dp*(*dt);
    GFS_VARIABLE (face->neighbor, gv[c]->i) += dp*GFS_FACE_FRACTION_RIGHT (face);
    break;
  case FTT_FINE_COARSE: {
    dp *= GFS_FACE_FRACTION_LEFT (face)/(FTT_CELLS/2);
    GFS_VARIABLE (face->neighbor, gv[c]->i) += dp;
    g_assert (GFS_FACE_FRACTION_RIGHT (face) > 0.);
    GFS_FACE_NORMAL_VELOCITY_RIGHT (face) -= dp/GFS_FACE_FRACTION_RIGHT (face)*(*dt);
    break;
  }
  default:
    g_assert_not_reached ();
  }
}

static void scale_gradients_weighted (FttCell * cell, gpointer * data)
{
  GfsVariable ** g = data[0];
  guint * dimension = data[1];
  FttComponent c;

  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * s = GFS_STATE (cell)->solid;

    for (c = 0; c < *dimension; c++) {
      g_assert (s->s[2*c] + s->s[2*c + 1] > 0.);
      GFS_VARIABLE (cell, g[c]->i) /= s->s[2*c] + s->s[2*c + 1];
    }
  }
  else {
    FttCellNeighbors n;

    ftt_cell_neighbors (cell, &n);
    for (c = 0; c < *dimension; c++) {
      FttCell * c1 = n.c[2*c], * c2 = n.c[2*c + 1];

      if (c1 && c2 && !GFS_CELL_IS_GRADIENT_BOUNDARY (c1) && !GFS_CELL_IS_GRADIENT_BOUNDARY (c2))
	GFS_VARIABLE (cell, g[c]->i) /= 2.;
    }
  }
}

/**
 * gfs_correct_normal_velocities_weighted:
 * @domain: a #GfsDomain.
 * @dimension: the number of dimensions (2 or 3).
 * @p: the pressure field.
 * @g: where to store the pressure gradient.
 * @dt: the timestep.
 * @weighted: whether to use fraction-weighting or not.
 *
 * Corrects the normal velocity field of @domain using @p and and @dt.
 *
 * Also allocates the @g variables and fills them with the centered gradient of @p.
 */
static void gfs_correct_normal_velocities_weighted (GfsDomain * domain,
						    guint dimension,
						    GfsVariable * p,
						    GfsVariable ** g,
						    gdouble dt,
						    gboolean weighted)
{
  if (!weighted)
    gfs_correct_normal_velocities (domain, dimension, p, g, dt);
  else {
    gpointer data[3];
    FttComponent c;
    
    g_return_if_fail (domain != NULL);
    g_return_if_fail (p != NULL);
    g_return_if_fail (g != NULL);
    
    for (c = 0; c < dimension; c++) {
      g[c] = gfs_temporary_variable (domain);
      gfs_variable_set_vector (g[c], c);
    }
    data[0] = g;
    data[1] = &dimension;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) reset_gradients, data);
    data[0] = p;
    data[1] = g;
    data[2] = &dt;
    gfs_domain_face_traverse (domain, dimension == 2 ? FTT_XY : FTT_XYZ,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) correct_normal_velocity_weighted, data);
    data[0] = g;
    data[1] = &dimension;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) scale_gradients_weighted, data);
    for (c = 0; c < dimension; c++)
      gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, g[c]);
  }
}

/* GfsOcean: Object */

static void ocean_destroy (GtsObject * object)
{
  guint i;
  GPtrArray * layer = GFS_OCEAN (object)->layer;

  for (i = 0; i < layer->len; i++) {
    GfsDomain * d = g_ptr_array_index (layer, i);
    d->allocated = g_array_new (FALSE, TRUE, sizeof (gboolean));
    gts_object_destroy (GTS_OBJECT (d));
  }
  g_ptr_array_free (layer, TRUE);

  (* GTS_OBJECT_CLASS (gfs_ocean_class ())->parent_class->destroy) (object);  
}

static void new_layer (GfsOcean * ocean)
{
  GfsDomain * domain = GFS_DOMAIN (ocean);
  GfsDomain * d = GFS_DOMAIN (gts_object_new (GTS_OBJECT_CLASS (gfs_domain_class ())));
  
  d->rootlevel = domain->rootlevel;
  d->refpos = domain->refpos;
  d->lambda = domain->lambda;
  g_array_free (d->allocated, TRUE);
  d->allocated = domain->allocated;
  g_ptr_array_add (ocean->layer, d);
}

static void add_layer (GfsBox * box, GfsDomain * domain)
{
#if FTT_2D
  gts_container_add (GTS_CONTAINER (g_ptr_array_index (GFS_OCEAN (domain)->layer, 0)),
		     GTS_CONTAINEE (box));
#else    /* 2D3 or 3D */
  if (box->neighbor[FTT_FRONT] == NULL || GFS_IS_BOUNDARY (box->neighbor[FTT_FRONT])) {
    GPtrArray * layer = GFS_OCEAN (domain)->layer;
    GtsObject * n;
    guint l = 0;

    gts_container_add (GTS_CONTAINER (g_ptr_array_index (layer, l++)), GTS_CONTAINEE (box));
    n = box->neighbor[FTT_BACK];
    while (GFS_IS_BOX (n)) {
      if (l == layer->len)
	new_layer (GFS_OCEAN (domain));
      gts_container_add (GTS_CONTAINER (g_ptr_array_index (layer, l++)), GTS_CONTAINEE (n));
      n = GFS_BOX (n)->neighbor[FTT_BACK];
    }
  }
#endif /* 2D3 or 3D */
}

static void ocean_post_read (GfsDomain * domain, GtsFile * fp)
{
  (* GFS_DOMAIN_CLASS (GTS_OBJECT_CLASS (gfs_ocean_class ())->parent_class)->post_read) 
    (domain, fp);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) add_layer, domain);
  g_assert (GFS_OCEAN (domain)->layer->len > 0);
  GFS_OCEAN (domain)->toplayer = g_ptr_array_index (GFS_OCEAN (domain)->layer, 0);
}

#if (!FTT_2D) /* 2D3 or 3D */
static void sum_coeff (FttCell * cell)
{
  FttCell * c = ftt_cell_neighbor (cell, FTT_BACK);
  guint level = ftt_cell_level (cell);
  FttDirection d;

  while (c) {
    g_assert (FTT_CELL_IS_LEAF (c) && ftt_cell_level (c) == level);
    for (d = 0; d < FTT_NEIGHBORS_2D; d++) {
      g_assert (ftt_cell_neighbor (c, d) || GFS_STATE (c)->f[d].v == 0.);
      g_assert (!GFS_IS_MIXED (c) || GFS_STATE (c)->solid->s[d] != 0. ||
		GFS_STATE (c)->f[d].v == 0.);
      GFS_STATE (cell)->f[d].v += GFS_STATE (c)->f[d].v;
    }
    c = ftt_cell_neighbor (c, FTT_BACK);
  }
}

static void face_coeff_from_below (FttCell * cell)
{
  FttDirection d;
  GfsFaceStateVector * f = GFS_STATE (cell)->f;

  for (d = 0; d < FTT_NEIGHBORS_2D; d++) {
    FttCellChildren child;
    guint i, n;

    f[d].v = 0.;
    n = ftt_cell_children_direction (cell, d, &child);
    for (i = 0; i < n; i++)
      if (child.c[i])
	f[d].v += GFS_STATE (child.c[i])->f[d].v;
    f[d].v /= n;
  }
}

static void sum_divergence (FttCell * cell, GfsVariable * div)
{
  FttCell * c = ftt_cell_neighbor (cell, FTT_BACK);
  guint level = ftt_cell_level (cell);

  while (c) {
    g_assert (FTT_CELL_IS_LEAF (c) && ftt_cell_level (c) == level);
    GFS_VARIABLE (cell, div->i) += GFS_VARIABLE (c, div->i);
    c = ftt_cell_neighbor (c, FTT_BACK);
  }
}

static void column_pressure (FttCell * cell, GfsVariable * p)
{
  FttCell * c = ftt_cell_neighbor (cell, FTT_BACK);
  guint level = ftt_cell_level (cell);

  while (c) {
    g_assert (FTT_CELL_IS_LEAF (c) && ftt_cell_level (c) == level);
    GFS_VARIABLE (c, p->i) = GFS_VARIABLE (cell, p->i);
    c = ftt_cell_neighbor (c, FTT_BACK);
  }
}

static void compute_w (FttCell * c, GfsVariable * W)
{
  guint level = ftt_cell_level (c);
  gdouble wf = 0., w = 0.;

  while (c) {
    GfsStateVector * s = GFS_STATE (c);

    g_assert (FTT_CELL_IS_LEAF (c) && ftt_cell_level (c) == level);
    s->f[FTT_BACK].un = w;
    wf += (s->f[FTT_LEFT].v*s->f[FTT_LEFT].un - s->f[FTT_RIGHT].v*s->f[FTT_RIGHT].un +
    	   s->f[FTT_BOTTOM].v*s->f[FTT_BOTTOM].un - s->f[FTT_TOP].v*s->f[FTT_TOP].un);
    if (GFS_IS_MIXED (c))
      s->f[FTT_FRONT].un = w = GFS_STATE (c)->solid->s[FTT_FRONT] > 0. ? 
	wf/GFS_STATE (c)->solid->s[FTT_FRONT] : 0.;
    else
      s->f[FTT_FRONT].un = w = wf;
    GFS_VARIABLE (c, W->i) = (s->f[FTT_BACK].un + s->f[FTT_FRONT].un)/2.;
    c = ftt_cell_neighbor (c, FTT_FRONT);
  }
}
#endif /* 2D3 or 3D */

#define THETA 0.5

typedef struct {
  GfsVariable * pn, * div, * divn, * dia;
  gdouble dt, G;
} FreeSurfaceParams;

static void normal_divergence (FttCell * cell, FreeSurfaceParams * p)
{
  gfs_normal_divergence_2D (cell, p->div);
  GFS_VARIABLE (cell, p->div->i) += (1. - THETA)*GFS_VARIABLE (cell, p->divn->i)/THETA;
}

static void scale_divergence_helmoltz (FttCell * cell, FreeSurfaceParams * p)
{
  gdouble h = ftt_cell_size (cell);
  gdouble c = 2.*h*h/(THETA*p->G*p->dt*p->dt);

  if (GFS_IS_MIXED (cell))
#if FTT_2D
    c *= GFS_STATE (cell)->solid->a;
#else  /* 2D3 or 3D */
    c *= GFS_STATE (cell)->solid->s[FTT_FRONT];
#endif /* 2D3 or 3D */

  GFS_VARIABLE (cell, p->dia->i) = c;
  GFS_VARIABLE (cell, p->div->i) = 2.*GFS_VARIABLE (cell, p->div->i)/p->dt -
    c*GFS_VARIABLE (cell, p->pn->i);
}

/**
 * gfs_free_surface_pressure:
 * @domain: a #GfsDomain.
 * @par: the multigrid paramaters.
 * @apar: the advection parameters.
 *
 */
static void gfs_free_surface_pressure (GfsDomain * domain,
				       GfsMultilevelParams * par,
				       GfsAdvectionParams * apar,
				       GfsVariable * p,
				       GfsVariable * divn,
				       GfsVariable * res,
				       gdouble G)
{
  GfsDomain * toplayer;
  FreeSurfaceParams fp;
  GfsVariable * res1, * g[2];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (apar != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (divn != NULL);
  g_return_if_fail (G > 0.);

  toplayer = GFS_OCEAN (domain)->toplayer;
  apar->v = gfs_variable_from_name (domain->variables, "U");

  fp.pn = p;
  fp.div = gfs_temporary_variable (domain);
  fp.dia = gfs_temporary_variable (toplayer);
  res1 = res ? res : gfs_temporary_variable (toplayer);
  fp.divn = divn;
  fp.dt = apar->dt;
  fp.G = G/GFS_OCEAN (domain)->layer->len;

  /* Initialize face coefficients */
#if (!FTT_2D) /* 2D3 or 3D */
  gfs_domain_cell_traverse (toplayer,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) sum_coeff, NULL);
  gfs_domain_cell_traverse (toplayer,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) face_coeff_from_below, NULL);
#endif /* 2D3 or 3D */

  /* compute MAC divergence */
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) normal_divergence, &fp);
#if (!FTT_2D)
  gfs_domain_cell_traverse (toplayer, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) sum_divergence, fp.div);
#endif /* 2D3 or 3D */
  gfs_domain_cell_traverse (toplayer, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
  			    (FttCellTraverseFunc) scale_divergence_helmoltz, &fp);
  
  /* solve for pressure */
  par->depth = gfs_domain_depth (toplayer);
  gfs_residual (toplayer, 2, FTT_TRAVERSE_LEAFS, -1, p, fp.div, fp.dia, res1);
  par->residual_before = par->residual = 
    gfs_domain_norm_residual (toplayer, FTT_TRAVERSE_LEAFS, -1, apar->dt, res1);
  par->niter = 0;
  par->dimension = 2;
  while (par->residual.infty > par->tolerance && par->niter < par->nitermax) {
#if 0
    fprintf (stderr, "%d bias: %g first: %g second: %g infty: %g\n",
	     par->niter, 
	     par->residual.bias, 
	     par->residual.first, 
	     par->residual.second, 
	     par->residual.infty);
#endif
    gfs_poisson_cycle (toplayer, par, p, fp.div, fp.dia, res1);
    par->residual = gfs_domain_norm_residual (toplayer, FTT_TRAVERSE_LEAFS, -1, apar->dt, res1);
    par->niter++;
  }

  if (!res)
    gts_object_destroy (GTS_OBJECT (res1));
  gts_object_destroy (GTS_OBJECT (fp.dia));
  gts_object_destroy (GTS_OBJECT (fp.div));

#if (!FTT_2D)
  gfs_domain_cell_traverse (toplayer, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) column_pressure, p);
  gfs_poisson_coefficients (toplayer, NULL, 0.);
#endif /* 2D3 or 3D */
  gfs_correct_normal_velocities_weighted (domain, 2, p, g, apar->dt/2., par->weighted);
#if (!FTT_2D)
  gfs_domain_cell_traverse_boundary (domain, FTT_BACK,
				     FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
  				     (FttCellTraverseFunc) compute_w,
				     gfs_variable_from_name (domain->variables, "W"));
#endif /* 2D3 or 3D */
  gfs_correct_centered_velocities (domain, 2, g, apar->dt/2.);
}

static void gfs_free_surface_divergence (GfsDomain * domain)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (div != NULL);

  gfs_domain_face_traverse (domain, FTT_XY,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
  gfs_domain_face_traverse (domain, FTT_XY,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) gfs_face_interpolated_normal_velocity, 
			    gfs_domain_velocity (domain));
}

static void ocean_run (GfsSimulation * sim)
{
  GfsVariable * p, * div, * res = NULL;
  GfsDomain * domain, * toplayer;
  GSList * i;

  domain = GFS_DOMAIN (sim);
  toplayer = GFS_OCEAN (sim)->toplayer;

  gfs_simulation_refine (sim);

  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_init, sim);
  gts_container_foreach (GTS_CONTAINER (sim->adapts), (GtsFunc) gfs_event_init, sim);

  gfs_set_merged (domain);
  i = domain->variables;
  while (i) {
    gfs_event_init (i->data, sim);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, i->data);
    if (GFS_IS_VARIABLE_RESIDUAL (i->data))
      res = i->data;
    i = i->next;
  }

  p = gfs_variable_from_name (domain->variables, "P");
  g_assert (p);

  div = gfs_temporary_variable (domain);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    GfsVariable * g[2];
    gdouble tstart;

    i = domain->variables;
    while (i) {
      gfs_event_do (GFS_EVENT (i->data), sim);
      i = i->next;
    }
    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    tstart = gfs_clock_elapsed (domain->timer);

    gfs_simulation_set_timestep (sim);

    gfs_free_surface_divergence (domain);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_normal_divergence_2D, div);

    gfs_predicted_face_velocities (domain, 2, &sim->advection_params);

    gfs_domain_timer_start (domain, "correct_normal_velocities");
    gfs_poisson_coefficients (domain, NULL, 0.);
    gfs_correct_normal_velocities_weighted (domain, 2, p, g, sim->advection_params.dt/2.,
					    sim->approx_projection_params.weighted);
#if (!FTT_2D)
    gfs_domain_cell_traverse_boundary (domain, FTT_BACK,
				       FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				       (FttCellTraverseFunc) compute_w, 
				       gfs_variable_from_name (domain->variables, "W"));
#endif /* 2D3 or 3D */
    gfs_domain_timer_stop (domain, "correct_normal_velocities");

    i = domain->variables;
    while (i) {
      if (GFS_IS_VARIABLE_TRACER (i->data)) {
	GfsVariableTracer * t = i->data;

	t->advection.dt = sim->advection_params.dt;
	switch (t->advection.scheme) {
	case GFS_GODUNOV: case GFS_NONE:
	  gfs_tracer_advection_diffusion (domain, &t->advection, &t->diffusion, NULL);
	  break;
	case GFS_VOF:
	  gfs_tracer_vof_advection (domain, &t->advection, NULL);
	  gfs_domain_variable_centered_sources (domain, i->data, i->data, t->advection.dt);
	  break;
	}
      }
      i = i->next;
    }

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_half_do, sim);

    gfs_centered_velocity_advection_diffusion (domain, 2,
					       &sim->advection_params,
					       &sim->diffusion_params,
					       g);

    gfs_poisson_coefficients (domain, NULL, 0.);
    gfs_correct_normal_velocities_weighted (domain, 2, p, g, 0., sim->approx_projection_params.weighted);
    if (gfs_has_source_coriolis (domain)) {
      gfs_correct_centered_velocities (domain, 2, g, sim->advection_params.dt);
      gfs_source_coriolis_implicit (domain, sim->advection_params.dt);
      gfs_correct_normal_velocities_weighted (domain, 2, p, g, 0., sim->approx_projection_params.weighted);
      gfs_correct_centered_velocities (domain, 2, g, -sim->advection_params.dt/2.);
    }
    else
      gfs_correct_centered_velocities (domain, 2, g, sim->advection_params.dt/2.);

    sim->time.t = sim->tnext;
    sim->time.i++;

    gfs_domain_timer_start (domain, "free_surface_pressure");
    gfs_free_surface_divergence (domain);
    gfs_free_surface_pressure (domain, &sim->approx_projection_params, &sim->advection_params,
			       p, div, res, sim->physical_params.g);
    gfs_domain_timer_stop (domain, "free_surface_pressure");

    gfs_simulation_adapt (sim, NULL);

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);

  gts_object_destroy (GTS_OBJECT (div));
}

static void gfs_ocean_class_init (GfsSimulationClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = ocean_destroy;
  GFS_DOMAIN_CLASS (klass)->post_read = ocean_post_read;
  klass->run = ocean_run;
}

static void gfs_ocean_init (GfsOcean * object)
{
  GFS_SIMULATION (object)->approx_projection_params.weighted = 1;
  object->layer = g_ptr_array_new ();
  new_layer (object);
}

GfsSimulationClass * gfs_ocean_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_ocean_info = {
      "GfsOcean",
      sizeof (GfsOcean),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) gfs_ocean_class_init,
      (GtsObjectInitFunc) gfs_ocean_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()), &gfs_ocean_info);
  }

  return klass;
}

#if (!FTT_2D) /* 2D3 or 3D */
static void hydrostatic_pressure (FttCell * cell, gpointer * data)
{
  GfsVariable * vp = data[0];
  GfsVariable * rho = data[1];
  gdouble * g = data[2];
  gdouble r = GFS_VARIABLE (cell, rho->i), p = (*g)*r/2., r1;
  FttCellFace f;
  
  GFS_VARIABLE (cell, vp->i) = p;
  f.cell = cell;
  f.d = FTT_BACK;
  f.neighbor = ftt_cell_neighbor (f.cell, f.d);
  while (f.neighbor) {
    g_assert (ftt_face_type (&f) == FTT_FINE_FINE);
    r1 = gfs_face_interpolated_value (&f, rho->i);
    //    g_assert (r1 >= r);
    r = r1;
    GFS_VARIABLE (f.neighbor, vp->i) = p = p + (*g)*r;
    f.cell = f.neighbor;
    f.neighbor = ftt_cell_neighbor (f.cell, f.d);
  }
}
#endif /* 2D3 or 3D */

/**
 * gfs_hydrostatic_pressure:
 * @domain: a #GfsDomain.
 * @p: the hydrostatic pressure.
 * @rho: the density.
 * @g: the acceleration.
 *
 * Computes the hydrostatic pressure @p in @domain using the density
 * @rho.
 */
void gfs_hydrostatic_pressure (GfsDomain * domain,
			       GfsVariable * p,
			       GfsVariable * rho,
			       gdouble g)
{
#if FTT_2D
  g_return_if_fail (domain != NULL);
  g_return_if_fail (p != NULL);

  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, p);
#else /* 2D3 or 3D */
  gpointer data[3];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (rho != NULL);
  g_return_if_fail (g >= 0.);

  g /= GFS_OCEAN (domain)->layer->len;
  data[0] = p;
  data[1] = rho;
  data[2] = &g;
  gfs_domain_cell_traverse_boundary (domain, FTT_FRONT,
				     FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				     (FttCellTraverseFunc) hydrostatic_pressure, data);
#endif /* 2D3 or 3D */
}

/* GfsSourceHydrostatic: Object */

static void gfs_source_hydrostatic_destroy (GtsObject * o)
{
  if (GFS_SOURCE_HYDROSTATIC (o)->ph1)
    gts_object_destroy (GTS_OBJECT (GFS_SOURCE_HYDROSTATIC (o)->ph1));

  (* GTS_OBJECT_CLASS (gfs_source_hydrostatic_class ())->parent_class->destroy) (o);
}


static void gfs_source_hydrostatic_read (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  GfsSourceHydrostatic * sh;

  if (GTS_OBJECT_CLASS (gfs_source_hydrostatic_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_hydrostatic_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  sh = GFS_SOURCE_HYDROSTATIC (*o);
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (rho)");
    return;
  }
  sh->rho = gfs_variable_from_name (domain->variables, fp->token->str);
  if (sh->rho == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (ph)");
    return;
  }
  if (!(sh->ph = gfs_variable_from_name (domain->variables, fp->token->str)) &&
      !(sh->ph = gfs_domain_add_variable (domain, fp->token->str))) {
    gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  sh->ph1 = gfs_temporary_variable (domain);
}

static void gfs_source_hydrostatic_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_hydrostatic_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s %s",
	   GFS_SOURCE_HYDROSTATIC (o)->rho->name, 
	   GFS_SOURCE_HYDROSTATIC (o)->ph->name);
}

static gdouble gfs_source_hydrostatic_mac_value (GfsSourceGeneric * s,
						 FttCell * cell,
						 GfsVariable * v)
{
  return - gfs_center_gradient (cell, v->component,
				GFS_SOURCE_HYDROSTATIC (s)->ph1->i)/ftt_cell_size (cell);
}

static gdouble gfs_source_hydrostatic_centered_value (GfsSourceGeneric * s,
						      FttCell * cell,
						      GfsVariable * v)
{
  GfsSourceHydrostatic * b = GFS_SOURCE_HYDROSTATIC (s);

  return - (gfs_center_gradient (cell, v->component, b->ph->i) + 
	    gfs_center_gradient (cell, v->component, b->ph1->i))/(2.*ftt_cell_size (cell));
}

static void copy_ph (FttCell * cell, GfsSourceHydrostatic * s)
{
  GFS_VARIABLE (cell, s->ph1->i) = GFS_VARIABLE (cell, s->ph->i);
}

static gboolean gfs_source_hydrostatic_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class)->event) 
      (event, sim)) {
    GfsSourceHydrostatic * s = GFS_SOURCE_HYDROSTATIC (event);

    if (s->not_first) {
      gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) copy_ph, s);
      gfs_domain_bc (GFS_DOMAIN (sim), FTT_TRAVERSE_LEAFS, -1, s->ph1);
    }
    else {
      gfs_hydrostatic_pressure (GFS_DOMAIN (sim), s->ph1, s->rho, sim->physical_params.g);
      gfs_domain_bc (GFS_DOMAIN (sim), FTT_TRAVERSE_LEAFS, -1, s->ph1);
      s->not_first = TRUE;
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_source_hydrostatic_event_half (GfsEvent * event, GfsSimulation * sim)
{
  GfsSourceHydrostatic * s = GFS_SOURCE_HYDROSTATIC (event);

  gfs_hydrostatic_pressure (GFS_DOMAIN (sim), s->ph, s->rho, sim->physical_params.g);
  gfs_domain_bc (GFS_DOMAIN (sim), FTT_TRAVERSE_LEAFS, -1, s->ph);
}

static void gfs_source_hydrostatic_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gfs_source_hydrostatic_destroy;
  GTS_OBJECT_CLASS (klass)->read = gfs_source_hydrostatic_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_source_hydrostatic_write;

  GFS_EVENT_CLASS (klass)->event = gfs_source_hydrostatic_event;
  GFS_EVENT_CLASS (klass)->event_half = gfs_source_hydrostatic_event_half;

  klass->mac_value = gfs_source_hydrostatic_mac_value;
  klass->centered_value = gfs_source_hydrostatic_centered_value;
}

GfsSourceGenericClass * gfs_source_hydrostatic_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_hydrostatic_info = {
      "GfsSourceHydrostatic",
      sizeof (GfsSourceHydrostatic),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_hydrostatic_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
				  &gfs_source_hydrostatic_info);
  }

  return klass;
}

/* GfsSourceFriction: Object */

static void gfs_source_friction_destroy (GtsObject * o)
{
  FttComponent c;

  for (c = 0; c <  FTT_DIMENSION; c++)
    if (GFS_SOURCE_FRICTION (o)->u[c])
      gts_object_destroy (GTS_OBJECT (GFS_SOURCE_FRICTION (o)->u[c]));

  (* GTS_OBJECT_CLASS (gfs_source_friction_class ())->parent_class->destroy) (o);
}

static void gfs_source_friction_read (GtsObject ** o, GtsFile * fp)
{
  FttComponent c;
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  (* GTS_OBJECT_CLASS (gfs_source_friction_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsVariable h)");
    return;
  }
  GFS_SOURCE_FRICTION (*o)->h = gfs_variable_from_name (domain->variables, fp->token->str);
  if (GFS_SOURCE_FRICTION (*o)->h == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (f)");
    return;
  }
  GFS_SOURCE_FRICTION (*o)->f = atof (fp->token->str);
  gts_file_next_token (fp);

  for (c = 0; c <  FTT_DIMENSION; c++)
    GFS_SOURCE_FRICTION (*o)->u[c] = gfs_temporary_variable (domain);
}

static void gfs_source_friction_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_friction_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s %g", GFS_SOURCE_FRICTION (o)->h->name, GFS_SOURCE_FRICTION (o)->f);
}

static gdouble gfs_source_friction_saved_value (GfsSourceGeneric * s, 
						FttCell * cell, 
						GfsVariable * v)
{
  gdouble H = GFS_VARIABLE (cell, GFS_SOURCE_FRICTION (s)->h->i);

  g_assert (H > 0.);
  return - GFS_SOURCE_FRICTION (s)->f*
    GFS_VARIABLE (cell, GFS_SOURCE_FRICTION (s)->u[v->component]->i)/H;
}

static void save_velocity (FttCell * cell, GfsSourceFriction * s)
{
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++)
    GFS_VARIABLE (cell, s->u[c]->i) = GFS_VARIABLE (cell, GFS_SOURCE_VELOCITY (s)->v[c]->i);
}

static gboolean gfs_source_friction_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class)->event)
      (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) save_velocity, event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_source_friction_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gfs_source_friction_destroy;
  GTS_OBJECT_CLASS (klass)->read = gfs_source_friction_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_source_friction_write;
  GFS_EVENT_CLASS (klass)->event = gfs_source_friction_event;
  klass->mac_value = klass->centered_value = gfs_source_friction_saved_value;
}

GfsSourceGenericClass * gfs_source_friction_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_friction_info = {
      "GfsSourceFriction",
      sizeof (GfsSourceFriction),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_friction_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
				  &gfs_source_friction_info);
  }

  return klass;
}

/* GfsBcFlather: Object */

static void bc_flather_write (GtsObject * o, FILE * fp)
{
  GfsBcFlather * bc = GFS_BC_FLATHER (o);

  (* GTS_OBJECT_CLASS (gfs_bc_flather_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s %s", bc->h->name, bc->p->name);
  if (bc->val)
    gfs_function_write (bc->val, fp);
}

static void set_gradient_boundary (FttCell * cell)
{
  cell->flags |= GFS_FLAG_GRADIENT_BOUNDARY;
}

static void bc_flather_read (GtsObject ** o, GtsFile * fp)
{
  GfsBcFlather * bc = GFS_BC_FLATHER (*o);
  GfsDomain * domain = gfs_box_domain (GFS_BC (bc)->b->box);

  (* GTS_OBJECT_CLASS (gfs_bc_flather_class ())->parent_class->read) (o, fp);

  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (h)");
    return;
  }
  bc->h = gfs_variable_from_name (domain->variables, fp->token->str);
  if (bc->h == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (p)");
    return;
  }
  bc->p = gfs_variable_from_name (domain->variables, fp->token->str);
  if (bc->p == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  if (bc->val == NULL)
    bc->val = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_read (bc->val, gfs_box_domain (GFS_BC (bc)->b->box), fp);

  ftt_cell_traverse (GFS_BC (bc)->b->root, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		     (FttCellTraverseFunc) set_gradient_boundary, NULL);
}

static void bc_flather_destroy (GtsObject * o)
{
  if (GFS_BC_FLATHER (o)->val)
    gts_object_destroy (GTS_OBJECT (GFS_BC_FLATHER (o)->val));

  (* GTS_OBJECT_CLASS (gfs_bc_flather_class ())->parent_class->destroy) (o);
}

static gdouble flather_value (FttCellFace * f, GfsBc * b)
{
  guint d, nb = 0;
  FttCellNeighbors n;
  gdouble H;

  ftt_cell_neighbors (f->neighbor, &n);
  for (d = 0; d < FTT_NEIGHBORS_2D; d++)
    if (n.c[d] != NULL && GFS_CELL_IS_BOUNDARY(n.c[d]) && nb++ > 0)
      /* if the boundary cell is bounded by more than one boundary -> no flux */
      return 0.;

  H = gfs_face_interpolated_value (f, GFS_BC_FLATHER (b)->h->i);
  if (H > 2e-3) { /* fixme: 2e-3 is an arbitrary constant which should be a parameter or sthg*/
    GfsSimulation * sim = GFS_SIMULATION (gfs_box_domain (b->b->box));
    gdouble cg = sqrt (sim->physical_params.g*H);
    
    return gfs_function_face_value (GFS_BC_VALUE (b)->val, f) +
      (FTT_FACE_DIRECT (f) ? -1. : 1.)*
      (GFS_VARIABLE (f->neighbor, GFS_BC_FLATHER (b)->p->i) - 
       gfs_function_face_value (GFS_BC_FLATHER (b)->val, f))*
      cg/sim->physical_params.g;
  }
  else
    return 0.;
}

static void flather (FttCellFace * f, GfsBc * b)
{
  g_assert (GFS_CELL_IS_GRADIENT_BOUNDARY (f->cell));
  GFS_VARIABLE (f->cell, b->v->i) = 2.*flather_value (f, b) - GFS_VARIABLE (f->neighbor, b->v->i);
}

static void homogeneous_flather (FttCellFace * f, GfsBc * b)
{
  g_assert (GFS_CELL_IS_GRADIENT_BOUNDARY (f->cell));
  GFS_VARIABLE (f->cell, b->v->i) = - GFS_VARIABLE (f->neighbor, b->v->i);
}

static void face_flather (FttCellFace * f, GfsBc * b)
{
  g_assert (GFS_CELL_IS_GRADIENT_BOUNDARY (f->cell));
  GFS_STATE (f->cell)->f[f->d].v = flather_value (f, b);
}

static void gfs_bc_flather_class_init (GtsObjectClass * klass)
{
  klass->write   = bc_flather_write;
  klass->read    = bc_flather_read;
  klass->destroy = bc_flather_destroy;
}

static void gfs_bc_flather_init (GfsBc * object)
{
  object->bc =             (FttFaceTraverseFunc) flather;
  object->homogeneous_bc = (FttFaceTraverseFunc) homogeneous_flather;
  object->face_bc =        (FttFaceTraverseFunc) face_flather;
}

GfsBcClass * gfs_bc_flather_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_flather_info = {
      "GfsBcFlather",
      sizeof (GfsBcFlather),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_flather_class_init,
      (GtsObjectInitFunc) gfs_bc_flather_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_value_class ()),
				  &gfs_bc_flather_info);
  }

  return klass;
}

/* GfsOcean1: Object */

static void ocean1_run (GfsSimulation * sim)
{
  GfsVariable * p, * div, * H, * res = NULL;
  GfsDomain * domain, * toplayer;
  GSList * i;

  domain = GFS_DOMAIN (sim);
  toplayer = GFS_OCEAN (sim)->toplayer;

  gfs_simulation_refine (sim);

  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_init, sim);
  gts_container_foreach (GTS_CONTAINER (sim->adapts), (GtsFunc) gfs_event_init, sim);

  gfs_set_merged (domain);
  i = domain->variables;
  while (i) {
    gfs_event_init (i->data, sim);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, i->data);
    if (GFS_IS_VARIABLE_RESIDUAL (i->data))
      res = i->data;
    i = i->next;
  }

  p = gfs_variable_from_name (domain->variables, "P");
  g_assert (p);
  H = gfs_variable_from_name (domain->variables, "H");
  g_assert (H);

  div = gfs_temporary_variable (domain);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    GfsVariable * g[2];
    gdouble tstart;

    i = domain->variables;
    while (i) {
      gfs_event_do (i->data, sim);
      i = i->next;
    }
    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    tstart = gfs_clock_elapsed (domain->timer);

    gfs_simulation_set_timestep (sim);

    gfs_free_surface_divergence (domain);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_normal_divergence_2D, div);

    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, p);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_half_do, sim);

    gfs_correct_normal_velocities_weighted (domain, 2, p, g, 0., FALSE); 
    // just there so that the call below 
    // has sthg to free
    gfs_centered_velocity_advection_diffusion (domain, 2,
					       &sim->advection_params,
					       &sim->diffusion_params,
					       g);

    gfs_poisson_coefficients (domain, H, 0.);
    gfs_correct_normal_velocities_weighted (domain, 2, p, g, 0., sim->approx_projection_params.weighted);
    if (gfs_has_source_coriolis (domain)) {
      gfs_correct_centered_velocities (domain, 2, g, sim->advection_params.dt);
      gfs_source_coriolis_implicit (domain, sim->advection_params.dt);
      gfs_correct_normal_velocities_weighted (domain, 2, p, g, 0., sim->approx_projection_params.weighted);
      gfs_correct_centered_velocities (domain, 2, g, -sim->advection_params.dt/2.);
    }
    else
      gfs_correct_centered_velocities (domain, 2, g, sim->advection_params.dt/2.);

    sim->time.t = sim->tnext;
    sim->time.i++;

    gfs_domain_timer_start (domain, "free_surface_pressure");
    gfs_free_surface_divergence (domain);
    gfs_free_surface_pressure (domain, &sim->approx_projection_params, &sim->advection_params,
			       p, div, res, sim->physical_params.g);
    gfs_domain_timer_stop (domain, "free_surface_pressure");

    gfs_simulation_adapt (sim, NULL);

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events),
			 (GtsFunc) gts_object_destroy, NULL);

  gts_object_destroy (GTS_OBJECT (div));
}

static void gfs_ocean1_class_init (GfsSimulationClass * klass)
{
  klass->run = ocean1_run;
}

GfsSimulationClass * gfs_ocean1_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_ocean_info = {
      "GfsOcean1",
      sizeof (GfsOcean),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) gfs_ocean1_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_ocean_class ()), &gfs_ocean_info);
  }

  return klass;
}
