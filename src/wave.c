/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2008 National Institute of Water and Atmospheric Research
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

#include "wave.h"
#include "adaptive.h"
#include "solid.h"

/* GfsWave: Object */

#define LENGTH 5000. /* box length: 5000 km */

static double frequency (int ik)
{
  double gamma = 1.1;
  double f0 = 0.04;
  return f0*pow(gamma, ik);
}
      
static double theta (guint ith, guint ntheta)
{
  return 2.*M_PI*ith/ntheta;
}

static void cg (int ik, int ith, FttVector * u, guint ntheta)
{
  double cg = 9.81/(4.*M_PI*frequency (ik))/1000./LENGTH*3600.;
  u->x = cg*cos (theta (ith, ntheta));
  u->y = cg*sin (theta (ith, ntheta));
  u->z = 0.;
}

static gdouble cell_E (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  GfsWave * wave = GFS_WAVE (domain);
  GfsVariable *** F = wave->F;
  guint ik, ith;
  gdouble E = 0.;
  for (ik = 0; ik < wave->nk - 1; ik++) {
    gdouble df = (frequency (ik + 1) - frequency (ik))/2.;
    for (ith = 0; ith < wave->ntheta; ith++)
      E += (GFS_VALUE (cell, F[ik + 1][ith]) + GFS_VALUE (cell, F[ik][ith]))*df;
  }
  return E*2.*M_PI/wave->ntheta;
}

static void set_group_velocity (const FttCellFace * face, FttVector * u)
{
  GFS_FACE_NORMAL_VELOCITY_RIGHT (face) = 
    GFS_FACE_NORMAL_VELOCITY_LEFT (face) = (&u->x)[face->d/2];
}

static void wave_run (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsWave * wave = GFS_WAVE (sim);

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    /* get global timestep */
    gfs_domain_face_traverse (domain, FTT_XYZ,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
    gfs_simulation_set_timestep (sim);
    gdouble dt = sim->advection_params.dt;
    gdouble tnext = sim->tnext;
    
    /* spatial advection */
    guint ik, ith;
    for (ik = 0; ik < wave->nk; ik++) {
      FttVector u;
      cg (ik, 0, &u, wave->ntheta);
      gfs_domain_face_traverse (domain, FTT_XYZ,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttFaceTraverseFunc) set_group_velocity, &u);
      gfs_simulation_set_timestep (sim);
      /* subcycling */
      guint n = rint (dt/sim->advection_params.dt);
      g_assert (fabs (sim->time.t + sim->advection_params.dt*n - tnext) < 1e-12);
      while (n--) {
	for (ith = 0; ith < wave->ntheta; ith++) {
	  FttVector u;
	  cg (ik, ith, &u, wave->ntheta);
	  gfs_domain_face_traverse (domain, FTT_XYZ,
				    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				    (FttFaceTraverseFunc) set_group_velocity, &u);
	  GfsVariable * t = GFS_WAVE (sim)->F[ik][ith];
	  sim->advection_params.v = t;
	  gfs_tracer_advection_diffusion (domain, &sim->advection_params);
	  gfs_domain_cell_traverse (domain,
				    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
				    (FttCellTraverseFunc) t->fine_coarse, t);
	}
	gts_container_foreach (GTS_CONTAINER (sim->adapts), (GtsFunc) gfs_event_redo, sim);
	gfs_simulation_adapt (sim);
      }
    }

    sim->advection_params.dt = dt;
    sim->time.t = sim->tnext = tnext;
    sim->time.i++;

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);
}

static void wave_destroy (GtsObject * object)
{
  if (GFS_WAVE (object)->F)
    gfs_matrix_free (GFS_WAVE (object)->F);
  (* GTS_OBJECT_CLASS (gfs_wave_class ())->parent_class->destroy) (object);
}

static void wave_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_wave_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsWave * wave = GFS_WAVE (*o);
  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_UINT, "nk",     TRUE},
      {GTS_UINT, "ntheta", TRUE},
      {GTS_NONE}
    };
    var[0].data = &wave->nk;
    var[1].data = &wave->ntheta;
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
  }

  GfsDomain * domain = GFS_DOMAIN (wave);
  guint ik, ith;
  wave->F = gfs_matrix_new (wave->nk, wave->ntheta, sizeof (GfsVariable *));
  for (ik = 0; ik < wave->nk; ik++)
    for (ith = 0; ith < wave->ntheta; ith++) {
      gchar * name = g_strdup_printf ("F%d_%d", ik, ith);
      gchar * description = g_strdup_printf ("Action density for f = %g Hz and theta = %g degrees",
					     frequency (ik), theta (ith, wave->ntheta)*180./M_PI);
      wave->F[ik][ith] = gfs_domain_get_or_add_variable (domain, name, description);
      g_assert (wave->F[ik][ith]);
      g_free (name);
      g_free (description);
    }
}

static void wave_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_wave_class ())->parent_class->write) (o, fp);

  GfsWave * wave = GFS_WAVE (o);
  fprintf (fp, " {\n"
	   "  nk = %d\n"
	   "  ntheta = %d\n"
	   "}",
	   wave->nk, wave->ntheta);
}

static void gfs_wave_class_init (GfsSimulationClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = wave_destroy;
  GTS_OBJECT_CLASS (klass)->read = wave_read;
  GTS_OBJECT_CLASS (klass)->write = wave_write;
  klass->run = wave_run;
}

static gdouble cell_hs (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  gdouble E = cell_E (cell, face, domain);
  return E > 0. ? 4.*sqrt (E) : 0.;
}

static gdouble cell_frequency (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  return frequency (GFS_WAVE (domain)->ik);
}

static gdouble cell_direction (FttCell * cell, FttCellFace * face, GfsDomain * domain)
{
  return theta (GFS_WAVE (domain)->ith, GFS_WAVE (domain)->ntheta);
}

static void wave_init (GfsWave * wave)
{
  wave->nk = 25;
  wave->ntheta = 24;

  GfsAdvectionParams * par = &GFS_SIMULATION (wave)->advection_params;
  par->gradient = gfs_center_van_leer_gradient;
  par->flux = gfs_face_advection_flux;
  par->use_centered_velocity = FALSE;  

  static GfsDerivedVariableInfo derived_variable[] = {
    { "Hs", "Significant wave height", cell_hs },
    { "Frequency", "Wave frequency", cell_frequency },
    { "Direction", "Wave direction (angle)", cell_direction },
    { NULL, NULL, NULL}
  };
  GfsDerivedVariableInfo * v = derived_variable;
  while (v->name) {
    g_assert (gfs_domain_add_derived_variable (GFS_DOMAIN (wave), *v));
    v++;
  }
}

GfsSimulationClass * gfs_wave_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_wave_info = {
      "GfsWave",
      sizeof (GfsWave),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) gfs_wave_class_init,
      (GtsObjectInitFunc) wave_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()), &gfs_wave_info);
  }

  return klass;
}

/* GfsInitWave: Object */

static void gfs_init_wave_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_wave_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!GFS_IS_WAVE (domain)) {
    gts_file_error (fp, "GfsInitWave can only be used within a GfsWave simulation");
    return;
  }
  
  gfs_function_read (GFS_INIT_WAVE (*o)->d, domain, fp);
  if (fp->type == GTS_ERROR)
    return;
  gfs_function_read (GFS_INIT_WAVE (*o)->hs, domain, fp);
}

static void gfs_init_wave_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_wave_class ())->parent_class->write) (o, fp);

  gfs_function_write (GFS_INIT_WAVE (o)->d, fp);
  gfs_function_write (GFS_INIT_WAVE (o)->hs, fp);
}

static void gfs_init_wave_destroy (GtsObject * object)
{
  gts_object_destroy (GTS_OBJECT (GFS_INIT_WAVE (object)->d));
  gts_object_destroy (GTS_OBJECT (GFS_INIT_WAVE (object)->hs));

  (* GTS_OBJECT_CLASS (gfs_init_wave_class ())->parent_class->destroy) (object);
}

static void init_action (FttCell * cell, GfsInitWave * event)
{
  GfsWave * wave = GFS_WAVE (gfs_object_simulation (event));
  for (wave->ik = 0; wave->ik < wave->nk; wave->ik++)
    for (wave->ith = 0; wave->ith < wave->ntheta; wave->ith++)
      GFS_VALUE (cell, wave->F[wave->ik][wave->ith]) = gfs_function_value (event->d, cell);

  gdouble E = cell_E (cell, NULL, GFS_DOMAIN (wave));
  gdouble Hs = gfs_function_value (event->hs, cell);
  gdouble scaling = Hs*Hs/(16.*E);
  guint ik, ith;
  for (ik = 0; ik < wave->nk; ik++)
    for (ith = 0; ith < wave->ntheta; ith++)
      GFS_VALUE (cell, wave->F[ik][ith]) *= scaling;
}

static gboolean gfs_init_wave_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_wave_class ())->parent_class)->event) 
      (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_action, event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_wave_class_init (GfsGenericInitClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_wave_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_wave_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_wave_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_init_wave_destroy;
}

static void gfs_init_wave_init (GfsInitWave * object)
{
  object->d = gfs_function_new (gfs_function_class (), 0.);
  object->hs = gfs_function_new (gfs_function_class (), 0.);
}

GfsGenericInitClass * gfs_init_wave_class (void)
{
  static GfsGenericInitClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_wave_info = {
      "GfsInitWave",
      sizeof (GfsInitWave),
      sizeof (GfsGenericInitClass),
      (GtsObjectClassInitFunc) gfs_init_wave_class_init,
      (GtsObjectInitFunc) gfs_init_wave_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_wave_info);
  }

  return klass;
}

