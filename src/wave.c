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

#define NK 25
#define NTHETA 24

static double frequency (int ik) 
{
  double gamma = 1.1;
  double f0 = 0.04;
  return f0*pow(gamma, ik);
}
      
static double gaussian (double f, double fmean, double fsigma) 
{
  return exp (-((f - fmean)*(f - fmean))/(fsigma*fsigma));
}

static double costheta (double theta, double thetam, double thetapower) 
{
  if (fabs (theta - thetam) > M_PI/2.) return 0.;
  double a = cos (theta - thetam);
  return pow (a, thetapower);
}
      
static double theta (int ith) 
{
  return 2.*M_PI*ith/NTHETA;
}

static void cg (int ik, int ith, FttVector * u) 
{
  double cg = 9.81/(4.*M_PI*frequency (ik))/1000./5000.*3600.;
  u->x = cg*cos (theta (ith));
  u->y = cg*sin (theta (ith));
  u->z = 0.;
}

static double action (int ik, int ith, double x, double y, double amp) 
{
  double xc = -0.5 + 500./5000.;
  double yc = -0.5 + 500./5000.;
  x -= xc;
  y -= yc;
  return amp*
    gaussian (frequency (ik), 0.1, 0.01)*
    costheta (theta (ith), 30.*M_PI/180., 2.)*
    gaussian (sqrt (x*x + y*y), 0., 150./5000.);
}

static void init_action (FttCell * cell, GfsVariable *** F)
{
  guint ik, ith;
  FttVector p;
  gfs_cell_cm (cell, &p);
  for (ik = 0; ik < NK; ik++)
    for (ith = 0; ith < NTHETA; ith++)
      GFS_VALUE (cell, F[ik][ith]) = action (ik, ith, p.x, p.y, 1.);
}

static void set_group_velocity (const FttCellFace * face, FttVector * u)
{
  GFS_FACE_NORMAL_VELOCITY_RIGHT (face) = 
    GFS_FACE_NORMAL_VELOCITY_LEFT (face) = (&u->x)[face->d/2];
}

static void wave_run (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);

  gfs_simulation_refine (sim);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) init_action, GFS_WAVE (sim)->F);
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
    for (ik = 0; ik < NK; ik++) {
      FttVector u;
      cg (ik, 0, &u);
      gfs_domain_face_traverse (domain, FTT_XYZ,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttFaceTraverseFunc) set_group_velocity, &u);
      gfs_simulation_set_timestep (sim);
      /* subcycling */
      guint n = ceil (dt/sim->advection_params.dt);
      g_assert (fabs (sim->time.t + sim->advection_params.dt*n - tnext) < 1e-12);
      while (n--) {
	for (ith = 0; ith < NTHETA; ith++) {
	  FttVector u;
	  cg (ik, ith, &u);
	  gfs_domain_face_traverse (domain, FTT_XYZ,
				    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				    (FttFaceTraverseFunc) set_group_velocity, &u);
	  GfsVariable * t = GFS_WAVE (sim)->F[ik][ith];
	  GFS_VARIABLE_TRACER (t)->advection.dt = sim->advection_params.dt;
	  gfs_tracer_advection_diffusion (domain, &GFS_VARIABLE_TRACER (t)->advection);
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
  gfs_matrix_free (GFS_WAVE (object)->F);
  (* GTS_OBJECT_CLASS (gfs_wave_class ())->parent_class->destroy) (object);
}

static void gfs_wave_class_init (GfsSimulationClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = wave_destroy;
  klass->run = wave_run;
}

static void wave_init (GfsWave * wave)
{
  guint ik, ith;
  wave->F = gfs_matrix_new (NK, NTHETA, sizeof (GfsVariable *));
  for (ik = 0; ik < NK; ik++)
    for (ith = 0; ith < NTHETA; ith++) {
      gchar * name = g_strdup_printf ("F%d_%d", ik, ith);
      gchar * description = g_strdup_printf ("Action density for f = %g Hz and theta = %g degrees",
					     frequency (ik), theta (ith)*180./M_PI);
      wave->F[ik][ith] = gfs_domain_add_variable (GFS_DOMAIN (wave), 
						  gfs_variable_tracer_class (), 
						  name, description);
      g_assert (wave->F[ik][ith]);
      GFS_VARIABLE_TRACER (wave->F[ik][ith])->advection.use_centered_velocity = FALSE;
      g_free (name);
      g_free (description);
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
