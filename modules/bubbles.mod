/* Gerris - The GNU Flow Solver			(-*-C-*-)
 * Copyright (C) 2010 National Institute of Water and Atmospheric Research
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

#include "simulation.h"
#include "source.h"
#include "adaptive.h"
#include "output.h"
#include "solid.h"
#include "particulates.h"
#include "ftt.h"

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>


/* GfsBubble: Header */

typedef struct _GfsBubble GfsBubble;

struct _GfsBubble {
  /*< private >*/
  GfsParticulate parent;
  
  /*< public >*/
  gdouble velR, p0, R0;
};

#define GFS_BUBBLE(obj)            GTS_OBJECT_CAST (obj, GfsBubble, gfs_bubble_class ())
#define GFS_IS_BUBBLE(obj)         (gts_object_is_from_class (obj, gfs_bubble_class ()))

static GfsEventClass * gfs_bubble_class  (void);

/* GfsBubble: Object */
/* The radius of each bubble varies according to the Rayleigh-Plesset equation */

typedef struct {
  gdouble liqpres, liqdens;
  GfsBubble * bubble;
} RPData;

int static func (double t, const double y[], double f[], void *params)
{
  RPData *rp = (RPData *)params;
  f[0] = y[1];
  /* interface acceleration- incompressible RP equation */
  gdouble pbubble = rp->bubble->p0*pow (rp->bubble->R0/y[0], 3.*1.4);
  f[1] = ((pbubble - rp->liqpres)/rp->liqdens - 3./2.*y[1]*y[1])/y[0];
  return GSL_SUCCESS;
}

/* jacobian matrix */
int static jac (double t, const double y[], double *dfdy, 
		double dfdt[], void *params)
{
  RPData *rp = (RPData *)params;
  gdouble pbubble = rp->bubble->p0*pow (rp->bubble->R0/y[0], 3.*1.4);
  gdouble dddRdR  = 2.*rp->liqpres-2.*(1. + 3.*1.4)*pbubble + 3.*rp->liqdens*y[1]*y[1];
  dddRdR  /= 2.*y[0]*y[0]*rp->liqdens;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix; 
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, dddRdR);
  gsl_matrix_set (m, 1, 1, - 3.*y[1]/y[0]);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

static gboolean gfs_bubble_event (GfsEvent * event, 
				  GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class)->event) 
      (event, sim)) {
    GfsParticle * p = GFS_PARTICLE (event);
    GfsParticulate * particulate = GFS_PARTICULATE (event);
    GfsBubble * bubble = GFS_BUBBLE (event);
    GfsDomain * domain = GFS_DOMAIN (sim);

    GfsVariable * liqpres = gfs_variable_from_name (domain->variables, "P");
  
    FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
    if (cell == NULL) 
      return TRUE;
    gdouble liq_rho = sim->physical_params.alpha ? 1./
      gfs_function_value (sim->physical_params.alpha, cell) : 1.;

    FttVector pos = p->pos;
    gfs_simulation_map (sim, &pos);

    gdouble point_pres = gfs_interpolate (cell, p->pos, liqpres);

    RPData rp = { point_pres, liq_rho, bubble };

    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;

    gsl_odeiv_step * s    = gsl_odeiv_step_alloc (T, 2);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-6, 0.0);
    gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (2);

    gsl_odeiv_system sys = {func, jac, 2, &rp};

    gdouble t = sim->time.t;
    gdouble t1 = t + sim->advection_params.dt;
    gdouble h = 1e-6; /* better criterion?? */
    /* variables R, dot{R} */
    gdouble y[2] = { pow(3./(4.*M_PI)*particulate->volume,1./3.) , bubble->velR };

    while (t < t1) {
      int status = gsl_odeiv_evolve_apply (e, c, s,
					   &sys, &t, t1, &h, y);
      if (status != GSL_SUCCESS) 
	g_error ("Error in the RK method");
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);

    bubble->velR = y[1];
    particulate->volume = 4./3.*M_PI*y[0]*y[0]*y[0];

    return TRUE;
  }
  return FALSE;
} 

static void gfs_bubble_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  GfsBubble * p = GFS_BUBBLE (*o);
  GfsParticulate * part = GFS_PARTICULATE (*o);
  gdouble L = gfs_object_simulation (*o)->physical_params.L;

  p->R0 = pow (part->volume*3./(4.*M_PI), 1./3.);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (radial velocity)");
    return;
  }
  p->velR = atof (fp->token->str)/L;
  gts_file_next_token (fp);
  
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (reference pressure)");
    return;
  }
  p->p0 = atof (fp->token->str)*L;
  gts_file_next_token (fp);
}

static void gfs_bubble_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class->write) (o, fp); 
  GfsBubble * p = GFS_BUBBLE (o);
  gdouble L = gfs_object_simulation (o)->physical_params.L;
  fprintf (fp, " %g %g", p->velR*L, p->p0/L);
}

static void gfs_bubble_class_init (GfsEventClass * klass)
{
  klass->event = gfs_bubble_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_bubble_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_bubble_write;
}


static GfsEventClass * gfs_bubble_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bubble_info = {
      "GfsBubble",
      sizeof (GfsBubble),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_bubble_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particulate_class ()),
				  &gfs_bubble_info);
  }
  return klass;
}

/* Initialize module */

const gchar gfs_module_name[] = "bubbles";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
  gfs_bubble_class ();
  return NULL; 
}
