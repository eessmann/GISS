/* Gerris - The GNU Flow Solver                       (-*-C-*-)
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

/* fixme: most of this stuff is obsolete (embedded functions can be
   used instead) */

#include <math.h>
#include <stdlib.h>

#include "event.h"
#include "domain.h"
#include "output.h"

/* GfsInitStationary: Header */

typedef struct _GfsInitStationary         GfsInitStationary;
typedef struct _GfsInitStationaryClass    GfsInitStationaryClass;

struct _GfsInitStationary {
  GfsEvent parent;

  gdouble m;
};

#define GFS_INIT_STATIONARY(obj)            GTS_OBJECT_CAST (obj,\
					           GfsInitStationary,\
					           gfs_init_stationary_class ())
#define GFS_INIT_STATIONARY_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						   GfsInitStationaryClass,\
						   gfs_init_stationary_class())
#define GFS_IS_INIT_STATIONARY(obj)         (gts_object_is_from_class (obj,\
						 gfs_init_stationary_class ()))

GfsEvent * gfs_init_stationary_class (void);

/* GfsInitStationary: Object */

static void init_velocity_stationary (FttCell * cell,
				      GfsInitStationary * init)
{
  FttVector pos;

  ftt_cell_pos (cell, &pos);
  GFS_STATE (cell)->u = 
    - cos (2.*init->m*M_PI*pos.x)*sin (2.*init->m*M_PI*pos.y);
  GFS_STATE (cell)->v =   
    sin (2.*init->m*M_PI*pos.x)*cos (2.*init->m*M_PI*pos.y);
}

static gboolean init_stationary_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* gfs_generic_init_class ()->event) (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_velocity_stationary,
			      event);
    return TRUE;
  }
  return FALSE;
}

static void init_stationary_write (GtsObject * object, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_stationary_class ())->parent_class->write) 
    (object, fp);
  fprintf (fp, " %g", GFS_INIT_STATIONARY (object)->m);
}

static void init_stationary_read (GtsObject ** o, GtsFile * fp)
{
  GfsInitStationary * init;

  (* GTS_OBJECT_CLASS (gfs_init_stationary_class ())->parent_class->read) 
    (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  init = GFS_INIT_STATIONARY (*o);
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (m)");
    return;
  }
  init->m = atof (fp->token->str);
  gts_file_next_token (fp);
}

static void gfs_init_stationary_class_init (GfsEventClass * klass)
{
  klass->event = init_stationary_event;

  GTS_OBJECT_CLASS (klass)->write = init_stationary_write;
  GTS_OBJECT_CLASS (klass)->read = init_stationary_read;
}

static void gfs_init_stationary_init (GfsInitStationary * object)
{
  object->m = 1.;
}

GfsEvent * gfs_init_stationary_class (void)
{
  static GfsEvent * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_stationary_info = {
      "GfsInitStationary",
      sizeof (GfsInitStationary),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_init_stationary_class_init,
      (GtsObjectInitFunc) gfs_init_stationary_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_stationary_info);
  }

  return klass;
}

/* GfsInitNonStationary: Header */

GfsEventClass * gfs_init_non_stationary_class  (void);

/* GfsInitNonStationary: Object */

static void init_velocity_non_stationary (FttCell * cell)
{
  FttVector pos;

  ftt_cell_pos (cell, &pos);
  GFS_STATE (cell)->u = 1. - 2.*cos (2.*M_PI*pos.x)*sin (2.*M_PI*pos.y);
  GFS_STATE (cell)->v = 1. + 2.*sin (2.*M_PI*pos.x)*cos (2.*M_PI*pos.y);
}

static gboolean init_non_stationary_event (GfsEvent * event, 
					   GfsSimulation * sim)
{
  if ((* gfs_generic_init_class ()->event) (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			     FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
       (FttCellTraverseFunc) init_velocity_non_stationary,
			     event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_non_stationary_class_init (GfsEventClass * klass)
{
  klass->event = init_non_stationary_event;
}

GfsEventClass * gfs_init_non_stationary_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_non_stationary_info = {
      "GfsInitNonStationary",
      sizeof (GfsEvent),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_init_non_stationary_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_non_stationary_info);
  }

  return klass;
}

/* GfsOutputNonStationaryError: Header */

GfsOutputClass * gfs_output_non_stationary_error_class  (void);

/* GfsOutputNonStationaryError: Object */

static void compute_non_stationary_error (FttCell * cell, gdouble * t)
{
  FttVector pos;

  ftt_cell_pos (cell, &pos);
  GFS_STATE (cell)->dp = 
    1. - 2.*cos (2.*M_PI*(pos.x - *t))*sin (2.*M_PI*(pos.y - *t))
    - GFS_STATE (cell)->u;
}

static gboolean non_stationary_error_event (GfsEvent * event, 
					    GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class ())->event) (event, sim)) {
    GfsNorm norm;
    
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			     FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
       (FttCellTraverseFunc) compute_non_stationary_error, 
			     &sim->time.t);
    norm = gfs_domain_norm_variable (GFS_DOMAIN (sim),
				     gfs_dp, FTT_TRAVERSE_LEAFS, -1);
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "domain: %10e %10e %10e %10e %10e %10e\n",
	     sim->time.t,
	     norm.first, norm.second, norm.infty, norm.bias, norm.w);
    norm = gfs_domain_norm_variable (GFS_DOMAIN (sim),
				     gfs_dp, FTT_TRAVERSE_LEVEL,
				     gfs_domain_depth (GFS_DOMAIN (sim)));
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "patch:  %10e %10e %10e %10e %10e %10e\n",
	     sim->time.t,
	     norm.first, norm.second, norm.infty, norm.bias, norm.w);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_non_stationary_error_class_init (GfsEventClass * klass)
{
  klass->event = non_stationary_error_event;
}

GfsOutputClass * gfs_output_non_stationary_error_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_non_stationary_error_info = {
      "GfsOutputNonStationaryError",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_non_stationary_error_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_non_stationary_error_info);
  }

  return klass;
}

/* GfsOutputStationaryError: Header */

typedef struct _GfsOutputStationaryError         GfsOutputStationaryError;

struct _GfsOutputStationaryError {
  GfsOutput parent;

  gdouble m;
};

#define GFS_OUTPUT_STATIONARY_ERROR(obj)            GTS_OBJECT_CAST (obj,\
					           GfsOutputStationaryError,\
					gfs_output_stationary_error_class ())
#define GFS_OUTPUT_STATIONARY_ERROR_CLASS(klass)  GTS_OBJECT_CLASS_CAST (klass,\
						GfsOutputStationaryErrorClass,\
					   gfs_output_stationary_error_class())
#define GFS_IS_OUTPUT_STATIONARY_ERROR(obj)  (gts_object_is_from_class (obj,\
					 gfs_output_stationary_error_class ()))


GfsOutputClass * gfs_output_stationary_error_class  (void);

/* GfsOutputStationaryError: Object */

static void compute_stationary_error (FttCell * cell, gdouble * m)
{
  FttVector pos;

  ftt_cell_pos (cell, &pos);
  GFS_STATE (cell)->dp = - cos (2.*(*m)*M_PI*pos.x)*sin (2.*(*m)*M_PI*pos.y)
    - GFS_STATE (cell)->u;
}

static gboolean stationary_error_event (GfsEvent * event, 
					GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class ())->event) (event, sim)) {
    GfsNorm norm;
    
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) compute_stationary_error, 
			      &GFS_OUTPUT_STATIONARY_ERROR (event)->m);
    norm = gfs_domain_norm_variable (GFS_DOMAIN (sim),
				     gfs_dp, FTT_TRAVERSE_LEAFS, -1);
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "%10e %10e %10e %10e %10e\n",
	     sim->time.t,
	     norm.first, norm.second, norm.infty, norm.bias);
    return TRUE;
  }
  return FALSE;
}

static void output_stationary_error_write (GtsObject * object, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_stationary_error_class ())->parent_class->write) (object, fp);
  fprintf (fp, " %g", GFS_OUTPUT_STATIONARY_ERROR (object)->m);
}

static void output_stationary_error_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputStationaryError * output;

  (* GTS_OBJECT_CLASS (gfs_output_stationary_error_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  output = GFS_OUTPUT_STATIONARY_ERROR (*o);
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (m)");
    return;
  }
  output->m = atof (fp->token->str);
  gts_file_next_token (fp);
}

static void gfs_output_stationary_error_class_init (GfsEventClass * klass)
{
  klass->event = stationary_error_event;

  GTS_OBJECT_CLASS (klass)->write = output_stationary_error_write;
  GTS_OBJECT_CLASS (klass)->read = output_stationary_error_read;
}

static void gfs_output_stationary_error_init (GfsOutputStationaryError * object)
{
  object->m = 1.;
}

GfsOutputClass * gfs_output_stationary_error_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_stationary_error_info = {
      "GfsOutputStationaryError",
      sizeof (GfsOutputStationaryError),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_stationary_error_class_init,
      (GtsObjectInitFunc) gfs_output_stationary_error_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_stationary_error_info);
  }

  return klass;
}

/* GfsInitShear: Header */

typedef struct _GfsInitShear         GfsInitShear;

struct _GfsInitShear {
  /*< private >*/
  GfsGenericInit parent;

  gdouble r;
  GfsVariable * c;
};

typedef struct _GfsInitShearClass    GfsInitShearClass;

struct _GfsInitShearClass {
  /*< private >*/
  GfsGenericInitClass parent_class;
};

#define GFS_INIT_SHEAR(obj)            GTS_OBJECT_CAST (obj,\
					         GfsInitShear,\
					         gfs_init_shear_class ())
#define GFS_INIT_SHEAR_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsInitShearClass,\
						 gfs_init_shear_class())
#define GFS_IS_INIT_SHEAR(obj)         (gts_object_is_from_class (obj,\
						 gfs_init_shear_class ()))

GfsInitShearClass * gfs_init_shear_class  (void);

/* GfsInitShear: Object */

static void gfs_init_shear_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_init_shear_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_shear_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (r)");
    return;
  }
  GFS_INIT_SHEAR (*o)->r = atof (fp->token->str);

  gts_file_next_token (fp);
}

static void gfs_init_shear_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_init_shear_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_init_shear_class ())->parent_class->write) 
      (o, fp);
  fprintf (fp, " %g", GFS_INIT_SHEAR (o)->r);
}

static void init_shear (FttCell * cell, GfsInitShear * s)
{
  FttVector p;

  ftt_cell_pos (cell, &p);
  GFS_STATE (cell)->u = p.y <= 0. ?
    tanh (s->r*(p.y + 0.25)) :
    tanh (s->r*(0.25 - p.y));
  GFS_STATE (cell)->v = 0.05*sin (2.*M_PI*(p.x + 0.5));
  if (s->c)
    GFS_VARIABLE (cell, s->c->i) = p.y <= 0. ?
      0.5 + 0.5*tanh (s->r*(p.y + 0.25)):
      0.5 + 0.5*tanh (s->r*(0.25 - p.y));
}

static gboolean gfs_init_shear_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_shear_class ())->parent_class)->event) (event, sim)) {
    GFS_INIT_SHEAR (event)->c = 
      gfs_variable_from_name (GFS_DOMAIN (sim)->variables, "C");
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_shear,
			      event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_shear_class_init (GfsInitShearClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_shear_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_shear_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_shear_write;
}

static void gfs_init_shear_init (GfsInitShear * object)
{
  object->r = 30.;
}

GfsInitShearClass * gfs_init_shear_class (void)
{
  static GfsInitShearClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_shear_info = {
      "GfsInitShear",
      sizeof (GfsInitShear),
      sizeof (GfsInitShearClass),
      (GtsObjectClassInitFunc) gfs_init_shear_class_init,
      (GtsObjectInitFunc) gfs_init_shear_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_shear_info);
  }

  return klass;
}

/* GfsInitWave: Header */

typedef struct _GfsInitWave         GfsInitWave;

struct _GfsInitWave {
  /*< private >*/
  GfsGenericInit parent;

  /*< public >*/
  gdouble a, h, g, y0, w, k, lambda;
};

#define GFS_INIT_WAVE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsInitWave,\
					         gfs_init_wave_class ())
#define GFS_IS_INIT_WAVE(obj)         (gts_object_is_from_class (obj,\
						 gfs_init_wave_class ()))

GfsGenericInitClass * gfs_init_wave_class  (void);

/* GfsInitWave: Object */

static void gfs_init_wave_read (GtsObject ** o, GtsFile * fp)
{
  GfsInitWave * w;
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "lambda", TRUE},
    {GTS_DOUBLE, "a",      TRUE},
    {GTS_DOUBLE, "h",      TRUE},
    {GTS_DOUBLE, "g",      TRUE},
    {GTS_DOUBLE, "y0",     TRUE},
    {GTS_DOUBLE, "w",      TRUE},
    {GTS_NONE}
  };

  if (GTS_OBJECT_CLASS (gfs_init_wave_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_wave_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  w = GFS_INIT_WAVE (*o);
  var[0].data = &w->lambda;
  var[1].data = &w->a;
  var[2].data = &w->h;
  var[3].data = &w->g;
  var[4].data = &w->y0;
  var[5].data = &w->w;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;
  if (w->g <= 0.) {
    gts_file_variable_error (fp, var, "g", "g must be strictly positive");
    return;
  }
  if (w->h <= 0.) {
    gts_file_variable_error (fp, var, "h", "h must be strictly positive");
    return;
  }
  if (w->lambda <= 0) {
    gts_file_variable_error (fp, var, "lambda", 
			     "lambda must be strictly positive");
    return;
  }
  w->k = 2.*M_PI/w->lambda;
}

static void gfs_init_wave_write (GtsObject * o, FILE * fp)
{
  GfsInitWave * w;

  if (GTS_OBJECT_CLASS (gfs_init_wave_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_init_wave_class ())->parent_class->write) 
      (o, fp);
  w = GFS_INIT_WAVE (o);
  fprintf (fp, " { lambda = %g a = %g h = %g g = %g y0 = %g w = %g }",
	   w->lambda, w->a, w->h, w->g, w->y0, w->w);
}

static void init_wave (FttCell * cell, GfsInitWave * w)
{
  FttVector p;
  gdouble sigma = sqrt (w->k*w->g), etha;
  gdouble u = w->k*w->g*w->a/(sigma*cosh (w->k*w->h));

  ftt_cell_pos (cell, &p);
  p.x += 0.5;
  etha = w->y0 + w->a*sin (w->k*p.x);
#if 0 /* fixme: all this stuff is obsolete anyway */
  GFS_STATE (cell)->c = (1. - tanh (w->w*(p.y - etha)))/2.;
  GFS_STATE (cell)->u = - GFS_STATE (cell)->c*u*
    cosh (w->k*(p.y - w->y0 + w->h))*sin (w->k*p.x);
  GFS_STATE (cell)->v = GFS_STATE (cell)->c*u*
    sinh (w->k*(p.y - w->y0 + w->h))*cos (w->k*p.x);
#endif
}

static gboolean gfs_init_wave_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_wave_class ())->parent_class)->event) (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_wave,
			      event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_wave_class_init (GfsGenericInitClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_wave_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_wave_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_wave_write;
}

static void gfs_init_wave_init (GfsInitWave * w)
{
  w->lambda = 1.;
  w->k = 2.*M_PI;
  w->g = 1.;
  w->h = 0.5;
  w->a = 0.05;
  w->w = 30.;
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

/* Initialize module */

const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_init_stationary_class ();
  gfs_output_stationary_error_class ();
  gfs_init_non_stationary_class ();
  gfs_output_non_stationary_error_class ();
  gfs_init_shear_class ();
  gfs_init_wave_class ();
  return NULL;
}
