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

#include <math.h>
#include <stdlib.h>

#include "refine.h"
#include "event.h"
#include "boundary.h"
#include "output.h"
#include "solid.h"

/* GfsRefineBox: Header */

#define GFS_IS_REFINE_BOX(obj)         (gts_object_is_from_class (obj,\
						 gfs_refine_box_class ()))

GfsRefineClass * gfs_refine_box_class  (void);

/* GfsRefineBox: Object */

static gboolean refine_box_maxlevel (FttCell * cell, guint * maxlevel)
{
  if (ftt_cell_level (cell) < *maxlevel) {
    FttVector pos;

    ftt_cell_pos (cell, &pos);
    if (pos.x < 0.25 && pos.x > -0.25 &&
	pos.y < 0.25 && pos.y > -0.25)
    return TRUE;
  }
  return FALSE;
}

static void refine_box (GfsBox * box, guint * maxlevel)
{
  ftt_cell_refine (box->root, 
		   (FttCellRefineFunc) refine_box_maxlevel, maxlevel,
		   (FttCellInitFunc) gfs_cell_init, gfs_box_domain (box));
}

static void box_refine (GfsRefine * refine, GfsSimulation * sim)
{
  gts_container_foreach (GTS_CONTAINER (sim),
			 (GtsFunc) refine_box,
			 &(refine->maxlevel));
}

static void gfs_refine_box_class_init (GfsRefineClass * klass)
{
  klass->refine = box_refine;
}

GfsRefineClass * gfs_refine_box_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_box_info = {
      "GfsRefineBox",
      sizeof (GfsRefine),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_box_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_class ()),
				  &gfs_refine_box_info);
  }

  return klass;
}

/* GfsRefineSphere: Header */

typedef struct _GfsRefineSphere         GfsRefineSphere;

struct _GfsRefineSphere {
  /*< private >*/
  GfsRefine parent;

  /*< public >*/
  gdouble x, y, z, r;
};

#define GFS_REFINE_SPHERE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsRefineSphere,\
					         gfs_refine_sphere_class ())
#define GFS_IS_REFINE_SPHERE(obj)         (gts_object_is_from_class (obj,\
						 gfs_refine_sphere_class ()))

GfsRefineClass * gfs_refine_sphere_class  (void);

/* GfsRefineSphere: Object */

static gboolean refine_sphere_maxlevel (FttCell * cell, GfsRefine * refine)
{
  if (ftt_cell_level (cell) < refine->maxlevel) {
    FttVector pos;
    GfsRefineSphere * s = GFS_REFINE_SPHERE (refine);

    ftt_cell_pos (cell, &pos);
    if ((pos.x - s->x)*(pos.x - s->x) + 
	(pos.y - s->y)*(pos.y - s->y) + 
	(pos.z - s->z)*(pos.z - s->z) <= s->r*s->r)
      return TRUE;
  }
  return FALSE;
}

static void refine_sphere (GfsBox * box, GfsRefine * refine)
{
  ftt_cell_refine (box->root, 
		   (FttCellRefineFunc) refine_sphere_maxlevel, refine,
		   (FttCellInitFunc) gfs_cell_init, gfs_box_domain (box));
}

static void sphere_refine (GfsRefine * refine, GfsSimulation * sim)
{
  gts_container_foreach (GTS_CONTAINER (sim),
			 (GtsFunc) refine_sphere, refine);
}

static void gfs_refine_sphere_read (GtsObject ** o, GtsFile * fp)
{  
  GfsRefineSphere * s = GFS_REFINE_SPHERE (*o);
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x",    TRUE},
    {GTS_DOUBLE, "y",   TRUE},
    {GTS_DOUBLE, "z",   TRUE},
    {GTS_DOUBLE, "r",   TRUE},
    {GTS_NONE}
  };
  
  if (GTS_OBJECT_CLASS (gfs_refine_sphere_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_refine_sphere_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  var[0].data = &s->x;
  var[1].data = &s->y;
  var[2].data = &s->z;
  var[3].data = &s->r;
  gts_file_assign_variables (fp, var);
}

static void gfs_refine_sphere_write (GtsObject * o, FILE * fp)
{
  GfsRefineSphere * s = GFS_REFINE_SPHERE (o);

  if (GTS_OBJECT_CLASS (gfs_refine_sphere_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_refine_sphere_class ())->parent_class->write) 
      (o, fp);
  fprintf (fp, " { x = %g y = %g z = %g r = %g }", s->x, s->y, s->z, s->r);
}

static void gfs_refine_sphere_class_init (GfsRefineClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_refine_sphere_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_refine_sphere_write;
  klass->refine = sphere_refine;
}

static void gfs_refine_sphere_init (GfsRefineSphere * object)
{
  object->x = object->y = object->z = 0.;
  object->r = 0.5;
}

GfsRefineClass * gfs_refine_sphere_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_sphere_info = {
      "GfsRefineSphere",
      sizeof (GfsRefineSphere),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_sphere_class_init,
      (GtsObjectInitFunc) gfs_refine_sphere_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_class ()),
				  &gfs_refine_sphere_info);
  }

  return klass;
}

/* GfsAddGaussianVortex: Header */

typedef struct _GfsAddGaussianVortex         GfsAddGaussianVortex;

struct _GfsAddGaussianVortex {
  /*< private >*/
  GfsGenericInit parent;

  /*< public >*/
  gdouble x, y, scale;
};

typedef struct _GfsAddGaussianVortexClass    GfsAddGaussianVortexClass;

struct _GfsAddGaussianVortexClass {
  /*< private >*/
  GfsGenericInitClass parent_class;
};

#define GFS_ADD_GAUSSIAN_VORTEX(obj)  GTS_OBJECT_CAST (obj,\
					  GfsAddGaussianVortex,\
					  gfs_add_gaussian_vortex_class ())
#define GFS_ADD_GAUSSIAN_VORTEX_CLASS(klass) GTS_OBJECT_CLASS_CAST (klass,\
					  GfsAddGaussianVortexClass,\
					  gfs_add_gaussian_vortex_class())
#define GFS_IS_ADD_GAUSSIAN_VORTEX(obj) (gts_object_is_from_class (obj,\
					 gfs_add_gaussian_vortex_class ()))

GfsAddGaussianVortexClass * gfs_add_gaussian_vortex_class  (void);

/* GfsAddGaussianVortex: Object */

static void add_gaussian_vortex (FttCell * cell, 
				    GfsAddGaussianVortex * g)
{
  FttVector pos;

  ftt_cell_pos (cell, &pos);
  GFS_STATE (cell)->div += 2.*M_PI*exp (- 2.*((pos.x - g->x)*(pos.x - g->x) + 
					      (pos.y - g->y)*(pos.y - g->y))/
					(g->scale*g->scale));
}

static gboolean gfs_add_gaussian_vortex_event (GfsEvent * event, 
						   GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_add_gaussian_vortex_class ())->parent_class)->event) (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) add_gaussian_vortex,
			      event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_add_gaussian_vortex_read (GtsObject ** o, GtsFile * fp)
{
  GfsAddGaussianVortex * g = GFS_ADD_GAUSSIAN_VORTEX (*o);
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x",    TRUE},
    {GTS_DOUBLE, "y",   TRUE},
    {GTS_DOUBLE, "scale",   TRUE},
    {GTS_NONE}
  };

  if (GTS_OBJECT_CLASS (gfs_add_gaussian_vortex_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_add_gaussian_vortex_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  var[0].data = &g->x;
  var[1].data = &g->y;
  var[2].data = &g->scale;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;
  if (g->scale <= 0.)
    gts_file_variable_error (fp, var, "scale",
			     "scale must be strictly positive");
}

static void gfs_add_gaussian_vortex_write (GtsObject * o, FILE * fp)
{
  GfsAddGaussianVortex * g = GFS_ADD_GAUSSIAN_VORTEX (o);

  if (GTS_OBJECT_CLASS (gfs_add_gaussian_vortex_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_add_gaussian_vortex_class ())->parent_class->write) 
      (o, fp);
  fprintf (fp, " { x = %g y = %g scale = %g }", g->x, g->y, g->scale);
}

static void 
gfs_add_gaussian_vortex_class_init (GfsAddGaussianVortexClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_add_gaussian_vortex_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_add_gaussian_vortex_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_add_gaussian_vortex_write;
}

static void 
gfs_add_gaussian_vortex_init (GfsAddGaussianVortex * object)
{
  object->x = object->y = 0.;
  object->scale = 1.;
}

GfsAddGaussianVortexClass * gfs_add_gaussian_vortex_class (void)
{
  static GfsAddGaussianVortexClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_add_gaussian_vortex_info = {
      "GfsAddGaussianVortex",
      sizeof (GfsAddGaussianVortex),
      sizeof (GfsAddGaussianVortexClass),
      (GtsObjectClassInitFunc) gfs_add_gaussian_vortex_class_init,
      (GtsObjectInitFunc) gfs_add_gaussian_vortex_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_add_gaussian_vortex_info);
  }

  return klass;
}

/* GfsAddAlmgrenVortex: Header */

typedef struct _GfsAddAlmgrenVortex         GfsAddAlmgrenVortex;

struct _GfsAddAlmgrenVortex {
  /*< private >*/
  GfsGenericInit parent;

  /*< public >*/
  gdouble x, y, r1, r2, strength;
};

typedef struct _GfsAddAlmgrenVortexClass    GfsAddAlmgrenVortexClass;

struct _GfsAddAlmgrenVortexClass {
  /*< private >*/
  GfsGenericInitClass parent_class;
};

#define GFS_ADD_ALMGREN_VORTEX(obj)  GTS_OBJECT_CAST (obj,\
					  GfsAddAlmgrenVortex,\
					  gfs_add_almgren_vortex_class ())
#define GFS_ADD_ALMGREN_VORTEX_CLASS(klass) GTS_OBJECT_CLASS_CAST (klass,\
					  GfsAddAlmgrenVortexClass,\
					  gfs_add_almgren_vortex_class())
#define GFS_IS_ADD_ALMGREN_VORTEX(obj) (gts_object_is_from_class (obj,\
					 gfs_add_almgren_vortex_class ()))

GfsAddAlmgrenVortexClass * gfs_add_almgren_vortex_class  (void);

/* GfsAddAlmgrenVortex: Object */

static void add_almgren_vortex (FttCell * cell, 
				GfsAddAlmgrenVortex * g)
{
  FttVector p;
  gdouble r;

  ftt_cell_pos (cell, &p);
  r = sqrt ((p.x - g->x)*(p.x - g->x) + (p.y - g->y)*(p.y - g->y));
  GFS_STATE (cell)->div += g->strength*(1. + tanh ((g->r1 - r)/g->r2))/2.;
}

static gboolean gfs_add_almgren_vortex_event (GfsEvent * event, 
						   GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_add_almgren_vortex_class ())->parent_class)->event) (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) add_almgren_vortex,
			      event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_add_almgren_vortex_read (GtsObject ** o, GtsFile * fp)
{
  GfsAddAlmgrenVortex * g = GFS_ADD_ALMGREN_VORTEX (*o);
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x",        TRUE},
    {GTS_DOUBLE, "y",        TRUE},
    {GTS_DOUBLE, "r1",       TRUE},
    {GTS_DOUBLE, "r2",       TRUE},
    {GTS_DOUBLE, "strength", TRUE},
    {GTS_NONE}
  };

  if (GTS_OBJECT_CLASS (gfs_add_almgren_vortex_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_add_almgren_vortex_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  var[0].data = &g->x;
  var[1].data = &g->y;
  var[2].data = &g->r1;
  var[3].data = &g->r2;
  var[4].data = &g->strength;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;
  if (g->r1 <= 0.)
    gts_file_variable_error (fp, var, "r1",
			     "r1 must be strictly positive");
  else if (g->r2 <= 0.)
    gts_file_variable_error (fp, var, "r2",
			     "r2 must be strictly positive");
}

static void gfs_add_almgren_vortex_write (GtsObject * o, FILE * fp)
{
  GfsAddAlmgrenVortex * g = GFS_ADD_ALMGREN_VORTEX (o);

  if (GTS_OBJECT_CLASS (gfs_add_almgren_vortex_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_add_almgren_vortex_class ())->parent_class->write) 
      (o, fp);
  fprintf (fp, " { x = %g y = %g r1 = %g r2 = %g strength = %g }", 
	   g->x, g->y, g->r1, g->r2, g->strength);
}

static void 
gfs_add_almgren_vortex_class_init (GfsAddAlmgrenVortexClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_add_almgren_vortex_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_add_almgren_vortex_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_add_almgren_vortex_write;
}

static void 
gfs_add_almgren_vortex_init (GfsAddAlmgrenVortex * object)
{
  object->x = object->y = 0.;
  object->r1 = 0.03;
  object->r2 = 0.01;
  object->strength = 1.;
}

GfsAddAlmgrenVortexClass * gfs_add_almgren_vortex_class (void)
{
  static GfsAddAlmgrenVortexClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_add_almgren_vortex_info = {
      "GfsAddAlmgrenVortex",
      sizeof (GfsAddAlmgrenVortex),
      sizeof (GfsAddAlmgrenVortexClass),
      (GtsObjectClassInitFunc) gfs_add_almgren_vortex_class_init,
      (GtsObjectInitFunc) gfs_add_almgren_vortex_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_add_almgren_vortex_info);
  }

  return klass;
}

/* GfsInitRayleigh: Header */

typedef struct _GfsInitRayleigh         GfsInitRayleigh;

struct _GfsInitRayleigh {
  /*< private >*/
  GfsGenericInit parent;

  gdouble r, dx;
  GfsVariable * c;
};

typedef struct _GfsInitRayleighClass    GfsInitRayleighClass;

struct _GfsInitRayleighClass {
  /*< private >*/
  GfsGenericInitClass parent_class;
};

#define GFS_INIT_RAYLEIGH(obj)            GTS_OBJECT_CAST (obj,\
					         GfsInitRayleigh,\
					         gfs_init_rayleigh_class ())
#define GFS_INIT_RAYLEIGH_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsInitRayleighClass,\
						 gfs_init_rayleigh_class())
#define GFS_IS_INIT_RAYLEIGH(obj)         (gts_object_is_from_class (obj,\
						 gfs_init_rayleigh_class ()))

GfsInitRayleighClass * gfs_init_rayleigh_class  (void);

/* GfsInitRayleigh: Object */

static void gfs_init_rayleigh_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_init_rayleigh_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_rayleigh_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (r)");
    return;
  }
  GFS_INIT_RAYLEIGH (*o)->r = atof (fp->token->str);

  gts_file_next_token (fp);
}

static void gfs_init_rayleigh_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_init_rayleigh_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_init_rayleigh_class ())->parent_class->write) 
      (o, fp);

  fprintf (fp, " %g", GFS_INIT_RAYLEIGH (o)->r);
}

static void init_rayleigh (FttCell * cell, GfsInitRayleigh * s)
{
  FttVector p;

  ftt_cell_pos (cell, &p);
#if 1
  GFS_STATE (cell)->v = 0.005*(cos (2.*M_PI*(p.x + s->dx + 0.5)) + 1.);
  if (s->c)
    GFS_VARIABLE (cell, s->c->i) = 0.5*(1. + tanh (s->r*(p.y + 0.5)));
#else
  if (s->c)
    GFS_VARIABLE (cell, s->c->i) = 
      0.5*(1. + tanh (s->r*(p.y + 0.025*sin (2.*M_PI*p.x))/2.));
#endif
}

static gboolean gfs_init_rayleigh_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_rayleigh_class ())->parent_class)->event) (event, sim)) {
    GfsInitRayleigh * rayleigh = GFS_INIT_RAYLEIGH (event);

    rayleigh->c = gfs_variable_from_name (GFS_DOMAIN (sim)->variables, "C");
    /* shift origin by half a cell to ensure symmetry of solution */
    rayleigh->dx = ftt_level_size (gfs_domain_depth (GFS_DOMAIN (sim)))/2.;
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_rayleigh,
			      rayleigh);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_rayleigh_class_init (GfsInitRayleighClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_rayleigh_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_rayleigh_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_rayleigh_write;
}

static void gfs_init_rayleigh_init (GfsInitRayleigh * object)
{
  object->r = 30.;
}

GfsInitRayleighClass * gfs_init_rayleigh_class (void)
{
  static GfsInitRayleighClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_rayleigh_info = {
      "GfsInitRayleigh",
      sizeof (GfsInitRayleigh),
      sizeof (GfsInitRayleighClass),
      (GtsObjectClassInitFunc) gfs_init_rayleigh_class_init,
      (GtsObjectInitFunc) gfs_init_rayleigh_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_rayleigh_info);
  }

  return klass;
}

/* GfsInitBubble: Header */

typedef struct _GfsInitBubble         GfsInitBubble;

struct _GfsInitBubble {
  /*< private >*/
  GfsGenericInit parent;

  gdouble x, y, z, r, t;
  GfsVariable * c;
};

typedef struct _GfsInitBubbleClass    GfsInitBubbleClass;

struct _GfsInitBubbleClass {
  /*< private >*/
  GfsGenericInitClass parent_class;
};

#define GFS_INIT_BUBBLE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsInitBubble,\
					         gfs_init_bubble_class ())
#define GFS_INIT_BUBBLE_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsInitBubbleClass,\
						 gfs_init_bubble_class())
#define GFS_IS_INIT_BUBBLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_init_bubble_class ()))

GfsInitBubbleClass * gfs_init_bubble_class  (void);

/* GfsInitBubble: Object */

static void gfs_init_bubble_read (GtsObject ** o, GtsFile * fp)
{
  GfsInitBubble * b = GFS_INIT_BUBBLE (*o);
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x",        TRUE},
    {GTS_DOUBLE, "y",        TRUE},
    {GTS_DOUBLE, "z",        TRUE},
    {GTS_DOUBLE, "r",        TRUE},
    {GTS_DOUBLE, "t",        TRUE},
    {GTS_NONE}
  };

  if (GTS_OBJECT_CLASS (gfs_init_bubble_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_bubble_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  var[0].data = &b->x;
  var[1].data = &b->y;
  var[2].data = &b->z;
  var[3].data = &b->r;
  var[4].data = &b->t;
  gts_file_assign_variables (fp, var);
}

static void gfs_init_bubble_write (GtsObject * o, FILE * fp)
{
  GfsInitBubble * b = GFS_INIT_BUBBLE (o);
  if (GTS_OBJECT_CLASS (gfs_init_bubble_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_init_bubble_class ())->parent_class->write) 
      (o, fp);

  fprintf (fp, " { x = %g y = %g z = %g r = %g t = %g }", 
	   b->x, b->y, b->z, b->r, b->t);
}

static void init_bubble (FttCell * cell, GfsInitBubble * b)
{
  FttVector p;

  ftt_cell_pos (cell, &p);
  GFS_VARIABLE (cell, b->c->i) = 
    (1. + tanh (b->t*(sqrt ((p.x - b->x)*(p.x - b->x) + 
			    (p.y - b->y)*(p.y - b->y) + 
			    (p.z - b->z)*(p.z - b->z)) - b->r)))/2.;
}

static gboolean gfs_init_bubble_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_bubble_class ())->parent_class)->event) (event, sim)) {
    GfsInitBubble * bubble = GFS_INIT_BUBBLE (event);

    bubble->c = gfs_variable_from_name (GFS_DOMAIN (sim)->variables, "C");
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_bubble,
			      bubble);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_bubble_class_init (GfsInitBubbleClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_bubble_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_bubble_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_bubble_write;
}

static void gfs_init_bubble_init (GfsInitBubble * object)
{
  object->r = 0.25;
  object->t = 50.;
}

GfsInitBubbleClass * gfs_init_bubble_class (void)
{
  static GfsInitBubbleClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_bubble_info = {
      "GfsInitBubble",
      sizeof (GfsInitBubble),
      sizeof (GfsInitBubbleClass),
      (GtsObjectClassInitFunc) gfs_init_bubble_class_init,
      (GtsObjectInitFunc) gfs_init_bubble_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_bubble_info);
  }

  return klass;
}

/* GfsInitGaussian: Header */

typedef struct _GfsInitGaussian         GfsInitGaussian;

struct _GfsInitGaussian {
  /*< private >*/
  GfsGenericInit parent;

  gdouble x, h, a;
  GfsVariable * c;
};

typedef struct _GfsInitGaussianClass    GfsInitGaussianClass;

struct _GfsInitGaussianClass {
  /*< private >*/
  GfsGenericInitClass parent_class;
};

#define GFS_INIT_GAUSSIAN(obj)            GTS_OBJECT_CAST (obj,\
					         GfsInitGaussian,\
					         gfs_init_gaussian_class ())
#define GFS_INIT_GAUSSIAN_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsInitGaussianClass,\
						 gfs_init_gaussian_class())
#define GFS_IS_INIT_GAUSSIAN(obj)         (gts_object_is_from_class (obj,\
						 gfs_init_gaussian_class ()))

GfsInitGaussianClass * gfs_init_gaussian_class  (void);

/* GfsInitGaussian: Object */

static void gfs_init_gaussian_read (GtsObject ** o, GtsFile * fp)
{
  GfsInitGaussian * b = GFS_INIT_GAUSSIAN (*o);
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x",        TRUE},
    {GTS_DOUBLE, "h",        TRUE},
    {GTS_DOUBLE, "a",        TRUE},
    {GTS_NONE}
  };

  if (GTS_OBJECT_CLASS (gfs_init_gaussian_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_gaussian_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  var[0].data = &b->x;
  var[1].data = &b->h;
  var[2].data = &b->a;
  gts_file_assign_variables (fp, var);
}

static void gfs_init_gaussian_write (GtsObject * o, FILE * fp)
{
  GfsInitGaussian * b = GFS_INIT_GAUSSIAN (o);
  if (GTS_OBJECT_CLASS (gfs_init_gaussian_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_init_gaussian_class ())->parent_class->write) 
      (o, fp);

  fprintf (fp, " { x = %g h = %g a = %g }", b->x, b->h, b->a);
}

static void init_gaussian (FttCell * cell, GfsInitGaussian * b)
{
  FttVector p;
  gdouble a = b->a*M_PI/180.;

  ftt_cell_pos (cell, &p);
  p.x = (p.x*cos (a) + p.y*sin (a) - b->x)/b->h;
  GFS_VARIABLE (cell, b->c->i) = exp (-p.x*p.x);
}

static gboolean gfs_init_gaussian_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_gaussian_class ())->parent_class)->event) (event, sim)) {
    GfsInitGaussian * gaussian = GFS_INIT_GAUSSIAN (event);

    gaussian->c = gfs_variable_from_name (GFS_DOMAIN (sim)->variables, "C");
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_gaussian,
			      gaussian);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_gaussian_class_init (GfsInitGaussianClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_gaussian_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_gaussian_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_gaussian_write;
}

static void gfs_init_gaussian_init (GfsInitGaussian * object)
{
  object->x = 0.0;
  object->h = 0.1;
  object->a = 0.0;
}

GfsInitGaussianClass * gfs_init_gaussian_class (void)
{
  static GfsInitGaussianClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_gaussian_info = {
      "GfsInitGaussian",
      sizeof (GfsInitGaussian),
      sizeof (GfsInitGaussianClass),
      (GtsObjectClassInitFunc) gfs_init_gaussian_class_init,
      (GtsObjectInitFunc) gfs_init_gaussian_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_gaussian_info);
  }

  return klass;
}

/* GfsInitErf: Header */

typedef struct _GfsInitErf         GfsInitErf;

struct _GfsInitErf {
  /*< private >*/
  GfsGenericInit parent;

  gdouble D, angle, scale;
  GfsVariable * c;
};

typedef struct _GfsInitErfClass    GfsInitErfClass;

struct _GfsInitErfClass {
  /*< private >*/
  GfsGenericInitClass parent_class;
};

#define GFS_INIT_ERF(obj)            GTS_OBJECT_CAST (obj,\
					         GfsInitErf,\
					         gfs_init_erf_class ())
#define GFS_INIT_ERF_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsInitErfClass,\
						 gfs_init_erf_class())
#define GFS_IS_INIT_ERF(obj)         (gts_object_is_from_class (obj,\
						 gfs_init_erf_class ()))

GfsInitErfClass * gfs_init_erf_class  (void);

/* GfsInitErf: Object */

static void gfs_init_erf_read (GtsObject ** o, GtsFile * fp)
{
  GfsInitErf * b = GFS_INIT_ERF (*o);

  if (GTS_OBJECT_CLASS (gfs_init_erf_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_erf_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (D)\n");
    return;
  }
  b->D = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (angle)\n");
    return;
  }
  b->angle = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (scale)\n");
    return;
  }
  b->scale = atof (fp->token->str);
  gts_file_next_token (fp);
}

static void gfs_init_erf_write (GtsObject * o, FILE * fp)
{
  GfsInitErf * b = GFS_INIT_ERF (o);
  if (GTS_OBJECT_CLASS (gfs_init_erf_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_init_erf_class ())->parent_class->write) 
      (o, fp);

  fprintf (fp, " %g %g %g", b->D, b->angle, b->scale);
}

static void init_erf (FttCell * cell, GfsInitErf * b)
{
  FttVector cm;
  gdouble psi;

  gfs_cell_cm (cell, &cm);
  psi = cos (b->angle)*cm.x + sin (b->angle)*cm.y;
  GFS_VARIABLE (cell, b->c->i) = 
    (erf ((0.75 - b->scale*psi)/sqrt (4.*b->D)) + 
     erf ((0.75 + b->scale*psi)/sqrt (4.*b->D)))/2.;
}

static gboolean gfs_init_erf_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_erf_class ())->parent_class)->event) (event, sim)) {
    GfsInitErf * erf = GFS_INIT_ERF (event);

    erf->c = gfs_variable_from_name (GFS_DOMAIN (sim)->variables, "C");
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_erf, erf);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_erf_class_init (GfsInitErfClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_erf_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_erf_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_erf_write;
}

static void gfs_init_erf_init (GfsInitErf * object)
{
  object->D = 0.01;
  object->angle = 0.;
  object->scale = 10.;
}

GfsInitErfClass * gfs_init_erf_class (void)
{
  static GfsInitErfClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_erf_info = {
      "GfsInitErf",
      sizeof (GfsInitErf),
      sizeof (GfsInitErfClass),
      (GtsObjectClassInitFunc) gfs_init_erf_class_init,
      (GtsObjectInitFunc) gfs_init_erf_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_erf_info);
  }

  return klass;
}

/* GfsOutputErfError: Header */

typedef struct _GfsOutputErfError         GfsOutputErfError;

struct _GfsOutputErfError {
  GfsOutput parent;

  GfsVariable * c;
  gdouble t, D, angle, scale, u;
};

#define GFS_OUTPUT_ERF_ERROR(obj)            GTS_OBJECT_CAST (obj,\
					           GfsOutputErfError,\
					gfs_output_erf_error_class ())
#define GFS_OUTPUT_ERF_ERROR_CLASS(klass)  GTS_OBJECT_CLASS_CAST (klass,\
						GfsOutputErfErrorClass,\
					   gfs_output_erf_error_class())
#define GFS_IS_OUTPUT_ERF_ERROR(obj)  (gts_object_is_from_class (obj,\
					 gfs_output_erf_error_class ()))


GfsOutputClass * gfs_output_erf_error_class  (void);

/* GfsOutputErfError: Object */

static void compute_erf_error (FttCell * cell, GfsOutputErfError * event)
{
  FttVector cm;
  gdouble psi;

  gfs_cell_cm (cell, &cm);
  psi = cos (event->angle)*cm.x + sin (event->angle)*cm.y -
    event->u*event->t;
  GFS_STATE (cell)->p = GFS_VARIABLE (cell, event->c->i)
    - (erf ((0.75 - event->scale*psi)/
	    sqrt (4.*event->D*(1. + event->t*event->scale*event->scale))) + 
       erf ((0.75 + event->scale*psi)/
	    sqrt (4.*event->D*(1. + event->t*event->scale*event->scale))))/2.;
}

static gboolean erf_error_event (GfsEvent * event, 
				 GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class ())->event) (event, sim)) {
    GfsNorm norm;
    
    GFS_OUTPUT_ERF_ERROR (event)->c = 
      gfs_variable_from_name (GFS_DOMAIN (sim)->variables, "C");
    GFS_OUTPUT_ERF_ERROR (event)->t = sim->time.t;
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) compute_erf_error, 
			      event);
    norm = gfs_domain_norm_variable (GFS_DOMAIN (sim),
				     gfs_p, FTT_TRAVERSE_LEAFS, -1);
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "%10e %10e %10e %10e %10e\n",
	     sim->time.t,
	     norm.first, norm.second, norm.infty, norm.bias);
    return TRUE;
  }
  return FALSE;
}

static void output_erf_error_write (GtsObject * object, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_erf_error_class ())->parent_class->write) (object, fp);
  fprintf (fp, " %g %g %g %g", 
	   GFS_OUTPUT_ERF_ERROR (object)->D, 
	   GFS_OUTPUT_ERF_ERROR (object)->angle,
	   GFS_OUTPUT_ERF_ERROR (object)->scale,
	   GFS_OUTPUT_ERF_ERROR (object)->u);
}

static void output_erf_error_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputErfError * output;

  (* GTS_OBJECT_CLASS (gfs_output_erf_error_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  output = GFS_OUTPUT_ERF_ERROR (*o);
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (D)");
    return;
  }
  output->D = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (angle)");
    return;
  }
  output->angle = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (scale)");
    return;
  }
  output->scale = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (u)");
    return;
  }
  output->u = atof (fp->token->str);
  gts_file_next_token (fp);
}

static void gfs_output_erf_error_class_init (GfsEventClass * klass)
{
  klass->event = erf_error_event;

  GTS_OBJECT_CLASS (klass)->write = output_erf_error_write;
  GTS_OBJECT_CLASS (klass)->read = output_erf_error_read;
}

static void gfs_output_erf_error_init (GfsOutputErfError * object)
{
  object->D = 0.01;
  object->angle = 0.;
  object->scale = 10.;
  object->u = 0.;
}

GfsOutputClass * gfs_output_erf_error_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_erf_error_info = {
      "GfsOutputErfError",
      sizeof (GfsOutputErfError),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_erf_error_class_init,
      (GtsObjectInitFunc) gfs_output_erf_error_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_erf_error_info);
  }

  return klass;
}

/* Initialize module */

const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_refine_box_class ();
  gfs_refine_sphere_class ();
  gfs_add_gaussian_vortex_class ();
  gfs_add_almgren_vortex_class ();
  gfs_init_rayleigh_class ();
  gfs_init_bubble_class ();
  gfs_init_gaussian_class ();
  gfs_init_erf_class ();
  gfs_output_erf_error_class ();
  return NULL;
}
