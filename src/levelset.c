/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2006 National Institute of Water and Atmospheric Research
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
#include "levelset.h"
#include "vof.h"

/* GfsVariableDistance: object */

static void variable_distance_read (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain;

  (* GTS_OBJECT_CLASS (gfs_variable_distance_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (c)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(GFS_VARIABLE_DISTANCE (*o)->v = 
	gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
}

static void variable_distance_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_distance_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s", GFS_VARIABLE_DISTANCE (o)->v->name);
}

static gdouble vof_distance2 (FttCell * cell, GtsPoint * t, gpointer v)
{
  gdouble f = GFS_VARIABLE (cell, GFS_VARIABLE1 (v)->i);
  
  if (GFS_IS_FULL (f))
    return G_MAXDOUBLE;
  if (!FTT_CELL_IS_LEAF (cell))
    return ftt_cell_point_distance2_min (cell, t);
  else {
    gdouble h = ftt_cell_size (cell), d2;
    GSList * l = gfs_vof_facet (cell, v);
    GtsPoint * p1 = l->data, * p2 = l->next->data;
    GtsSegment s;
    FttVector p;

#if FTT_3D
  g_assert_not_implemented ();
#endif

    ftt_cell_pos (cell, &p);
    p1->x = p.x + h*p1->x; p1->y = p.y + h*p1->y;
    p2->x = p.x + h*p2->x; p2->y = p.y + h*p2->y;
    s.v1 = (GtsVertex *) p1; s.v2 = (GtsVertex *) p2;
    d2 = gts_point_segment_distance2 (t, &s);
    gts_object_destroy (GTS_OBJECT (p1));
    gts_object_destroy (GTS_OBJECT (p2));
    g_slist_free (l);
    return d2;
  }
}

static void distance (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsVariable * s2 = data[2];

  if (GFS_VARIABLE (cell, s2->i)) {
    GfsVariableDistance * l = GFS_VARIABLE_DISTANCE (v);
    GtsPoint p;
    gdouble d2;
    
    ftt_cell_pos (cell, (FttVector *) &p.x);
    d2 = gfs_domain_cell_point_distance2 (v->domain, &p, vof_distance2, l->v, NULL);
    GFS_VARIABLE (cell, v->i) = GFS_VARIABLE (cell, l->v->i) > 0.5 ? sqrt (d2) : -sqrt (d2);
  }
  else
    GFS_VARIABLE (cell, v->i) = G_MAXDOUBLE;
}

static void stencil_interpolate (FttCell * cell, gpointer * data)
{
  GfsVariableDistance * v = data[0];
  gdouble f = GFS_VARIABLE (cell, v->v->i);
  
  if (!GFS_IS_FULL (f))
    gfs_interpolate_stencil (cell, data[1]);
}

static void stencil_gradient (FttCell * cell, gpointer * data)
{
  GfsVariable * s1 = data[1];
  GfsVariable * s2 = data[2];
  FttComponent c;

  if (GFS_VARIABLE (cell, s1->i))
    for (c = 0; c < FTT_DIMENSION; c++)
      gfs_center_gradient_stencil (cell, c, s2->i);
}

static void variable_distance_event_half (GfsEvent * event, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsVariableDistance * v = GFS_VARIABLE_DISTANCE (event);
  gpointer data[3], tmp;

  gfs_domain_timer_start (domain, "distance");

  data[0] = v;
  data[1] = gfs_temporary_variable (domain);
  data[2] = gfs_temporary_variable (domain);
  gfs_domain_cell_traverse (domain,  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, data[1]);
  gfs_domain_cell_traverse (domain,  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) stencil_interpolate, data);
  gfs_domain_cell_traverse (domain,  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, data[2]);
  gfs_domain_cell_traverse (domain,  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) stencil_gradient, data);
  tmp = data[1]; data[1] = data[2]; data[2] = tmp;
  gfs_domain_cell_traverse (domain,  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) stencil_gradient, data);

  gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
  			    (FttCellTraverseFunc) v->v->fine_coarse, v->v);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) distance, data);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, GFS_VARIABLE1 (event));

  gts_object_destroy (data[1]);
  gts_object_destroy (data[2]);

  gfs_domain_timer_stop (domain, "distance");
}

static gboolean variable_distance_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_distance_class ())->parent_class)->event)
      (event, sim)) {
    if (!GFS_VARIABLE_DISTANCE (event)->first_done) {
      variable_distance_event_half (event, sim);
      GFS_VARIABLE_DISTANCE (event)->first_done = TRUE;
    }
    return TRUE;
  }
  return FALSE;
}

static void variable_distance_class_init (GtsObjectClass * klass)
{
  klass->read = variable_distance_read;
  klass->write = variable_distance_write;
  GFS_EVENT_CLASS (klass)->event = variable_distance_event;
  GFS_EVENT_CLASS (klass)->event_half = variable_distance_event_half;
}

GfsVariableClass * gfs_variable_distance_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_distance_info = {
      "GfsVariableDistance",
      sizeof (GfsVariableDistance),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_distance_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_distance_info);
  }

  return klass;
}

/* GfsVariableCurvature: object */

static void variable_curvature_read (GtsObject ** o, GtsFile * fp)
{
  GfsVariableCurvature * v = GFS_VARIABLE_CURVATURE (*o);
  GfsDomain * domain;

  (* GTS_OBJECT_CLASS (gfs_variable_curvature_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (f)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(v->f = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
}

static void variable_curvature_write (GtsObject * o, FILE * fp)
{
  GfsVariableCurvature * v = GFS_VARIABLE_CURVATURE (o);

  (* GTS_OBJECT_CLASS (gfs_variable_curvature_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s", v->f->name);
}

static void curvature (FttCell * cell, GfsVariable * v)
{
  GfsVariable * t = GFS_VARIABLE_CURVATURE (v)->f;
  gdouble f = GFS_VARIABLE (cell, t->i);

  if (GFS_IS_FULL (f))
    GFS_VARIABLE (cell, v->i) = G_MAXDOUBLE;
  else
    GFS_VARIABLE (cell, v->i) = gfs_height_curvature (cell, t);
}

#define THETA 0.5

static void filter (FttCell * cell, gpointer * data)
{
  GfsVariable * tmp = data[0], * v = data[1];
  GfsVariable * t = GFS_VARIABLE_CURVATURE (v)->f;
  gdouble f = GFS_VARIABLE (cell, t->i);

  if (GFS_IS_FULL (f))
    GFS_VARIABLE (cell, tmp->i) = G_MAXDOUBLE;
  else {
    gdouble h = ftt_cell_size (cell);
    guint level = ftt_cell_level (cell);
    FttVector p;
    gint x, y;
    gdouble w = 0., st = 0.;
    
    ftt_cell_pos (cell, &p);
    for (x = -1; x <= 1; x++)
      for (y = -1; y <= 1; y++) 
	if (x != 0 || y != 0) {
	  FttVector o;
	  o.x = p.x + h*x; o.y = p.y + h*y; o.z = 0.;
	  FttCell * neighbor = gfs_domain_locate (v->domain, o, level);
	  g_assert (neighbor);
	  g_assert (ftt_cell_level (neighbor) == level);
	  f = GFS_VARIABLE (neighbor, t->i);
	  if (!GFS_IS_FULL (f)) {
	    w += f*(1. - f);
	    st += f*(1. - f)*GFS_VARIABLE (neighbor, v->i);
	  }
	}
    g_assert (w > 0.);
    GFS_VARIABLE (cell, tmp->i) = (1. - THETA)*GFS_VARIABLE (cell, v->i) + THETA*st/w;
  }
}

static void filter_curvature (GfsDomain * domain, GfsVariable * k)
{
  gpointer data[2];

  data[0] = gfs_temporary_variable (domain);
  data[1] = k;
  guint i = 1;
  while (i--) {
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) filter, data);
    gfs_variables_swap (data[1], data[0]);
  }
  gts_object_destroy (data[0]);  
}

static void variable_curvature_event_half (GfsEvent * event, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) curvature, event);
  // filter_curvature (domain, GFS_VARIABLE1 (event));
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, GFS_VARIABLE1 (event));
}

static gboolean variable_curvature_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_curvature_class ())->parent_class)->event)
      (event, sim)) {
    if (!GFS_VARIABLE_CURVATURE (event)->first_done) {
      variable_curvature_event_half (event, sim);
      GFS_VARIABLE_CURVATURE (event)->first_done = TRUE;
    }
    return TRUE;
  }
  return FALSE;
}

static void variable_curvature_class_init (GtsObjectClass * klass)
{
  klass->read = variable_curvature_read;
  klass->write = variable_curvature_write;
  GFS_EVENT_CLASS (klass)->event = variable_curvature_event;
  GFS_EVENT_CLASS (klass)->event_half = variable_curvature_event_half;
}

static void curvature_coarse_fine (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  guint n;

  ftt_cell_children (parent, &child);
  for (n = 0; n < FTT_CELLS; n++) {
    GfsVariableCurvature * k = GFS_VARIABLE_CURVATURE (v);
    gdouble f = GFS_VARIABLE (child.c[n], k->f->i);

    GFS_VARIABLE (child.c[n], v->i) = GFS_VARIABLE (parent, v->i);
    if (!GFS_IS_FULL (f))
      g_assert (fabs (GFS_VARIABLE (child.c[n], v->i)) < G_MAXDOUBLE);
  }
}

static void curvature_fine_coarse (FttCell * parent, GfsVariable * v)
{
  GfsVariableCurvature * k = GFS_VARIABLE_CURVATURE (v);
  GfsVariable * t = k->f;
  FttCellChildren child;
  gdouble val = 0., sa = 0.;
  guint i;

  ftt_cell_children (parent, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i] && !GFS_IS_FULL (GFS_VARIABLE (child.c[i], t->i))) {
      gdouble a = GFS_IS_MIXED (child.c[i]) ? GFS_STATE (child.c[i])->solid->a : 1.;

      val += GFS_VARIABLE (child.c[i], v->i)*a;
      sa += a;
    }
  if (sa > 0.)
    GFS_VARIABLE (parent, v->i) = val/sa;
  else {
    GFS_VARIABLE (parent, v->i) = G_MAXDOUBLE;
    g_assert (GFS_IS_FULL (GFS_VARIABLE (parent, t->i)));
  }
}

static void variable_curvature_init (GfsVariableCurvature * v)
{
  GFS_VARIABLE1 (v)->coarse_fine = curvature_coarse_fine;
  GFS_VARIABLE1 (v)->fine_coarse = curvature_fine_coarse;
}

GfsVariableClass * gfs_variable_curvature_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_curvature_info = {
      "GfsVariableCurvature",
      sizeof (GfsVariableCurvature),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_curvature_class_init,
      (GtsObjectInitFunc) variable_curvature_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_curvature_info);
  }

  return klass;
}
