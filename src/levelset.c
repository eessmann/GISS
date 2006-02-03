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

static void distance (FttCell * cell, GfsVariable * v)
{
  GfsVariableDistance * l = GFS_VARIABLE_DISTANCE (v);
  GtsPoint p;
  gdouble d2;
  
  ftt_cell_pos (cell, (FttVector *) &p.x);
  d2 = gfs_domain_cell_point_distance2 (v->domain, &p, vof_distance2, l->v, NULL);
  GFS_VARIABLE (cell, v->i) = GFS_VARIABLE (cell, l->v->i) > 0.5 ? sqrt (d2) : -sqrt (d2);
}

static void variable_distance_event_half (GfsEvent * event, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsVariableDistance * v = GFS_VARIABLE_DISTANCE (event);

  gfs_domain_timer_start (domain, "distance");

  gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
  			    (FttCellTraverseFunc) v->v->fine_coarse, v->v);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) distance, event);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, GFS_VARIABLE1 (event));

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
  GfsDomain * domain;

  (* GTS_OBJECT_CLASS (gfs_variable_curvature_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (d)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(GFS_VARIABLE_CURVATURE (*o)->d = 
	gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  if (!GFS_IS_VARIABLE_DISTANCE (GFS_VARIABLE_CURVATURE (*o)->d)) {
    gts_file_error (fp, "variable `%s' is not a GfsVariableDistance", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_DOUBLE, "sigma", TRUE},
      {GTS_DOUBLE, "theta", TRUE},
      {GTS_NONE}
    };

    var[0].data = &GFS_VARIABLE_CURVATURE (*o)->sigma;
    var[1].data = &GFS_VARIABLE_CURVATURE (*o)->theta;
    gts_file_assign_variables (fp, var);
  }
}

static void variable_curvature_write (GtsObject * o, FILE * fp)
{
  GfsVariableCurvature * v = GFS_VARIABLE_CURVATURE (o);

  (* GTS_OBJECT_CLASS (gfs_variable_curvature_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s", v->d->name);
  if (v->sigma != 1. || v->theta != 0.5)
    fprintf (fp, " { sigma = %g theta = %g }", v->sigma, v->theta);
}

static void normal (FttCell * cell, gpointer * data)
{
  GfsVariable ** nv = data[0];
  GfsVariable * d = GFS_VARIABLE_CURVATURE (data[1])->d;
  GtsVector n = { 0., 0., 0. };
  FttComponent c;

  gfs_youngs_normal (cell, d, (FttVector *) n);
  gts_vector_normalize (n);
  for (c = 0; c < FTT_DIMENSION; c++)
    GFS_VARIABLE (cell, nv[c]->i) = n[c];
}

static void curvature (FttCell * cell, gpointer * data)
{
  GfsVariable ** nv = data[0];
  gdouble kappa = 0.;
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++)
    kappa += gfs_center_gradient (cell, c, nv[c]->i);
  GFS_VARIABLE (cell, nv[FTT_DIMENSION]->i) = kappa/ftt_cell_size (cell);
}

static void interface_curvature (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[1];
  GfsVariableCurvature * k = GFS_VARIABLE_CURVATURE (v);
  gdouble f = GFS_VARIABLE (cell, GFS_VARIABLE_DISTANCE (k->d)->v->i);

  if (GFS_IS_FULL (f))
    GFS_VARIABLE (cell, v->i) = G_MAXDOUBLE;
  else {
    GfsVariable ** nv = data[0];
    FttComponent c;
    FttVector p;
    gdouble kappa;

    ftt_cell_pos (cell, &p);
    for (c = 0; c < FTT_DIMENSION; c++)
      (&p.x)[c] -= GFS_VARIABLE (cell, k->d->i)*GFS_VARIABLE (cell, nv[c]->i);
    kappa = k->sigma*gfs_interpolate (cell, p, nv[FTT_DIMENSION]);
    if (GFS_VARIABLE (cell, v->i) < G_MAXDOUBLE)
      GFS_VARIABLE (cell, v->i) = (k->a*kappa + (1. - k->a)*GFS_VARIABLE (cell, v->i));
    else
      GFS_VARIABLE (cell, v->i) = kappa;
  }
}

static void variable_curvature_event_half (GfsEvent * event, GfsSimulation * sim)
{
  GfsVariable * n[FTT_DIMENSION + 1];
  GfsDomain * domain = GFS_DOMAIN (sim);
  gpointer data[2];
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION + 1; c++) {
    n[c] = gfs_temporary_variable (domain);
    gfs_variable_set_vector (n[c], c);
  }
  data[0] = n;
  data[1] = event;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) normal, data);
  for (c = 0; c < FTT_DIMENSION; c++)
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, n[c]);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) curvature, data);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) interface_curvature, data);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, GFS_VARIABLE1 (event));
  for (c = 0; c < FTT_DIMENSION + 1; c++)
    gts_object_destroy (GTS_OBJECT (n[c]));
}

static gboolean variable_curvature_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_curvature_class ())->parent_class)->event)
      (event, sim)) {
    GFS_VARIABLE_CURVATURE (event)->a = 1.;
    variable_curvature_event_half (event, sim);
    GFS_VARIABLE_CURVATURE (event)->a = GFS_VARIABLE_CURVATURE (event)->theta;
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

static void variable_curvature_init (GfsVariableCurvature * v)
{
  v->sigma = 1.;
  v->theta = 0.5;
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
