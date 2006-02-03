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

/* GfsVariableLevelSet: object */

static void variable_levelset_read (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain;

  (* GTS_OBJECT_CLASS (gfs_variable_levelset_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (c)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(GFS_VARIABLE_LEVELSET (*o)->v = 
	gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (level)");
    return;
  }
  GFS_VARIABLE_LEVELSET (*o)->level = atof (fp->token->str);
  gts_file_next_token (fp);
}

static void variable_levelset_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_levelset_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s %g", GFS_VARIABLE_LEVELSET (o)->v->name, GFS_VARIABLE_LEVELSET (o)->level);
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

static void levelset (FttCell * cell, GfsVariable * v)
{
  GfsVariableLevelSet * l = GFS_VARIABLE_LEVELSET (v);
  GtsPoint p;
  gdouble d2;
  
  ftt_cell_pos (cell, (FttVector *) &p.x);
  d2 = gfs_domain_cell_point_distance2 (v->domain, &p, vof_distance2, l->v, NULL);
  GFS_VARIABLE (cell, v->i) = GFS_VARIABLE (cell, l->v->i) > 0.5 ? sqrt (d2) : -sqrt (d2);
}

static void variable_levelset_event_half (GfsEvent * event, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsVariableLevelSet * v = GFS_VARIABLE_LEVELSET (event);

  gfs_domain_timer_start (domain, "levelset");

  gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
  			    (FttCellTraverseFunc) v->v->fine_coarse, v->v);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) levelset, event);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, GFS_VARIABLE1 (event));

  gfs_domain_timer_stop (domain, "levelset");
}

static gboolean variable_levelset_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_levelset_class ())->parent_class)->event)
      (event, sim)) {
    if (!GFS_VARIABLE_LEVELSET (event)->first_done) {
      variable_levelset_event_half (event, sim);
      GFS_VARIABLE_LEVELSET (event)->first_done = TRUE;
    }
    return TRUE;
  }
  return FALSE;
}

static void variable_levelset_class_init (GtsObjectClass * klass)
{
  klass->read = variable_levelset_read;
  klass->write = variable_levelset_write;
  GFS_EVENT_CLASS (klass)->event = variable_levelset_event;
  GFS_EVENT_CLASS (klass)->event_half = variable_levelset_event_half;
}

GfsVariableClass * gfs_variable_levelset_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_levelset_info = {
      "GfsVariableLevelSet",
      sizeof (GfsVariableLevelSet),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_levelset_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_levelset_info);
  }

  return klass;
}
