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

#include "tension.h"
#include "vof.h"

/* GfsSourceTension: Object */

static void gfs_source_tension_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceTension * s = GFS_SOURCE_TENSION (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  (* GTS_OBJECT_CLASS (gfs_source_tension_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable (C)");
    return;
  }
  if ((s->c = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  
  if ((s->t[0] = gfs_variable_from_name (domain->variables, "_gfs_source_tension_x")) == NULL) 
    s->t[0] = gfs_domain_add_variable (domain, "_gfs_source_tension_x");
  if ((s->t[1] = gfs_variable_from_name (domain->variables, "_gfs_source_tension_y")) == NULL) 
    s->t[1] = gfs_domain_add_variable (domain, "_gfs_source_tension_y");
#if (!FTT_2D)
  if ((s->t[2] = gfs_variable_from_name (domain->variables, "_gfs_source_tension_z")) == NULL) 
    s->t[2] = gfs_domain_add_variable (domain, "_gfs_source_tension_z");
#endif /* 3D */

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (sigma)");
    return;
  }
  s->sigma = atof (fp->token->str);

  gts_file_next_token (fp);
}

static void gfs_source_tension_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_tension_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s %g", GFS_SOURCE_TENSION (o)->c->name, GFS_SOURCE_TENSION (o)->sigma);
}

static void foreach_cell_normal (FttCell * cell, GfsVariable * v)
{
  FttVector n;
  gdouble nn = 0.;
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++) {
    (&n.x)[c] = gfs_youngs_gradient (cell, c, v);
    nn += (&n.x)[c]*(&n.x)[c];
  }
  nn = sqrt (nn + 1e-50);
  GFS_STATE (cell)->g[0] = n.x*n.x/nn;
  GFS_STATE (cell)->g[1] = n.y*n.y/nn;
  GFS_STATE (cell)->div = n.x*n.y/nn;
}

static void foreach_cell_tension (FttCell * cell, GfsSourceTension * s)
{
  gdouble h = ftt_cell_size (cell);
  gdouble sigh = s->sigma/(h*h*h);

  GFS_VARIABLE (cell, s->t[0]->i) = sigh*(gfs_youngs_gradient (cell, FTT_X, gfs_gy) -
					  gfs_youngs_gradient (cell, FTT_Y, gfs_div));
  GFS_VARIABLE (cell, s->t[1]->i) = sigh*(gfs_youngs_gradient (cell, FTT_Y, gfs_gx) -
					  gfs_youngs_gradient (cell, FTT_X, gfs_div));
}

static gboolean gfs_source_tension_event (GfsEvent * event, 
					  GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_source_tension_class ())->parent_class)->event) 
      (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim),
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) foreach_cell_normal, 
			      GFS_SOURCE_TENSION (event)->c);
    /* fixme: boundary conditions for normal */
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) foreach_cell_tension, event);
    return TRUE;
  }
  return FALSE;
}

static gdouble gfs_source_tension_value (GfsSourceGeneric * s, 
					 FttCell * cell,
					 GfsVariable * v)
{
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++)
    if (GFS_SOURCE_VECTOR (s)->v[c] == v)
      return GFS_VARIABLE (cell, GFS_SOURCE_TENSION (s)->t[c]->i);
  g_assert_not_reached ();
  return 0;
}

static void gfs_source_tension_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read =  gfs_source_tension_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_source_tension_write;
  GFS_EVENT_CLASS (klass)->event =  gfs_source_tension_event;
  klass->centered_value =           gfs_source_tension_value;
}

GfsSourceGenericClass * gfs_source_tension_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_tension_info = {
      "GfsSourceTension",
      sizeof (GfsSourceTension),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_tension_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = 
      gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_vector_class ()),
			    &gfs_source_tension_info);
  }

  return klass;
}

