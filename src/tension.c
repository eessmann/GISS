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
  FttComponent c;

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

  for (c = 0; c < FTT_DIMENSION; c++) {
    static gchar * name[3] = {"_Tx", "_Ty", "_Tz"};
    if ((s->t[c] = gfs_variable_from_name (domain->variables, name[c])) == NULL)
      s->t[c] = gfs_domain_add_variable (domain, name[c]);
  }

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

static void foreach_cell_normal (FttCell * cell, GfsSourceTension * s)
{
  FttVector n;
  gdouble nn = 0.;
  gdouble sigh = s->sigma/ftt_cell_size (cell);
  FttComponent c;

  gfs_youngs_normal (cell, s->c, &n);
  for (c = 0; c < FTT_DIMENSION; c++)
    nn += (&n.x)[c]*(&n.x)[c];
  nn = sqrt (nn + 1e-50);
  GFS_STATE (cell)->g[0] = sigh*n.x*n.x/nn;
  GFS_STATE (cell)->g[1] = sigh*n.y*n.y/nn;
  GFS_STATE (cell)->div  = sigh*n.x*n.y/nn;
}

static void foreach_cell_tension (FttCell * cell, GfsSourceTension * s)
{
  gdouble h = ftt_cell_size (cell);
  FttVector nx, ny, nxy;

  gfs_youngs_normal (cell, gfs_gx, &nx);
  gfs_youngs_normal (cell, gfs_gy, &ny);
  gfs_youngs_normal (cell, gfs_div, &nxy);

  GFS_VARIABLE (cell, s->t[0]->i) = (ny.x - nxy.y)/h;
  GFS_VARIABLE (cell, s->t[1]->i) = (nx.y - nxy.x)/h;
}

static void gfs_source_tension_event (GfsEvent * event, 
				      GfsSimulation * sim)
{
  gfs_domain_cell_traverse (GFS_DOMAIN (sim),
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) foreach_cell_normal, event);
  /* fixme: boundary conditions for normal */
  gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) foreach_cell_tension, event);
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
  GTS_OBJECT_CLASS (klass)->read =       gfs_source_tension_read;
  GTS_OBJECT_CLASS (klass)->write =      gfs_source_tension_write;
  GFS_EVENT_CLASS (klass)->event_half =  gfs_source_tension_event;
  klass->centered_value =                gfs_source_tension_value;
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

