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

static void surface_add_box (GtsSurface * s,
			     gdouble x1, gdouble y1, gdouble z1,
			     gdouble x2, gdouble y2, gdouble z2)
{
  GtsVertex * v0 = gts_vertex_new (s->vertex_class, x1, y1, z1);
  GtsVertex * v1 = gts_vertex_new (s->vertex_class, x1, y1, z2);
  GtsVertex * v2 = gts_vertex_new (s->vertex_class, x1, y2, z2);
  GtsVertex * v3 = gts_vertex_new (s->vertex_class, x1, y2, z1);
  GtsVertex * v4 = gts_vertex_new (s->vertex_class, x2, y1, z1);
  GtsVertex * v5 = gts_vertex_new (s->vertex_class, x2, y1, z2);
  GtsVertex * v6 = gts_vertex_new (s->vertex_class, x2, y2, z2);
  GtsVertex * v7 = gts_vertex_new (s->vertex_class, x2, y2, z1);

  GtsEdge * e1 = gts_edge_new (s->edge_class, v0, v1);
  GtsEdge * e2 = gts_edge_new (s->edge_class, v1, v2);
  GtsEdge * e3 = gts_edge_new (s->edge_class, v2, v3);
  GtsEdge * e4 = gts_edge_new (s->edge_class, v3, v0);
  GtsEdge * e5 = gts_edge_new (s->edge_class, v0, v2);

  GtsEdge * e6 = gts_edge_new (s->edge_class, v4, v5);
  GtsEdge * e7 = gts_edge_new (s->edge_class, v5, v6);
  GtsEdge * e8 = gts_edge_new (s->edge_class, v6, v7);
  GtsEdge * e9 = gts_edge_new (s->edge_class, v7, v4);
  GtsEdge * e10 = gts_edge_new (s->edge_class, v4, v6);
  
  GtsEdge * e11 = gts_edge_new (s->edge_class, v3, v7);
  GtsEdge * e12 = gts_edge_new (s->edge_class, v2, v6);
  GtsEdge * e13 = gts_edge_new (s->edge_class, v1, v5);
  GtsEdge * e14 = gts_edge_new (s->edge_class, v0, v4);

  GtsEdge * e15 = gts_edge_new (s->edge_class, v1, v6);
  GtsEdge * e16 = gts_edge_new (s->edge_class, v2, v7);
  GtsEdge * e17 = gts_edge_new (s->edge_class, v3, v4);
  GtsEdge * e18 = gts_edge_new (s->edge_class, v0, v5);

  GtsFaceClass * klass = gts_face_class ();

  gts_surface_add_face (s, gts_face_new (klass, e1, e2, e5));
  gts_surface_add_face (s, gts_face_new (klass, e5, e3, e4));
  gts_surface_add_face (s, gts_face_new (klass, e6, e10, e7));
  gts_surface_add_face (s, gts_face_new (klass, e10, e9, e8));
  gts_surface_add_face (s, gts_face_new (klass, e2, e15, e12));
  gts_surface_add_face (s, gts_face_new (klass, e15, e13, e7));
  gts_surface_add_face (s, gts_face_new (klass, e3, e16, e11));
  gts_surface_add_face (s, gts_face_new (klass, e16, e12, e8));
  gts_surface_add_face (s, gts_face_new (klass, e17, e14, e4));
  gts_surface_add_face (s, gts_face_new (klass, e17, e11, e9));
  gts_surface_add_face (s, gts_face_new (klass, e18, e13, e1));
  gts_surface_add_face (s, gts_face_new (klass, e18, e14, e6));
}

static GtsBBox * bbox_cell (GtsBBoxClass * klass, FttCell * cell)
{
  FttVector p;
  gdouble size;

  ftt_cell_pos (cell, &p);
  size = ftt_cell_size (cell)/2.;
  return gts_bbox_new (klass, cell, 
		       p.x - size, p.y - size, p.z - size,
		       p.x + size, p.y + size, p.z + size);
}

static void Knds (GtsTriangle * t, FttVector * te)
{
  GtsVertex * v1, * v2, * v3;
  FttComponent c;
  gdouble a = gts_triangle_area (t);

  gts_triangle_vertices (t, &v1, &v2, &v3);
  for (c = 0; c < FTT_DIMENSION; c++)
    (&te->x)[c] += a*(GTS_VERTEX_NORMAL (v1)->n[c] +
		      GTS_VERTEX_NORMAL (v2)->n[c] +
		      GTS_VERTEX_NORMAL (v3)->n[c]);
}

static void cell_tension (FttCell * cell, 
			  GtsSurface * s,
			  GNode * stree,
			  FttVector * t)
{
  GtsBBox * bbox;

  t->x = t->y = t->z = 0.;
  bbox = bbox_cell (gts_bbox_class (), cell);
  if (gts_bb_tree_is_overlapping (stree, bbox)) {
    GtsSurface * cs;
    GNode * ctree;
    GtsSurfaceInter * si;
    
    if (GFS_IS_MIXED (cell))
      g_assert_not_implemented ();

    cs = gts_surface_new (gts_surface_class (),
			  gts_face_class (),
			  gts_edge_class (),
			  gts_vertex_class ());
    surface_add_box (cs, 
		     bbox->x1, bbox->y1, bbox->z1,
		     bbox->x2, bbox->y2, bbox->z2);    
    ctree = gts_bb_tree_surface (cs);
    si = gts_surface_inter_new (gts_surface_inter_class (),
				s, cs, stree, ctree, FALSE, FALSE);
    if (si->edges != NULL) {
      GtsSurface * inter = gts_surface_new (gts_surface_class (),
					    gts_face_class (),
					    gts_edge_class (),
					    gts_vertex_class ());
      FttComponent c;
      gdouble h = ftt_cell_size (cell);
      gdouble vol = 6.*h*h*h;

      gts_surface_inter_boolean (si, inter, GTS_1_IN_2);
      gts_surface_foreach_face (inter, (GtsFunc) Knds, t);
      for (c = 0; c < FTT_DIMENSION; c++)
	(&t->x)[c] /= vol;
#if 0
      {
	GtsVector p;

	gts_surface_center_of_area (inter, p);
	printf ("VECT 1 2 0 2 0 %g %g %g %g %g %g\n",
		p[0], p[1], p[2],
		p[0] - 3.*t->x,
		p[1] - 3.*t->y,
		p[2] - 3.*t->z);
      }
#endif
      gts_object_destroy (GTS_OBJECT (inter));
    }
    gts_object_destroy (GTS_OBJECT (si));
    gts_bb_tree_destroy (ctree, TRUE);
    gts_object_destroy (GTS_OBJECT (cs));
  }
  gts_object_destroy (GTS_OBJECT (bbox));
}

static void foreach_cell_tension (FttCell * cell, gpointer * data)
{
  FttVector t;
  GfsSimulation * sim = data[0];
  GfsSourceTension * s = data[1];
  FttComponent c;

  cell_tension (cell, sim->interface, sim->itree, &t);
  for (c = 0; c < FTT_DIMENSION; c++)
    GFS_VARIABLE (cell, s->t[c]->i) = s->sigma*(&t.x)[c];
}

static void vertex_normal (GtsVertex * v, GtsSurface * s)
{
  gts_vertex_mean_curvature_normal (v, s, GTS_VERTEX_NORMAL (v)->n);
}

/* GfsSourceTension: Object */

static void gfs_source_tension_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceTension * s = GFS_SOURCE_TENSION (*o);

  if (GTS_OBJECT_CLASS (gfs_source_tension_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_tension_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  s->t[0] = gfs_domain_add_variable (GFS_DOMAIN (gfs_object_simulation (*o)),
				     "_gfs_source_tension_x");
  if (s->t[0] == NULL) {
    gts_file_error (fp, "only one GfsSourceTension is allowed");
    return;
  }
  s->t[1] = gfs_domain_add_variable (GFS_DOMAIN (gfs_object_simulation (*o)),
				     "_gfs_source_tension_y");
  g_assert (s->t[1]);
#if (!FTT_2D)
  s->t[2] = gfs_domain_add_variable (GFS_DOMAIN (gfs_object_simulation (*o)),
				     "_gfs_source_tension_z");
  g_assert (s->t[2]);
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
  if (GTS_OBJECT_CLASS (gfs_source_tension_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_source_tension_class ())->parent_class->write) 
      (o, fp);
  fprintf (fp, " %g", GFS_SOURCE_TENSION (o)->sigma);
}

static gboolean gfs_source_tension_event (GfsEvent * event, 
					  GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_source_tension_class ())->parent_class)->event) (event, sim)) {
    if (sim->interface) {
      gpointer data[2];

      gts_surface_foreach_vertex (sim->interface, 
				  (GtsFunc) vertex_normal, sim->interface);
      data[0] = sim;
      data[1] = event;
      gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) foreach_cell_tension, 
				data);
    }
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

