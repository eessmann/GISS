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


#include "simulation.h"
#include "surface.h"

/* GfsSurface: Object */

static void check_solid_surface (GtsSurface * s, 
				 const gchar * fname,
				 GtsFile * fp)
{
  GString * name = g_string_new ("surface");

  if (fname) {
    g_string_append (name, " `");
    g_string_append (name, fname);
    g_string_append_c (name, '\'');
  }

  if (!gts_surface_is_orientable (s))
    gts_file_error (fp, "%s is not orientable", name->str);
  g_string_free (name, TRUE);
}

static void surface_read (GtsObject ** o, GtsFile * fp)
{
  GfsSurface * surface = GFS_SURFACE (*o);

  if (fp->type == '(') { /* implicit surface */
    gts_file_next_token (fp);
    if (surface->f)
      gts_object_destroy (GTS_OBJECT (surface->f));
    surface->f = gfs_function_new (gfs_function_spatial_class (), 0.);
    gfs_function_read (surface->f, gfs_object_simulation (*o), fp);
    if (fp->type == GTS_ERROR)
      return;
    if (fp->type != ')') {
      gts_file_error (fp, "expecting a closing bracket");
      return;
    }
  }
  else if (fp->type == '{') { /* embedded surface */
    fp->scope_max++;
    gts_file_next_token (fp);
    if (surface->s)
      gts_object_destroy (GTS_OBJECT (surface->s));
    surface->s = gts_surface_new (gts_surface_class (), 
				  gts_face_class (), 
				  gts_edge_class (), 
				  gts_vertex_class ());
    if (gts_surface_read (surface->s, fp))
      return;
    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    check_solid_surface (surface->s, NULL, fp);
    if (fp->type == GTS_ERROR)
      return;
    fp->scope_max--;
  }
  else { /* surface file name */
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a string (filename)");
      return;
    }
    FILE * fptr = fopen (fp->token->str, "rt");
    if (fptr == NULL) {
      gts_file_error (fp, "cannot open file `%s'", fp->token->str);
      return;
    }
    GtsFile * fp1 = gts_file_new (fptr);
    surface->s = gts_surface_new (gts_surface_class (), 
				  gts_face_class (), 
				  gts_edge_class (), 
				  gts_vertex_class ());
    if (gts_surface_read (surface->s, fp1)) {
      gts_file_error (fp, 
		      "file `%s' is not a valid GTS file\n"
		      "%s:%d:%d: %s",
		      fp->token->str, fp->token->str,
		      fp1->line, fp1->pos, fp1->error);
      gts_file_destroy (fp1);
      fclose (fptr);
      return;
    }
    gts_file_destroy (fp1);
    fclose (fptr);
    check_solid_surface (surface->s, fp->token->str, fp);
    if (fp->type == GTS_ERROR)
      return;
  }
  gts_file_next_token (fp);

  if (fp->type == '{') {
    gdouble scale = 1.;
    GtsFileVariable var[] = {
      {GTS_DOUBLE, "rx", TRUE},
      {GTS_DOUBLE, "ry", TRUE},
      {GTS_DOUBLE, "rz", TRUE},
      {GTS_DOUBLE, "sx", TRUE},
      {GTS_DOUBLE, "sy", TRUE},
      {GTS_DOUBLE, "sz", TRUE},
      {GTS_DOUBLE, "tx", TRUE},
      {GTS_DOUBLE, "ty", TRUE},
      {GTS_DOUBLE, "tz", TRUE},
      {GTS_DOUBLE, "scale", TRUE},
      {GTS_INT,    "flip", TRUE},
      {GTS_INT,    "twod", TRUE},
      {GTS_NONE}
    };
    GtsFileVariable * v = var;

    (v++)->data = &surface->rotate[0];
    (v++)->data = &surface->rotate[1];
    (v++)->data = &surface->rotate[2];

    (v++)->data = &surface->scale[0];
    (v++)->data = &surface->scale[1];
    (v++)->data = &surface->scale[2];

    (v++)->data = &surface->translate[0];
    (v++)->data = &surface->translate[1];
    (v++)->data = &surface->translate[2];

    (v++)->data = &scale;

    (v++)->data = &surface->flip;
    (v++)->data = &surface->twod;

    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;

    if (var[9].set) {
      surface->scale[0] *= scale;
      surface->scale[1] *= scale;
      surface->scale[2] *= scale;
    }

    GtsMatrix * m = gts_matrix_translate (NULL, surface->translate);
    gint c;
    for (c = 2; c >= 0; c--)
      if (surface->rotate[c] != 0.) {
	GtsVector r = {0.,0.,0.};
	r[c] = 1.;
	GtsMatrix * mr = gts_matrix_rotate (NULL, r, surface->rotate[c]*M_PI/180.);
	GtsMatrix * m1 = gts_matrix_product (m, mr);
	gts_matrix_destroy (m);
	gts_matrix_destroy (mr);
	m = m1;
      }
    GtsMatrix * ms = gts_matrix_scale (NULL, surface->scale);
    if (surface->m)
      gts_matrix_destroy (surface->m);
    surface->m = gts_matrix_product (m, ms);
    gts_matrix_destroy (m);
    gts_matrix_destroy (ms);

    if (surface->s) {
      gts_surface_foreach_vertex (surface->s, (GtsFunc) gts_point_transform, surface->m);
      gts_matrix_destroy (surface->m);
      surface->m = NULL;
      if (surface->flip)
	gts_surface_foreach_face (surface->s, (GtsFunc) gts_triangle_revert, NULL);
    }
    else {
      GtsMatrix * i = gts_matrix_inverse (surface->m);
      gts_matrix_destroy (surface->m);
      surface->m = i;
    }
  }
}

static void surface_write (GtsObject * o, FILE * fp)
{
  GfsSurface * surface = GFS_SURFACE (o);
  if (surface->s) {
    fputs (" { ", fp);
    GtsSurface * s = surface->s;
    if (GFS_DOMAIN (gfs_object_simulation (o))->binary) {
      gboolean binary = GTS_POINT_CLASS (s->vertex_class)->binary;
      GTS_POINT_CLASS (s->vertex_class)->binary = TRUE;
      gts_surface_write (s, fp);
      GTS_POINT_CLASS (s->vertex_class)->binary = binary;
    }
    else
      gts_surface_write (s, fp);
    fputc ('}', fp);
  }
  else if (surface->f) {
    fputs (" (", fp);
    gfs_function_write (surface->f, fp);
    fputs (" )", fp);
  }
  if (surface->m) {
    fputs (" {\n", fp);
    if (gts_vector_norm (surface->translate) > 0.)
      fprintf (fp, "  tx = %g ty = %g tz = %g\n",
	       surface->translate[0], surface->translate[1], surface->translate[2]);
    if (surface->scale[0] != 1. || surface->scale[1] != 1. || surface->scale[2] != 1.)
      fprintf (fp, "  sx = %g sy = %g sz = %g\n",
	       surface->scale[0], surface->scale[1], surface->scale[2]);
    if (surface->rotate[0] != 0.)
      fprintf (fp, "  rx = %g\n", surface->rotate[0]);
    if (surface->rotate[1] != 0.)
      fprintf (fp, "  ry = %g\n", surface->rotate[1]);
    if (surface->rotate[2] != 0.)
      fprintf (fp, "  rz = %g\n", surface->rotate[2]);
    if (surface->flip)
      fputs ("  flip = 1\n", fp);
    if (surface->twod)
      fputs ("  twod = 1\n", fp);
    fputc ('}', fp);
  }
}

static void surface_destroy (GtsObject * object)
{
  GfsSurface * s = GFS_SURFACE (object);
  if (s->s)
    gts_object_destroy (GTS_OBJECT (s->s));
  if (s->f)
    gts_object_destroy (GTS_OBJECT (s->f));
  if (s->m)
    gts_matrix_destroy (s->m);

  (* GTS_OBJECT_CLASS (gfs_surface_class ())->parent_class->destroy) (object);
}

static void gfs_surface_class_init (GtsObjectClass * klass)
{
  klass->read = surface_read;
  klass->write = surface_write;
  klass->destroy = surface_destroy;
}

static void gfs_surface_init (GfsSurface * s)
{
  s->scale[0] = 1.; s->scale[1] = 1.; s->scale[2] = 1.;
  s->flip = FALSE;
}

GtsObjectClass * gfs_surface_class (void)
{
  static GtsObjectClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_surface_info = {
      "GfsSurface",
      sizeof (GfsSurface),
      sizeof (GtsObjectClass),
      (GtsObjectClassInitFunc) gfs_surface_class_init,
      (GtsObjectInitFunc) gfs_surface_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (gts_object_class (), &gfs_surface_info);
  }

  return klass;
}

/**
 * gfs_surface_read:
 * @s: a #GfsSurface.
 * @sim: a #GfsSimulation.
 * @fp: a #GtsFile.
 * 
 * Calls the read() method of @s.
 */
void gfs_surface_read (GfsSurface * s, gpointer sim, GtsFile * fp)
{
  GtsObject * o = (GtsObject *) s;

  g_return_if_fail (s != NULL);
  g_return_if_fail (fp != NULL);

  o->reserved = sim;
  (* GTS_OBJECT (s)->klass->read) (&o, fp);
}

/**
 * gfs_surface_write:
 * @s: a #GfsSurface.
 * @sim: a #GfsSimulation.
 * @fp: a file pointer.
 * 
 * Calls the write() method of @s.
 */
void gfs_surface_write (GfsSurface * s, gpointer sim, FILE * fp)
{
  g_return_if_fail (s != NULL);
  g_return_if_fail (fp != NULL);

  GTS_OBJECT (s)->reserved = sim;
  (* GTS_OBJECT (s)->klass->write) (GTS_OBJECT (s), fp);
}

/**
 * gfs_surface_implicit_value:
 * @s: an (implicit) #GfsSurface.
 * @p: a #GtsPoint.
 *
 * Returns: the value of the implicit surface a location @p.
 */
gdouble gfs_surface_implicit_value (GfsSurface * s, GtsPoint p)
{
  g_return_val_if_fail (s != NULL, 0.);
  g_return_val_if_fail (s->f != NULL, 0.);

  if (s->m)
    gts_point_transform (&p, s->m);
  return (s->flip ? -1. : 1.)*gfs_function_spatial_value (s->f, (FttVector *)&p.x);
}

static gdouble segment_triangle_intersection (GtsPoint * E, GtsPoint * D,
					      GtsTriangle * t,
					      gboolean * inside)
{
  GtsVertex * vA, * vB, * vC;
  GtsPoint * A, * B, * C;
  gint ABCE, ABCD, ADCE, ABDE, BCDE;
  GtsEdge * AB, * BC, * CA;
  gdouble a, b;
  gboolean reversed = FALSE;

  gts_triangle_vertices_edges (t, NULL, &vA, &vB, &vC, &AB, &BC, &CA);
  A = GTS_POINT (vA);
  B = GTS_POINT (vB);
  C = GTS_POINT (vC);
  ABCE = gts_point_orientation_3d_sos (A, B, C, E);
  ABCD = gts_point_orientation_3d_sos (A, B, C, D);
  if (ABCE < 0 || ABCD > 0) {
    GtsPoint * tmpp;
    gint tmp;

    tmpp = E; E = D; D = tmpp;
    tmp = ABCE; ABCE = ABCD; ABCD = tmp;
    reversed = TRUE;
  }
  if (ABCE < 0 || ABCD > 0)
    return -1.;
  ADCE = gts_point_orientation_3d_sos (A, D, C, E);
  if (ADCE < 0)
    return -1.;
  ABDE = gts_point_orientation_3d_sos (A, B, D, E);
  if (ABDE < 0)
    return -1.;
  BCDE = gts_point_orientation_3d_sos (B, C, D, E);
  if (BCDE < 0)
    return -1.;
  *inside = reversed ? (ABCD < 0) : (ABCE < 0);
  a = gts_point_orientation_3d (A, B, C, E);
  b = gts_point_orientation_3d (A, B, C, D);
  if (a != b)
    return reversed ? 1. - a/(a - b) : a/(a - b);
  /* D and E are contained within ABC */
  g_assert (a == 0.);
  return 0.5;
}

static void triangle_face_intersection (GtsTriangle * t, GfsSegment * I)
{
  gboolean ins;
  gdouble x = segment_triangle_intersection (I->E, I->D, t, &ins);
  
  if (x != -1.) {
    I->x += x; I->n++;
    I->inside += ins ? 1 : -1;
  }
}

static GtsPoint segment_intersection (GfsSegment * I)
{
  GtsPoint p;
  p.x = I->E->x + I->x*(I->D->x - I->E->x);
  p.y = I->E->y + I->x*(I->D->y - I->E->y);
  p.z = I->E->z + I->x*(I->D->z - I->E->z);
  return p;
}

/**
 * gfs_surface_segment_intersection:
 * @s: a #GfsSurface.
 * @I: a GfsSegment.
 *
 * Fills @I with the intersection of @s and @I.
 *
 * Returns: the number of times @s intersects @I.
 */
guint gfs_surface_segment_intersection (GfsSurface * s,
					GfsSegment * I)
{
  g_return_val_if_fail (s != NULL, 0);
  g_return_val_if_fail (I != NULL, 0);

  I->n = 0;
  I->x = 0.;
  I->inside = 0;

  if (s->s)
    gts_surface_foreach_face (s->s, (GtsFunc) triangle_face_intersection, I);
  else {
    gdouble vE = gfs_surface_implicit_value (s, *I->E);
    gdouble vD = gfs_surface_implicit_value (s, *I->D);
    
    if ((vE > 0. && vD <= 0.) || (vE <= 0. && vD > 0.)) {
      I->n = 1;
      I->inside = vE > 0. ? -1 : 1;

      /* secant-bisection root-finding */
      gdouble t, t1, t2, v1, v2;
      if (vE > vD) {
	v1 = vD; t1 = 1.;
	v2 = vE; t2 = 0.;
      }
      else {
	v1 = vE; t1 = 0.;
	v2 = vD; t2 = 1.;
      }
      I->x = (v1*t2 - v2*t1)/(v1 - v2);
      guint n = 0;
      do {
	t = I->x;
	gdouble v = gfs_surface_implicit_value (s, segment_intersection (I));
	if (v < 0.) {
	  v1 = v; t1 = I->x;
	}
	else {
	  v2 = v; t2 = I->x;
	}
	if (v2 > v1)
	  I->x = (v1*t2 - v2*t1)/(v1 - v2);
	n++;
      }
      while (fabs (t - I->x) > 1e-3 && n < 10);
      if (fabs (t - I->x) > 1e-3)
	g_warning ("gfs_surface_segment_intersection(): convergence could not be reached\n"
		   "after %d iterations", n);
    }
  }
  return I->n;
}

static void surface_normal (GtsTriangle * t, GtsVector n)
{
  GtsVector m;
  gts_triangle_normal (t, &m[0], &m[1], &m[2]);
  n[0] -= m[0];
  n[1] -= m[1];
  n[2] -= m[2];
}

/**
 * gfs_surface_segment_normal:
 * @s: a #GfsSurface.
 * @I: a GfsSegment.
 * @n: a #GtsVector.
 *
 * Fills @n with the normal to @s at the intersection of @s and @I.
 */
void gfs_surface_segment_normal (GfsSurface * s,
				 GfsSegment * I,
				 GtsVector n)
{
  g_return_if_fail (s != NULL);
  g_return_if_fail (I != NULL);
  g_return_if_fail (I->n > 0);
  g_return_if_fail (n != NULL);
  
  if (s->s) {
    n[0] = n[1] = n[2] = 0.;
    gts_surface_foreach_face (s->s, (GtsFunc) surface_normal, n);
  }
  else {
    FttComponent c;
    GtsPoint p = segment_intersection (I);
    for (c = 0; c < FTT_DIMENSION; c++) {
      GtsPoint p1 = p;
      (&p1.x)[c] -= 1e-4;
      gdouble v1 = gfs_surface_implicit_value (s, p1);
      (&p1.x)[c] += 2e-4;
      gdouble v2 = gfs_surface_implicit_value (s, p1);
      n[c] = v2 - v1;
    }
  }
}

static void face_overlaps_box (GtsTriangle * t, gpointer * data)
{
  GtsBBox * bb = data[0];
  GtsSurface ** s1 = data[1];

  if (gts_bbox_overlaps_triangle (bb, t)) {
    if (*s1 == NULL)
      *s1 = gts_surface_new (gts_surface_class (),
			     gts_face_class (),
			     gts_edge_class (),
			     gts_vertex_class ());
    gts_surface_add_face (*s1, GTS_FACE (t));
  }
}

/**
 * gfs_cell_is_cut:
 * @cell: a #FttCell.
 * @s: a #GfsSurface.
 * @flatten: if set to %TRUE, @cell is flattened in the z direction.
 *
 * Returns: a (possibly new) #GfsSurface containing a subset of @s which may
 * intersect @cell or %NULL if @s does not intersect @cell.
 */
GfsSurface * gfs_cell_is_cut (FttCell * cell, GfsSurface * s, gboolean flatten)
{
  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (s != NULL, NULL);

  if (s->s) {
    GtsSurface * s1 = NULL;
    gpointer data[2];
    GtsBBox bb;
    ftt_cell_bbox (cell, &bb);
    if (flatten)
      bb.z1 = bb.z2 = 0.;
    data[0] = &bb;
    data[1] = &s1;
    gts_surface_foreach_face (s->s, (GtsFunc) face_overlaps_box, data);
    if (s1 == NULL)
      return NULL;
    GfsSurface * s2 = GFS_SURFACE (gts_object_new (gfs_surface_class ()));
    s2->s = s1;
    return s2;
  }
  else if (s->f)
    return s;
  g_assert_not_reached ();
  return NULL;
}

static void cell_traverse_cut (FttCell * cell,
			       GfsSurface * s,
			       FttTraverseType order,
			       FttTraverseFlags flags,
			       FttCellTraverseCutFunc func,
			       gpointer data,
			       gboolean flatten)
{
  GfsSurface * s1 = gfs_cell_is_cut (cell, s, flatten && FTT_CELL_IS_LEAF (cell));

  if (s1 == NULL)
    return;
  if (order == FTT_PRE_ORDER &&
      (flags == FTT_TRAVERSE_ALL ||
       ((flags & FTT_TRAVERSE_LEAFS) != 0 && FTT_CELL_IS_LEAF (cell)) ||
       ((flags & FTT_TRAVERSE_NON_LEAFS) != 0 && !FTT_CELL_IS_LEAF (cell))))
    (* func) (cell, s1, data);
  if (!FTT_CELL_IS_LEAF (cell)) {
    struct _FttOct * children = cell->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);

      if (!FTT_CELL_IS_DESTROYED (c))
	cell_traverse_cut (c, s1, order, flags, func, data, flatten);
    }
  }
  if (order == FTT_POST_ORDER &&
      (flags == FTT_TRAVERSE_ALL ||
       ((flags & FTT_TRAVERSE_LEAFS) != 0 && FTT_CELL_IS_LEAF (cell)) ||
       ((flags & FTT_TRAVERSE_NON_LEAFS) != 0 && !FTT_CELL_IS_LEAF (cell))))
    (* func) (cell, s1, data);
  if (s1 != s)
    gts_object_destroy (GTS_OBJECT (s1));
}

/**
 * gfs_cell_traverse_cut:
 * @root: the root #FttCell of the tree to traverse.
 * @s: a #GfsSurface.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 * 
 * Traverses a cell tree starting at the given root #FttCell. Calls
 * the given function for each cell cut by @s.
 */
void gfs_cell_traverse_cut (FttCell * root,
			    GfsSurface * s,
			    FttTraverseType order,
			    FttTraverseFlags flags,
			    FttCellTraverseCutFunc func,
			    gpointer data)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (func != NULL);

  cell_traverse_cut (root, s, order, flags, func, data, FALSE);
}

/**
 * gfs_cell_traverse_cut_2D:
 * @root: the root #FttCell of the tree to traverse.
 * @s: a #GfsSurface.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 * 
 * Traverses a cell tree starting at the given root #FttCell. Calls
 * the given function for each cell cut by @s.
 *
 * The cells are "flattened" in the z-direction.
 */
void gfs_cell_traverse_cut_2D (FttCell * root,
			       GfsSurface * s,
			       FttTraverseType order,
			       FttTraverseFlags flags,
			       FttCellTraverseCutFunc func,
			       gpointer data)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (func != NULL);

  cell_traverse_cut (root, s, order, flags, func, data, TRUE);
}
