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

#include "surface.h"
#include "simulation.h"

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

  if (fp->type == '{') {
    fp->scope_max++;
    gts_file_next_token (fp);
    g_assert (!surface->s);
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
  else {
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
    GtsVector r = {0.,0.,0.}, s = {1.,1.,1.}, t = {0.,0.,0.};
    gdouble angle = 0., scale = 1.;
    gboolean flip = FALSE;
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
      {GTS_DOUBLE, "angle", TRUE},
      {GTS_INT,  "flip", TRUE},
      {GTS_NONE}
    };
    GtsFileVariable * v = var;

    (v++)->data = &r[0];
    (v++)->data = &r[1];
    (v++)->data = &r[2];

    (v++)->data = &s[0];
    (v++)->data = &s[1];
    (v++)->data = &s[2];

    (v++)->data = &t[0];
    (v++)->data = &t[1];
    (v++)->data = &t[2];

    (v++)->data = &scale;
    (v++)->data = &angle;

    (v++)->data = &flip;

    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;

    if (var[9].set)
      s[0] = s[1] = s[2] = scale;
    if (var[10].set && gts_vector_norm (r) == 0.) {
      gts_file_variable_error (fp, var, "angle",
			       "a non-zero rotation vector must be specified");
      return;
    }
    
    GtsMatrix * m = gts_matrix_translate (NULL, t);
    if (angle != 0.) {
      GtsMatrix * mr = gts_matrix_rotate (NULL, r, angle*M_PI/180.);
      GtsMatrix * m1 = gts_matrix_product (m, mr);
      gts_matrix_destroy (m);
      gts_matrix_destroy (mr);
      m = m1;
    }
    GtsMatrix * ms = gts_matrix_scale (NULL, s);
    GtsMatrix * M = gts_matrix_product (m, ms);
    gts_matrix_destroy (m);
    gts_matrix_destroy (ms);
    gts_surface_foreach_vertex (surface->s, (GtsFunc) gts_point_transform, M);
    gts_matrix_destroy (M);

    if (flip)
      gts_surface_foreach_face (surface->s, (GtsFunc) gts_triangle_revert, NULL);
  }
}

static void surface_write (GtsObject * o, FILE * fp)
{
  fputs (" { ", fp);
  GtsSurface * s = GFS_SURFACE (o)->s;
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

static void surface_destroy (GtsObject * object)
{
  if (GFS_SURFACE (object)->s)
    gts_object_destroy (GTS_OBJECT (GFS_SURFACE (object)->s));

  (* GTS_OBJECT_CLASS (gfs_surface_class ())->parent_class->destroy) (object);
}

static void gfs_surface_class_init (GtsObjectClass * klass)
{
  klass->read = surface_read;
  klass->write = surface_write;
  klass->destroy = surface_destroy;
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
      (GtsObjectInitFunc) NULL,
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
 * @fp: a #GtsFile.
 * 
 * Calls the read() method of @s.
 */
void gfs_surface_read (GfsSurface * s, GtsFile * fp)
{
  GtsObject * o = (GtsObject *) s;

  g_return_if_fail (s != NULL);
  g_return_if_fail (fp != NULL);

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
 * Returns: a new #GfsSurface containing a subset of @s which may
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

  g_assert_not_implemented ();
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
