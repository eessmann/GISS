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

#include <stdlib.h>
#include <math.h>

#include "config.h"
#include "boundary.h"
#include "simulation.h"
#include "adaptive.h"

static FttVector rpos[FTT_NEIGHBORS] = {
#if FTT_2D
  {1.,0.,0.}, {-1.,0.,0.}, {0.,1.,0.}, {0.,-1.,0.}
#else  /* FTT_3D */
  {1.,0.,0.}, {-1.,0.,0.}, {0.,1.,0.}, {0.,-1.,0.}, {0.,0.,1.}, {0.,0.,-1.}
#endif /* FTT_3D */
};

/* GfsBc: Object */

static void symmetry (FttCellFace * f, GfsBc * b)
{
  if (b->v->i == GFS_VELOCITY_INDEX (f->d/2) ||
      b->v->i == GFS_GRADIENT_INDEX (f->d/2))
    GFS_VARIABLE (f->cell, b->v->i) = - GFS_VARIABLE (f->neighbor, b->v->i);
  else
    GFS_VARIABLE (f->cell, b->v->i) =   GFS_VARIABLE (f->neighbor, b->v->i);
}

static void face_symmetry (FttCellFace * f, GfsBc * b)
{
  if (b->v->i == GFS_VELOCITY_INDEX (f->d/2) ||
      b->v->i == GFS_GRADIENT_INDEX (f->d/2))
    GFS_STATE (f->cell)->f[f->d].v = 
      GFS_STATE (f->neighbor)->f[FTT_OPPOSITE_DIRECTION (f->d)].v = 0.;
  else
    GFS_STATE (f->cell)->f[f->d].v = 
      GFS_STATE (f->neighbor)->f[FTT_OPPOSITE_DIRECTION (f->d)].v;
}

static void bc_write (GtsObject * o, FILE * fp)
{
  g_assert (GFS_BC (o)->v);
  fprintf (fp, "%s %s", o->klass->info.name, GFS_BC (o)->v->name);
}

static void bc_read (GtsObject ** o, GtsFile * fp)
{
  GfsBc * bc = GFS_BC (*o);

  g_assert (bc->b);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (klass)");
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  bc->v = gfs_variable_from_name (gfs_box_domain (bc->b->box)->variables, 
				  fp->token->str);
  if (bc->v == NULL)
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
  else
    gts_file_next_token (fp);
}

static void gfs_bc_class_init (GtsObjectClass * klass)
{
  klass->write = bc_write;
  klass->read =  bc_read;
}

static void gfs_bc_init (GfsBc * object)
{
  object->bc =             (FttFaceTraverseFunc) symmetry;
  object->homogeneous_bc = (FttFaceTraverseFunc) symmetry;
  object->face_bc =        (FttFaceTraverseFunc) face_symmetry;
}

GfsBcClass * gfs_bc_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_info = {
      "GfsBc",
      sizeof (GfsBc),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_class_init,
      (GtsObjectInitFunc) gfs_bc_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &gfs_bc_info);
  }

  return klass;
}

GfsBc * gfs_bc_new (GfsBcClass * k, GfsVariable * v, gboolean extra)
{
  GfsBc * b;

  g_return_val_if_fail (k != NULL, NULL);

  b = GFS_BC (gts_object_new (GTS_OBJECT_CLASS (k)));
  b->v = v;
  b->extra = extra;

  return b;
}

/* GfsBcValue: Object */

static void bc_value_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_bc_value_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_bc_value_class ())->parent_class->write) 
      (o, fp);
  if (GFS_BC_VALUE (o)->val)
    gfs_function_write (GFS_BC_VALUE (o)->val, fp);
}

static void bc_value_read (GtsObject ** o, GtsFile * fp)
{
  GfsBcValue * bc = GFS_BC_VALUE (*o);

  if (GTS_OBJECT_CLASS (gfs_bc_value_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_bc_value_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  if (bc->val == NULL)
    bc->val = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_read (GFS_BC_VALUE (*o)->val, fp);
}

static void bc_value_destroy (GtsObject * o)
{
  if (GFS_BC_VALUE (o)->val)
    gts_object_destroy (GTS_OBJECT (GFS_BC_VALUE (o)->val));

  (* GTS_OBJECT_CLASS (gfs_bc_value_class ())->parent_class->destroy) (o);
}

static void gfs_bc_value_class_init (GtsObjectClass * klass)
{
  klass->write   = bc_value_write;
  klass->read    = bc_value_read;
  klass->destroy = bc_value_destroy;
}

GfsBcClass * gfs_bc_value_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_value_info = {
      "GfsBcValue",
      sizeof (GfsBcValue),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) gfs_bc_value_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_class ()),
				  &gfs_bc_value_info);
  }

  return klass;
}

GfsBc * gfs_bc_value_new (GfsBcClass * k,
			  GfsVariable * v,
			  GfsFunction * val,
			  gboolean extra)
{
  GfsBcValue * bc = GFS_BC_VALUE (gfs_bc_new (k, v, extra));

  if (val == NULL)
    bc->val = gfs_function_new (gfs_function_class (), 0.);
  else
    bc->val = val;

  return GFS_BC (bc);
}

/* GfsBcDirichlet: Object */

static void dirichlet (FttCellFace * f, GfsBc * b)
{
  GFS_VARIABLE (f->cell, b->v->i) = 
    2.*gfs_function_face_value (GFS_BC_VALUE (b)->val, f, 
		     GFS_SIMULATION (gfs_box_domain (b->b->box))->time.t)
    - GFS_VARIABLE (f->neighbor, b->v->i);
}

static void homogeneous_dirichlet (FttCellFace * f, GfsBc * b)
{
  GFS_VARIABLE (f->cell, b->v->i) = - GFS_VARIABLE (f->neighbor, b->v->i);
}

static void face_dirichlet (FttCellFace * f, GfsBc * b)
{
  gdouble v = gfs_function_face_value (GFS_BC_VALUE (b)->val, f,
		   GFS_SIMULATION (gfs_box_domain (b->b->box))->time.t);
  
  GFS_STATE (f->cell)->f[f->d].v = 
    GFS_STATE (f->neighbor)->f[FTT_OPPOSITE_DIRECTION (f->d)].v = v;
}

static void gfs_bc_dirichlet_init (GfsBc * object)
{
  object->bc =             (FttFaceTraverseFunc) dirichlet;
  object->homogeneous_bc = (FttFaceTraverseFunc) homogeneous_dirichlet;
  object->face_bc =        (FttFaceTraverseFunc) face_dirichlet;
}

GfsBcClass * gfs_bc_dirichlet_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_dirichlet_info = {
      "GfsBcDirichlet",
      sizeof (GfsBcValue),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_bc_dirichlet_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_value_class ()),
				  &gfs_bc_dirichlet_info);
  }

  return klass;
}

/* GfsBcNeumann: Object */

static void neumann (FttCellFace * f, GfsBc * b)
{
  GFS_VARIABLE (f->cell, b->v->i) = 
    GFS_VARIABLE (f->neighbor, b->v->i) +
    gfs_function_face_value (GFS_BC_VALUE (b)->val, f, 
		  GFS_SIMULATION (gfs_box_domain (b->b->box))->time.t)
    *ftt_cell_size (f->cell);
}

static void homogeneous_neumann (FttCellFace * f, GfsBc * b)
{
  GFS_VARIABLE (f->cell, b->v->i) = GFS_VARIABLE (f->neighbor, b->v->i);
}

static void face_neumann (FttCellFace * f, GfsBc * b)
{
  gdouble v = 
    GFS_VARIABLE (f->neighbor, b->v->i) +
    gfs_function_face_value (GFS_BC_VALUE (b)->val, f, 
	         GFS_SIMULATION (gfs_box_domain (b->b->box))->time.t)
    *ftt_cell_size (f->cell)/2.;

  GFS_STATE (f->cell)->f[f->d].v = 
    GFS_STATE (f->neighbor)->f[FTT_OPPOSITE_DIRECTION (f->d)].v = v;
}

static void gfs_bc_neumann_init (GfsBc * object)
{
  object->bc =             (FttFaceTraverseFunc) neumann;
  object->homogeneous_bc = (FttFaceTraverseFunc) homogeneous_neumann;
  object->face_bc =        (FttFaceTraverseFunc) face_neumann;
}

GfsBcClass * gfs_bc_neumann_class (void)
{
  static GfsBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bc_neumann_info = {
      "GfsBcNeumann",
      sizeof (GfsBcValue),
      sizeof (GfsBcClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_bc_neumann_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_bc_value_class ()),
				  &gfs_bc_neumann_info);
  }

  return klass;
}

/* GfsBoundary: Object */

static void destroy_bc (GfsVariable * v, GtsObject * o)
{
  gts_object_destroy (o);
}

static void gfs_boundary_destroy (GtsObject * object)
{
  GfsBoundary * boundary = GFS_BOUNDARY (object);

  if (boundary->root)
    ftt_cell_destroy (boundary->root, 
		      (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
  boundary->box->neighbor[FTT_OPPOSITE_DIRECTION (boundary->d)] = NULL;

  gts_object_destroy (GTS_OBJECT (boundary->default_bc));
  if (boundary->bc) {
    g_hash_table_foreach (boundary->bc, (GHFunc) destroy_bc, NULL);
    g_hash_table_destroy (boundary->bc);
  }

  (* GTS_OBJECT_CLASS (gfs_boundary_class ())->parent_class->destroy) (object);
}

static void match (FttCell * cell, GfsBoundary * boundary)
{
  FttCell * neighbor = ftt_cell_neighbor (cell, boundary->d);
  guint level = ftt_cell_level (cell);

  cell->flags |= GFS_FLAG_BOUNDARY;
  if (neighbor == NULL || ftt_cell_level (neighbor) < level) {
    if (FTT_CELL_IS_ROOT (cell))
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	     "root cell is entirely outside of the fluid domain\n"
	     "the solid surface orientation may be incorrect");
    ftt_cell_destroy (cell, (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
    boundary->changed = TRUE;
    return;
  }
  if (ftt_cell_level (neighbor) == level) {
    GfsSolidVector * s = GFS_STATE (neighbor)->solid;

    if (s && s->s[FTT_OPPOSITE_DIRECTION (boundary->d)] == 0.) {
      if (FTT_CELL_IS_ROOT (cell))
	g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	       "root cell is entirely outside of the fluid domain\n"
	       "the solid surface orientation may be incorrect");
      ftt_cell_destroy (cell, (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
      boundary->changed = TRUE;
      return;
    }
    if (s) {
      FttDirection d;
      FttComponent c;
      GfsSolidVector * t;

      if (GFS_STATE (cell)->solid == NULL)
	GFS_STATE (cell)->solid = g_malloc0 (sizeof (GfsSolidVector));
      t = GFS_STATE (cell)->solid;
      t->a = s->a;
      for (d = 0; d < FTT_NEIGHBORS; d++)
	if (d/2 == boundary->d/2)
	  t->s[d] = s->s[FTT_OPPOSITE_DIRECTION (d)];
	else
	  t->s[d] = s->s[d];
      for (c = 0; c < FTT_DIMENSION; c++)
	if (c == boundary->d/2) {
	  FttVector p1, p2;
	  ftt_cell_pos (cell, &p1);
	  ftt_cell_pos (cell, &p2);
	  (&t->cm.x)[c] = (&p1.x)[c] + (&p2.x)[c] - (&s->cm.x)[c];
	  (&t->ca.x)[c] = (&p1.x)[c] + (&p2.x)[c] - (&s->ca.x)[c];
	}
	else {
	  (&t->cm.x)[c] = (&s->cm.x)[c];
	  (&t->ca.x)[c] = (&s->ca.x)[c];
	}
    }
    else if (GFS_STATE (cell)->solid != NULL) {
      g_free (GFS_STATE (cell)->solid);
      GFS_STATE (cell)->solid = NULL;
    }      
    if (FTT_CELL_IS_LEAF (cell) && !FTT_CELL_IS_LEAF (neighbor)) {
      ftt_cell_refine_single (cell, (FttCellInitFunc) gfs_cell_fine_init,
			      gfs_box_domain (boundary->box));
      boundary->changed = TRUE;
    }
  }
  else
    g_assert_not_reached ();
  if (!FTT_CELL_IS_LEAF (cell))
    level++;
  if (level > boundary->depth)
    boundary->depth = level;
}

static void boundary_match (GfsBoundary * boundary)
{
  guint l = ftt_cell_level (boundary->root);

  boundary->changed = FALSE;
  boundary->depth = l;
  while (l <= boundary->depth) {
    ftt_cell_traverse_boundary (boundary->root, boundary->d,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, l,
				(FttCellTraverseFunc) match, boundary);
    l++;
  }
  if (boundary->changed)
    ftt_cell_flatten (boundary->root, boundary->d, (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
}

static void is_extra (GfsVariable * v, GfsBc * bc, gboolean * extra)
{
  if (bc->extra)
    *extra = TRUE;
}

static void write_extra (GfsVariable * v, GfsBc * bc, FILE * fp)
{
  if (bc->extra) {
    g_assert (GTS_OBJECT (bc)->klass->write);
    (* GTS_OBJECT (bc)->klass->write) (GTS_OBJECT (bc), fp);
    fputc ('\n', fp);
  }
}

static void gfs_boundary_write (GtsObject * o, FILE * fp)
{
  GfsBoundary * b = GFS_BOUNDARY (o);
  gboolean any_extra = FALSE;

  if (GTS_OBJECT_CLASS (gfs_boundary_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_boundary_class ())->parent_class->write) 
      (o, fp);

  g_hash_table_foreach (b->bc, (GHFunc) is_extra, &any_extra);
  if (any_extra) {
    fputs (" {\n", fp);
    g_hash_table_foreach (b->bc, (GHFunc) write_extra, fp);
    fputc ('}', fp);
  }
}

static gboolean boundary_read_extra_bc (GfsBoundary * b, GtsFile * fp)
{
  gboolean ret = FALSE;

  if (fp->type != '{')
    return ret;

  fp->scope_max++;
  gts_file_next_token (fp);

  while (fp->type != GTS_ERROR && fp->type != '}') {
    if (fp->type == '\n') {
      gts_file_next_token (fp);
      continue;
    }
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a keyword");
      return ret;
    }
    else {
      GtsObjectClass * klass = gfs_object_class_from_name (fp->token->str);
      GtsObject * object;
      
      if (klass == NULL) {
	gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
	return ret;
      }
      else if (!gts_object_class_is_from_class (klass, gfs_bc_class ())) {
	gts_file_error (fp, "`%s' is not a GfsBc", fp->token->str);
	return ret;
      }

      object = gts_object_new (klass);
      g_assert (klass->read);
      GFS_BC (object)->b = b;
      GFS_BC (object)->extra = TRUE;
      (* klass->read) (&object, fp);
      if (fp->type == GTS_ERROR) {
	gts_object_destroy (object);
	return ret;
      }

      gfs_boundary_add_bc (b, GFS_BC (object));
      ret = TRUE;
    }
  }
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return ret;
  }
  fp->scope_max--;
  gts_file_next_token (fp);

  return ret;
}

static void gfs_boundary_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_boundary_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_boundary_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  boundary_read_extra_bc (GFS_BOUNDARY (*o), fp);
}

static void gfs_boundary_class_init (GfsBoundaryClass * klass)
{
  klass->match = boundary_match;

  GTS_OBJECT_CLASS (klass)->write =   gfs_boundary_write;
  GTS_OBJECT_CLASS (klass)->read =    gfs_boundary_read;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_boundary_destroy;
}

static void gfs_boundary_init (GfsBoundary * b)
{
  b->type = GFS_BOUNDARY_CENTER_VARIABLE;
  b->bc = g_hash_table_new (g_str_hash, g_str_equal);
  gfs_boundary_set_default_bc (b, gfs_bc_new (gfs_bc_class (), NULL, FALSE));
}

GfsBoundaryClass * gfs_boundary_class (void)
{
  static GfsBoundaryClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_boundary_info = {
      "GfsBoundary",
      sizeof (GfsBoundary),
      sizeof (GfsBoundaryClass),
      (GtsObjectClassInitFunc) gfs_boundary_class_init,
      (GtsObjectInitFunc) gfs_boundary_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (gts_object_class (), &gfs_boundary_info);
  }

  return klass;
}

/**
 * gfs_boundary_new:
 * @klass: a #GfsBoundaryClass.
 * @box: a #GfsBox.
 * @d: a direction.
 *
 * Creates a new boundary of type @klass for @box in direction @d.
 *
 * This function fails if @box has already a boundary in direction @d.
 *
 * Returns: a new #GfsBoundary.
 */
GfsBoundary * gfs_boundary_new (GfsBoundaryClass * klass,
				GfsBox * box,
				FttDirection d)
{
  GfsBoundary * boundary;
  GfsDomain * domain;
  FttVector pos;
  gdouble size;

  g_return_val_if_fail (box != NULL, NULL);
  g_return_val_if_fail (d < FTT_NEIGHBORS, NULL);
  g_return_val_if_fail (box->neighbor[d] == NULL, NULL);

  boundary = GFS_BOUNDARY (gts_object_new (GTS_OBJECT_CLASS (klass)));
  boundary->box = box;
  box->neighbor[d] = GTS_OBJECT (boundary);
  boundary->d = FTT_OPPOSITE_DIRECTION (d);
  if (box->root) {
    domain = gfs_box_domain (box);
    boundary->root = ftt_cell_new ((FttCellInitFunc) gfs_cell_init, domain);
    ftt_cell_set_level (boundary->root, ftt_cell_level (box->root));
    ftt_cell_set_neighbor_match (boundary->root, box->root, boundary->d, 
				 (FttCellInitFunc) gfs_cell_init, domain);
    ftt_cell_pos (box->root, &pos);
    size = ftt_cell_size (box->root);
    pos.x += rpos[d].x*size;
    pos.y += rpos[d].y*size;
    pos.z += rpos[d].z*size;
    ftt_cell_set_pos (boundary->root, &pos);

    boundary_match (boundary);
  }

  return boundary;
}

/**
 * gfs_boundary_send:
 * @boundary: a #GfsBoundary.
 *
 * Calls the @send() method of @boundary.
 */
void gfs_boundary_send (GfsBoundary * boundary)
{
  g_return_if_fail (boundary != NULL);

  if (GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->send)
    (* GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->send) (boundary);
}

/**
 * gfs_boundary_receive:
 * @boundary: a #GfsBoundary.
 * @flags: the traversal flags.
 * @max_depth: the maximum depth of the traversal.
 *
 * Calls the @receive() method of @boundary.
 */
void gfs_boundary_receive (GfsBoundary * boundary,
			   FttTraverseFlags flags,
			   gint max_depth)
{
  g_return_if_fail (boundary != NULL);

  if (GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->receive)
    (* GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->receive)
      (boundary, flags, max_depth);
}

/**
 * gfs_boundary_synchronize:
 * @boundary: a #GfsBoundary.
 *
 * Calls the @synchronize() method of @boundary.
 */
void gfs_boundary_synchronize (GfsBoundary * boundary)
{
  g_return_if_fail (boundary != NULL);

  if (GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->synchronize)
    (* GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass)->synchronize) (boundary);
}

/**
 * gfs_boundary_lookup:
 * @b: a #GfsBoundary.
 * @v: a #GfsVariable.
 *
 * Returns: the #GfsBc associated with @b and @v.
 */
GfsBc * gfs_boundary_lookup_bc (GfsBoundary * b, GfsVariable * v)
{
  GfsBc * bv;

  g_return_val_if_fail (b != NULL, NULL);
  g_return_val_if_fail (v != NULL, NULL);

  if (!v->name || !(bv = g_hash_table_lookup (b->bc, v->name))) {
    bv = b->default_bc;
    bv->v = v;
  }
  return bv;
}

/**
 * gfs_boundary_set_default_bc:
 * @b: a #GfsBoundary.
 * @bc: a #GfsBc.
 *
 * Sets the default boundary condition for @b to @bc.
 */
void gfs_boundary_set_default_bc (GfsBoundary * b, GfsBc * bc)
{
  g_return_if_fail (b != NULL);
  g_return_if_fail (bc != NULL);
  g_return_if_fail (bc->b == NULL || bc->b == b);

  if (b->default_bc)
    gts_object_destroy (GTS_OBJECT (b->default_bc));
  b->default_bc = bc;
  bc->b = b;
}

/**
 * gfs_boundary_add_bc:
 * @b: a #GfsBoundary.
 * @bc: a #GfsBc.
 *
 * Adds boundary condition @bc to @b.
 */
void gfs_boundary_add_bc (GfsBoundary * b, GfsBc * bc)
{
  GfsBc * old;
  
  g_return_if_fail (b != NULL);
  g_return_if_fail (bc != NULL);
  g_return_if_fail (bc->v != NULL);
  g_return_if_fail (bc->v->name != NULL);
  g_return_if_fail (bc->b == NULL || bc->b == b);

  old = g_hash_table_lookup (b->bc, bc->v->name);
  if (!old || !old->extra) {
    if (old) gts_object_destroy (GTS_OBJECT (old));
    g_hash_table_insert (b->bc, bc->v->name, bc);
    bc->b = b;
  }
  else
    gts_object_destroy (GTS_OBJECT (bc));
}

/* GfsBoundaryInflowConstant: Object */

static GtsColor inflow_color (GtsObject * o)
{
  GtsColor c = { 0., 0., 1. }; /* blue */

  return c;
}

static void inflow_constant_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_boundary_inflow_constant_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_boundary_inflow_constant_class ())->parent_class->write) 
      (o, fp);

  gfs_function_write (GFS_BOUNDARY_INFLOW_CONSTANT (o)->un, fp);
}

static void inflow_constant_read (GtsObject ** o, GtsFile * fp)
{
  GfsBoundary * b = GFS_BOUNDARY (*o);  
  FttComponent c;
  GfsVariable * v;
  GfsFunction * un = GFS_BOUNDARY_INFLOW_CONSTANT (*o)->un;

  if (GTS_OBJECT_CLASS (gfs_boundary_inflow_constant_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_boundary_inflow_constant_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  gfs_function_read (un, fp);

  v = gfs_variable_from_name (gfs_box_domain (b->box)->variables, "U");
  for (c = 0; c < FTT_DIMENSION; c++, v = v->next)
    if (c == b->d/2)
      gfs_boundary_add_bc (b, gfs_bc_value_new (gfs_bc_dirichlet_class (),
						v, un, FALSE));
    else
      gfs_boundary_add_bc (b, gfs_bc_value_new (gfs_bc_dirichlet_class (),
						v, NULL, FALSE));
}

static void gfs_boundary_inflow_constant_class_init (GtsObjectClass * klass)
{
  klass->read    = inflow_constant_read;
  klass->write   = inflow_constant_write;
  klass->color   = inflow_color;
}

static void gfs_boundary_inflow_constant_init 
  (GfsBoundaryInflowConstant * object)
{
  object->un = gfs_function_new (gfs_function_class (), 0.);
}

GfsBoundaryInflowConstantClass * gfs_boundary_inflow_constant_class (void)
{
  static GfsBoundaryInflowConstantClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_boundary_inflow_constant_info = {
      "GfsBoundaryInflowConstant",
      sizeof (GfsBoundaryInflowConstant),
      sizeof (GfsBoundaryInflowConstantClass),
      (GtsObjectClassInitFunc) gfs_boundary_inflow_constant_class_init,
      (GtsObjectInitFunc) gfs_boundary_inflow_constant_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_boundary_class ()),
				  &gfs_boundary_inflow_constant_info);
  }

  return klass;
}

/* GfsBoundaryOutflow: Object */

static GtsColor outflow_color (GtsObject * o)
{
  GtsColor c = { 0., 1., 0. }; /* green */

  return c;
}

static void outflow_read (GtsObject ** o, GtsFile * fp)
{
  GfsBoundary * b = GFS_BOUNDARY (*o);  
  FttComponent c;
  GfsVariable * v;

  if (GTS_OBJECT_CLASS (gfs_boundary_outflow_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_boundary_outflow_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  v = gfs_variable_from_name (gfs_box_domain (b->box)->variables, "U");
  for (c = 0; c < b->d/2; c++, v = v->next) 
    ;

  gfs_boundary_add_bc (b, gfs_bc_value_new (gfs_bc_neumann_class (), v, NULL,
					    FALSE));
  gfs_boundary_add_bc (b, gfs_bc_value_new (gfs_bc_dirichlet_class (),
					    gfs_p, NULL, FALSE));
}

static void gfs_boundary_outflow_class_init (GfsBoundaryClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read   = outflow_read;
  GTS_OBJECT_CLASS (klass)->color  = outflow_color;
}

GfsBoundaryOutflowClass * gfs_boundary_outflow_class (void)
{
  static GfsBoundaryOutflowClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_boundary_outflow_info = {
      "GfsBoundaryOutflow",
      sizeof (GfsBoundary),
      sizeof (GfsBoundaryOutflowClass),
      (GtsObjectClassInitFunc) gfs_boundary_outflow_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_boundary_class ()),
				  &gfs_boundary_outflow_info);
  }

  return klass;
}

/* GfsGEdge: Object */

static void gfs_gedge_write (GtsObject * object, FILE * fp)
{
  fprintf (fp, " %s", ftt_direction_name [GFS_GEDGE (object)->d]);
}

static void gfs_gedge_read (GtsObject ** o, GtsFile * fp)
{
  GfsGEdge * e = GFS_GEDGE (*o);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (direction)");
    return;
  }
  e->d = ftt_direction_from_name (fp->token->str);
  if (e->d >= FTT_NEIGHBORS) {
    gts_file_error (fp, "unknown direction `%s'", fp->token->str);
    e->d = 0;
    return;
  }
  gts_file_next_token (fp);
}

static void gfs_gedge_class_init (GtsObjectClass * klass)
{
  klass->write = gfs_gedge_write;
  klass->read = gfs_gedge_read;
}

static void gfs_gedge_init (GfsGEdge * object)
{
  object->d = -1;
}

GfsGEdgeClass * gfs_gedge_class (void)
{
  static GfsGEdgeClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gedge_info = {
      "GfsGEdge",
      sizeof (GfsGEdge),
      sizeof (GfsGEdgeClass),
      (GtsObjectClassInitFunc) gfs_gedge_class_init,
      (GtsObjectInitFunc) gfs_gedge_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_gedge_class ()),
				  &gfs_gedge_info);
  }

  return klass;
}

/**
 * gfs_gedge_link_boxes:
 * @edge: a #GfsGEdge.
 *
 * Links the two boxes connected by @edge. The boxes are set as their
 * respective neighbors in the direction defined by @edge and their
 * relative positions are set accordingly.
 */
void gfs_gedge_link_boxes (GfsGEdge * edge)
{
  GfsBox * b1, * b2;

  g_return_if_fail (edge != NULL);
  g_return_if_fail (GTS_GEDGE (edge)->n1 != NULL);
  g_return_if_fail (GTS_GEDGE (edge)->n2 != NULL);
  g_return_if_fail (edge->d >= 0 && edge->d < FTT_NEIGHBORS);

  b1 = GFS_BOX (GTS_GEDGE (edge)->n1);
  b2 = GFS_BOX (GTS_GEDGE (edge)->n2);

  g_return_if_fail (b1->neighbor[edge->d] == NULL);
  g_return_if_fail (b2->neighbor[FTT_OPPOSITE_DIRECTION (edge->d)] == NULL);

  ftt_cell_set_neighbor (b1->root, b2->root, edge->d, 
			 (FttCellInitFunc) gfs_cell_init, gfs_box_domain (b1));
  b1->neighbor[edge->d] = GTS_OBJECT (b2);
  b2->neighbor[FTT_OPPOSITE_DIRECTION (edge->d)] = GTS_OBJECT (b1);
  if (b1 != b2)
    gfs_box_set_relative_pos (b2, b1, edge->d);
}

/**
 * gfs_gedge_new:
 * @klass: a #GfsGEdgeClass.
 * @b1: a #GfsBox.
 * @b2: another #GfsBox.
 * @d: a direction.
 *
 * Returns: a new #GfsGEdge linking @b1 to @b2 in direction @d. The
 * boxes are linked using gfs_gedge_link_boxes().
 */
GfsGEdge * gfs_gedge_new (GfsGEdgeClass * klass,
			  GfsBox * b1, GfsBox * b2,
			  FttDirection d)
{
  GfsGEdge * edge;

  g_return_val_if_fail (klass != NULL, NULL);
  g_return_val_if_fail (b1 != NULL, NULL);
  g_return_val_if_fail (b2 != NULL, NULL);
  g_return_val_if_fail (d >= 0 && d < FTT_NEIGHBORS, NULL);

  edge = GFS_GEDGE (gts_gedge_new (GTS_GEDGE_CLASS (klass),
				  GTS_GNODE (b1), GTS_GNODE (b2)));
  edge->d = d;
  
  gfs_gedge_link_boxes (edge);

  return edge;
}

/* GfsBox: Object */

static void gfs_box_destroy (GtsObject * object)
{
  GfsBox * box = GFS_BOX (object);
  FttDirection d;

  if (box->root)
    ftt_cell_destroy (box->root, (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d]))
      gts_object_destroy (box->neighbor[d]);
    else if (GFS_IS_BOX (box->neighbor[d])) {
      g_assert (GFS_BOX (box->neighbor[d])->neighbor[FTT_OPPOSITE_DIRECTION (d)] == GTS_OBJECT (box));
      GFS_BOX (box->neighbor[d])->neighbor[FTT_OPPOSITE_DIRECTION (d)] = NULL;
    }

  (* GTS_OBJECT_CLASS (gfs_box_class ())->parent_class->destroy) (object);
}

static void box_size (FttCell * cell, guint * size)
{
  (*size)++;
}

static void gfs_box_write (GtsObject * object, FILE * fp)
{
  GfsBox * box = GFS_BOX (object);
  FttDirection d;
  guint size = 0;
  GfsDomain * domain = gfs_box_domain (box);

  ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) box_size, &size);
  fprintf (fp, "%s { id = %u pid = %d size = %u", 
	   object->klass->info.name, box->id, box->pid, size);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      fprintf (fp, " %s = %s",
	       ftt_direction_name[d],
	       box->neighbor[d]->klass->info.name);
      if (box->neighbor[d]->klass->write)
	(* box->neighbor[d]->klass->write) (box->neighbor[d], fp);
    }
  fputs (" }", fp);
  if (domain != NULL && domain->max_depth_write > -2) {
    fputs (" {\n", fp);
    if (domain->binary)
      ftt_cell_write_binary (box->root, domain->max_depth_write, fp, 
			     (FttCellWriteFunc) gfs_cell_write_binary, domain->variables_io);
    else
      ftt_cell_write (box->root, domain->max_depth_write, fp, 
		      (FttCellWriteFunc) gfs_cell_write, domain->variables_io);
    fputc ('}', fp);
  }
}

static void gfs_box_read (GtsObject ** o, GtsFile * fp)
{
  GfsBox * b = GFS_BOX (*o);
  GtsObjectClass * klass;
  gboolean class_changed = FALSE;
  GtsFileVariable var[] = {
    {GTS_UINT, "id",     TRUE},
    {GTS_INT,  "pid",    TRUE},
    {GTS_UINT, "size",   TRUE},
    {GTS_FILE, "right",  TRUE},
    {GTS_FILE, "left",   TRUE},
    {GTS_FILE, "top",    TRUE},
    {GTS_FILE, "bottom", TRUE},
#if (!FTT_2D)
    {GTS_FILE, "front",  TRUE},
    {GTS_FILE, "back",   TRUE},
#endif /* 3D */
    {GTS_NONE}
  };
  GtsFileVariable * v;
  gfloat weight;
  GfsDomain * domain;

  g_assert (GTS_SLIST_CONTAINEE (b)->containers &&
	    !GTS_SLIST_CONTAINEE (b)->containers->next);
  domain = GFS_DOMAIN (GTS_SLIST_CONTAINEE (b)->containers->data);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsBoxClass)");
    return;
  }
  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_box_class ())) {
    gts_file_error (fp, "`%s' is not a GfsBox", fp->token->str);
    return;
  }
  if (klass != (*o)->klass) {
    *o = gts_object_new (klass);
    gts_object_destroy (GTS_OBJECT (b));
    b = GFS_BOX (*o);
    gts_container_add (GTS_CONTAINER (domain), GTS_CONTAINEE (b));
    class_changed = TRUE;
  }
  gts_file_next_token (fp);

  g_assert (b->root == NULL);
  b->root = ftt_cell_new ((FttCellInitFunc) gfs_cell_init, domain);

  weight = gts_gnode_weight (GTS_GNODE (b));
  var[0].data = &b->id;
  var[1].data = &b->pid;
  var[2].data = &b->size;
  gts_file_assign_start (fp, var);
  if (fp->type == GTS_ERROR)
    return;
  while ((v = gts_file_assign_next (fp, var)))
    if (v->type == GTS_FILE) {
      GtsObjectClass * boundary_class = gfs_object_class_from_name (fp->token->str);
      GtsObject * boundary;
	
      if (boundary_class == NULL) {
	gts_file_error (fp, "unknown class `%s'", fp->token->str);
	return;
      }
      if (!gts_object_class_is_from_class (boundary_class, gfs_boundary_class ())) {
	gts_file_error (fp, "`%s' is not a GfsBoundary", fp->token->str);
	return;
      }
      boundary = GTS_OBJECT (gfs_boundary_new (GFS_BOUNDARY_CLASS (boundary_class),
					       b, ftt_direction_from_name (v->name)));
      gts_file_next_token (fp);
      if (boundary_class->read)
	(* boundary_class->read) (&boundary, fp);
    }
  
  if (fp->type == '{') {
    FttDirection d;

    ftt_cell_destroy (b->root, (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
    fp->scope_max++;
    if (domain->binary) {
      if (gts_file_getc (fp) != '\n') {
      	gts_file_error (fp, "expecting a newline");
      	return;
      }
      b->root = ftt_cell_read_binary (fp, (FttCellReadFunc) gfs_cell_read_binary, domain);
      if (fp->type == GTS_ERROR)
	return;
      gts_file_next_token (fp);
    }
    else {
      gts_file_first_token_after (fp, '\n');
      b->root = ftt_cell_read (fp, (FttCellReadFunc) gfs_cell_read, domain);
    }
    fp->scope_max--;
    if (fp->type == GTS_ERROR)
      return;
    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    gts_file_next_token (fp);

    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (GFS_IS_BOUNDARY (b->neighbor[d])) {
	GfsBoundary * boundary = GFS_BOUNDARY (b->neighbor[d]);

	ftt_cell_set_neighbor_match (boundary->root, b->root, boundary->d, 
				     (FttCellInitFunc) gfs_cell_init, domain);
      }
  }

  if (ftt_cell_level (b->root) != domain->rootlevel) {
    FttDirection d;

    ftt_cell_set_level (b->root, domain->rootlevel);
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (GFS_IS_BOUNDARY (b->neighbor[d]))
	ftt_cell_set_level (GFS_BOUNDARY (b->neighbor[d])->root, domain->rootlevel);
  }
  /* updates weight of domain */
  GTS_WGRAPH (domain)->weight += gts_gnode_weight (GTS_GNODE (b)) - weight;

  if (class_changed && klass->read)
    (* klass->read) (o, fp);
}

static gfloat gfs_box_weight (GtsGNode * node)
{
  GfsBox * box = GFS_BOX (node);

  if (box->size >= 0)
    return box->size;
  else {
    guint size = 0;

    if (box->root)
      ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			 (FttCellTraverseFunc) box_size, &size);
    return size;
  }
}

static void gfs_box_class_init (GfsBoxClass * klass)
{
  GTS_GNODE_CLASS (klass)->weight = gfs_box_weight;

  GTS_OBJECT_CLASS (klass)->destroy = gfs_box_destroy;
  GTS_OBJECT_CLASS (klass)->write = gfs_box_write;
  GTS_OBJECT_CLASS (klass)->read = gfs_box_read;
}

static void gfs_box_init (GfsBox * box)
{
  static guint id = 1;

  box->id = id++;
  box->pid = -1;
  box->size = -1;
}

GfsBoxClass * gfs_box_class (void)
{
  static GfsBoxClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_box_info = {
      "GfsBox",
      sizeof (GfsBox),
      sizeof (GfsBoxClass),
      (GtsObjectClassInitFunc) gfs_box_class_init,
      (GtsObjectInitFunc) gfs_box_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_gnode_class ()),
				  &gfs_box_info);
  }

  return klass;
}

GfsBox * gfs_box_new (GfsBoxClass * klass)
{
  GfsBox * object;

  object = GFS_BOX (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

/* GfsBoxNotAdapt: Object */

static void gfs_box_not_adapt_read (GtsObject ** o, GtsFile * fp)
{
  GfsSimulation * sim;
  GfsAdaptNotBox * a;
  GfsBoxNotAdapt * b = GFS_BOX_NOT_ADAPT (*o);
  GfsDomain * domain = gfs_box_domain (GFS_BOX (b));

  g_assert (GFS_IS_SIMULATION (domain));
  sim = GFS_SIMULATION (domain);

  g_assert (gts_container_size (GTS_CONTAINER (b->c)) == 0);
  a = gfs_adapt_not_box_new (gfs_adapt_not_box_class (), GFS_BOX (b));
  gts_container_add (GTS_CONTAINER (b->c), GTS_CONTAINEE (a));
  gts_container_add (GTS_CONTAINER (sim->adapts), GTS_CONTAINEE (a));
}

static void gfs_box_not_adapt_destroy (GtsObject * object)
{
  GfsBoxNotAdapt * b = GFS_BOX_NOT_ADAPT (object);

  gts_container_foreach (GTS_CONTAINER (b->c),
			 (GtsFunc) gts_object_destroy, NULL);
  gts_object_destroy (GTS_OBJECT (b->c));

  (* GTS_OBJECT_CLASS (gfs_box_not_adapt_class ())->parent_class->destroy) 
    (object);
}

static void gfs_box_not_adapt_class_init (GfsBoxClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_box_not_adapt_read;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_box_not_adapt_destroy;
}

static void gfs_box_not_adapt_init (GfsBoxNotAdapt * object)
{
  object->c = GTS_SLIST_CONTAINER (gts_container_new 
				   (GTS_CONTAINER_CLASS 
				    (gts_slist_container_class ())));
}

GfsBoxClass * gfs_box_not_adapt_class (void)
{
  static GfsBoxClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_box_not_adapt_info = {
      "GfsBoxNotAdapt",
      sizeof (GfsBoxNotAdapt),
      sizeof (GfsBoxClass),
      (GtsObjectClassInitFunc) gfs_box_not_adapt_class_init,
      (GtsObjectInitFunc) gfs_box_not_adapt_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_box_class ()),
				  &gfs_box_not_adapt_info);
  }

  return klass;
}

static void box_set_pos (GfsBox * box, FttVector * pos, 
			 GHashTable * set,
			 FttDirection dold)
{
  FttVector p;
  static FttDirection id[FTT_NEIGHBORS][FTT_NEIGHBORS] = 
#if FTT_2D
  {
    {0,1,2,3},
    {1,0,3,2},
    {2,3,1,0},
    {3,2,0,1},
  };
#else  /* 3D */
  {
    {0,1,2,3,5,4},
    {1,0,3,2,4,5},
    {2,3,1,0,5,4},
    {3,2,0,1,4,5},
    {4,5,2,3,0,1},
    {5,4,3,2,1,0}
  };
#endif /* 3D */
  FttDirection i;
  gdouble size;

  if (g_hash_table_lookup (set, box))
    return;
  g_hash_table_insert (set, box, box);

  size = ftt_cell_size (box->root);
  ftt_cell_set_pos (box->root, pos);
  for (i = 0; i < FTT_NEIGHBORS; i++) {
    FttDirection d = id[dold][i];
    
    p.x = pos->x + rpos[d].x*size;
    p.y = pos->y + rpos[d].y*size;
    p.z = pos->z + rpos[d].z*size;
    if (GFS_IS_BOX (box->neighbor[d]))
      box_set_pos (GFS_BOX (box->neighbor[d]), &p, set, d);
    else if (GFS_IS_BOUNDARY (box->neighbor[d]))
      ftt_cell_set_pos (GFS_BOUNDARY (box->neighbor[d])->root, &p);
  }
}

/**
 * gfs_box_set_pos:
 * @box: a #GfsBox.
 * @pos: the new position of the center of the box.
 *
 * Recursively sets the position of the center of @box and of its
 * neighbors.  
 */
void gfs_box_set_pos (GfsBox * box, FttVector * pos)
{
  GHashTable * set;

  g_return_if_fail (box != NULL);
  g_return_if_fail (pos != NULL);

  set = g_hash_table_new (NULL, NULL);
  box_set_pos (box, pos, set, FTT_RIGHT);
  g_hash_table_destroy (set);
}

/**
 * gfs_box_set_relative_pos:
 * @box: a #GfsBox.
 * @reference: a reference #GfsBox.
 * @d: the direction in which @box is found relative to @reference.
 *
 * Recursively sets the position of the center of @box and of its
 * neighbors relative to the position of @reference in direction @d.
 */
void gfs_box_set_relative_pos (GfsBox * box, GfsBox * reference, FttDirection d)
{
  FttVector pos;
  gdouble size;

  g_return_if_fail (box != NULL);
  g_return_if_fail (reference != NULL);
  g_return_if_fail (d >= 0 && d < FTT_NEIGHBORS);

  ftt_cell_pos (reference->root, &pos);
  size = ftt_cell_size (reference->root);
  pos.x += rpos[d].x*size;
  pos.y += rpos[d].y*size;
  pos.z += rpos[d].z*size;
  gfs_box_set_pos (box, &pos);
}

#ifndef G_CAN_INLINE
/**
 * gfs_box_domain:
 * @box: a #GfsBox.
 *
 * Returns: the #GfsDomain to which @box belongs or %NULL if @box does not
 * belong to any domain.
 */
GfsDomain * gfs_box_domain (GfsBox * box)
{
  GfsDomain * d;

  g_return_val_if_fail (box != NULL, NULL);
  
  d = GTS_OBJECT (box)->reserved;
  if (GTS_SLIST_CONTAINEE (box)->containers) {
    GSList * i = GTS_SLIST_CONTAINEE (box)->containers;

    while (i->next)
      i = i->next;
    d = i->data;
  }
  g_assert (GFS_IS_DOMAIN (d));
  return d;
}
#endif /* not G_CAN_INLINE */
