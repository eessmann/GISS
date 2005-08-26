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

#ifndef __FTT_H__
#define __FTT_H__

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <gts.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "gfsconfig.h"

#define FTT_MAINTAINER "s.popinet@niwa.cri.nz"

#if (FTT_2D || FTT_2D3)
# define FTT_CELLS     4
#else  /* FTT_3D */
# define FTT_CELLS     8
#endif /* FTT_3D */

typedef struct _FttCell          FttCell;
typedef struct _FttCellFace      FttCellFace;
typedef struct _FttCellNeighbors FttCellNeighbors;
typedef struct _FttCellChildren  FttCellChildren;

typedef struct _FttVector        FttVector;

struct _FttVector {
  gdouble x, y, z;
};

#if FTT_2D
# define ftt_vector_norm(v) (sqrt((v)->x*(v)->x + (v)->y*(v)->y))
#else  /* 3D */
# define ftt_vector_norm(v) (sqrt((v)->x*(v)->x + (v)->y*(v)->y + (v)->z*(v)->z))
#endif /* 3D */

typedef enum
{
  FTT_TRAVERSE_LEAFS          = 1 << 0,
  FTT_TRAVERSE_NON_LEAFS      = 1 << 1,
  FTT_TRAVERSE_LEVEL          = 1 << 2,
  FTT_TRAVERSE_BOUNDARY_FACES = 1 << 3,
  FTT_TRAVERSE_ALL            = FTT_TRAVERSE_LEAFS | FTT_TRAVERSE_NON_LEAFS
} FttTraverseFlags;

typedef enum
{
  FTT_PRE_ORDER,
  FTT_POST_ORDER
} FttTraverseType;

typedef enum
{
  FTT_RIGHT = 0,
  FTT_LEFT,
  FTT_TOP,
  FTT_BOTTOM,
#if (!FTT_2D)
  FTT_FRONT,
  FTT_BACK,
#endif /* FTT_3D || FTT_2D3 */
  FTT_NEIGHBORS
} FttDirection;

#define FTT_NEIGHBORS_2D (FTT_BOTTOM + 1)

#if FTT_2D3
# define FTT_CELLS_DIRECTION(d) ((d) < FTT_NEIGHBORS_2D ? FTT_CELLS/2 : FTT_CELLS)
#else  /* 2D && 3D */
# define FTT_CELLS_DIRECTION(d) (FTT_CELLS/2)
#endif /* 2D && 3D */

GTS_C_VAR gchar * ftt_direction_name[FTT_NEIGHBORS]; /* defined in ftt.c */

typedef enum
{
  FTT_X = 0,
  FTT_Y,
#if (!FTT_2D)
  FTT_Z,
#endif /* FTT_3D || FTT_2D3 */
  FTT_DIMENSION,
  FTT_XY,
#if FTT_2D
  FTT_XYZ = FTT_XY
#else  /* FTT_3D || FTT_2D3 */
  FTT_XYZ
#endif /* FTT_3D || FTT_2D3 */
} FttComponent;

typedef enum {
  FTT_FLAG_ID        = 7,
  FTT_FLAG_DESTROYED = 1 << 3,
  FTT_FLAG_LEAF      = 1 << 4,        /* used only for I/O operations */
  FTT_FLAG_TRAVERSED = FTT_FLAG_LEAF, /* used for face traversal */
  FTT_FLAG_USER      =      5         /* user flags start here */
} FttCellFlags;

typedef void      (* FttCellTraverseFunc)            (FttCell * cell,
						      gpointer data);
typedef void      (* FttCellInitFunc)                (FttCell * cell,
						      gpointer data);

struct _FttCellNeighbors {
  /* right, left, top, bottom, front, back */
  FttCell * c[FTT_NEIGHBORS];
};

struct _FttCellChildren {
  FttCell * c[FTT_CELLS];
};

struct _FttCell {
  /*< public >*/
  guint flags;
  gpointer data;

  /*< private >*/
  struct _FttOct * parent, * children;
};

struct _FttRootCell {
  FttCell cell;

  FttCellNeighbors neighbors;
  FttVector pos;
  guint level;
#if FTT_2D3
  gdouble dz;
#endif
};

struct _FttOct {
  guint level;
  FttCell * parent;
  FttCellNeighbors neighbors;
  FttVector pos;
#if FTT_2D3
  gdouble dz;
#endif

  FttCell cell[FTT_CELLS];
};

struct _FttCellFace {
  FttCell * cell, * neighbor;
  FttDirection d;
};

#define  FTT_CELL_ID(c)           ((c)->flags & FTT_FLAG_ID)
#define  FTT_CELL_IS_LEAF(c)      ((c)->children == NULL)
#define  FTT_CELL_IS_ROOT(c)      ((c)->parent == NULL)
#define  FTT_CELL_IS_DESTROYED(c) (((c)->flags & FTT_FLAG_DESTROYED) != 0)

typedef enum {
  FTT_BOUNDARY,
  FTT_FINE_FINE,
  FTT_FINE_COARSE
} FttFaceType;

#define  FTT_FACE_DIRECT(f)       ((f)->d % 2 == 0)
#define  FTT_FACE_REVERSE(dst, src) \
   ((dst)->cell = (src)->neighbor,\
    (dst)->neighbor = (src)->cell,\
    (dst)->d = FTT_OPPOSITE_DIRECTION((src)->d))

GTS_C_VAR
gint                 ftt_opposite_direction[FTT_NEIGHBORS];

#define FTT_OPPOSITE_DIRECTION(d)     (ftt_opposite_direction[d])
#define FTT_ORTHOGONAL_COMPONENT(c)   (((c) + 1) % FTT_DIMENSION)

#ifdef G_DISABLE_ASSERT

#define g_assert_not_implemented()

#else /* !G_DISABLE_ASSERT */

#ifdef __GNUC__

#define g_assert_not_implemented()      G_STMT_START{		\
     g_log (G_LOG_DOMAIN,					\
	    G_LOG_LEVEL_ERROR,					\
	    "file %s: line %d (%s): not implemented (yet)",	\
	    __FILE__,						\
	    __LINE__,						\
	    __PRETTY_FUNCTION__);	}G_STMT_END

#else /* !__GNUC__ */

#define g_assert_not_implemented()	G_STMT_START{	\
     g_log (G_LOG_DOMAIN,				\
	    G_LOG_LEVEL_ERROR,				\
	    "file %s: line %d: not implemented (yet)",	\
	    __FILE__,					\
	    __LINE__);		}G_STMT_END

#endif /* __GNUC__ */

#endif /* !G_DISABLE_ASSERT */

FttCell *            ftt_cell_new                    (FttCellInitFunc init,
						      gpointer data);
#define              ftt_cell_level(c)  ((c)->parent ?\
                                         (c)->parent->level + 1 :\
                                         ((struct _FttRootCell *) c)->level)
#define              ftt_cell_parent(c) ((c)->parent ?\
                                         (c)->parent->parent : NULL)
#ifdef FTT_2D3
# define             ftt_cell_dz(c)     ((c)->parent ?\
                                         (c)->parent->dz :\
                                         ((struct _FttRootCell *) c)->dz)
#else  /* 2D or 3D */
# define             ftt_cell_dz(c)     (1.)
#endif /* 2D or 3D */
guint                ftt_cell_depth                  (const FttCell * root);
G_INLINE_FUNC
gdouble              ftt_level_size                  (guint level);
G_INLINE_FUNC
gdouble              ftt_cell_size                   (const FttCell * cell);
G_INLINE_FUNC
gdouble              ftt_cell_volume                 (const FttCell * cell);
G_INLINE_FUNC
void                 ftt_cell_children           (const FttCell * cell,
						  FttCellChildren * children);
G_INLINE_FUNC
guint                ftt_cell_children_direction (const FttCell * cell,
						  FttDirection d,
						  FttCellChildren * children);
G_INLINE_FUNC
FttCell *            ftt_cell_child_corner       (const FttCell * cell,
						  FttDirection 
						  d[FTT_DIMENSION]);
G_INLINE_FUNC
void                 ftt_cell_neighbors_not_cached (const FttCell * cell,
						 FttCellNeighbors * neighbors);
G_INLINE_FUNC
FttCell *            ftt_cell_neighbor_not_cached (const FttCell * cell,
						   FttDirection d);
G_INLINE_FUNC
void                 ftt_cell_neighbors         (const FttCell * cell,
						 FttCellNeighbors * neighbors);
G_INLINE_FUNC
FttCell *            ftt_cell_neighbor          (const FttCell * cell,
						 FttDirection d);
G_INLINE_FUNC
FttCellFace          ftt_cell_face              (FttCell * cell,
						 FttDirection d);
G_INLINE_FUNC
FttFaceType          ftt_face_type              (const FttCellFace * face);

G_INLINE_FUNC
gboolean             ftt_cell_neighbor_is_brother (FttCell * cell, 
						   FttDirection d);

#ifdef G_CAN_INLINE
G_INLINE_FUNC
gdouble ftt_level_size (guint level)
{
  gdouble size = 1.;

  while (level) {
    size /= 2.;
    level--;
  }

  return size;
}

G_INLINE_FUNC
gdouble ftt_cell_size (const FttCell * cell)
{
  g_return_val_if_fail (cell != NULL, 0.);

  return ftt_level_size (ftt_cell_level (cell));
}

G_INLINE_FUNC
gdouble ftt_cell_volume (const FttCell * cell)
{
  gdouble size;

  g_return_val_if_fail (cell != NULL, 0.);

  size = ftt_level_size (ftt_cell_level (cell));
#if (FTT_2D || FTT_2D3)
  return size*size;
#else  /* FTT_3D */
  return size*size*size;
#endif /* FTT_3D */
}

G_INLINE_FUNC
void ftt_cell_children (const FttCell * cell,
			FttCellChildren * children)
{
  struct _FttOct * oct;
  guint i;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (!FTT_CELL_IS_LEAF (cell));
  g_return_if_fail (children != NULL);

  oct = cell->children;
  for (i = 0; i < FTT_CELLS; i++)
    children->c[i] = FTT_CELL_IS_DESTROYED (&(oct->cell[i])) ? NULL : &(oct->cell[i]);
}

G_INLINE_FUNC
guint ftt_cell_children_direction (const FttCell * cell,
				  FttDirection d,
				  FttCellChildren * children)
{
  struct _FttOct * oct;
  guint i;
#if (FTT_2D || FTT_2D3)
  static gint index[FTT_NEIGHBORS_2D][FTT_CELLS/2] =
  {{1, 3},
   {0, 2},
   {0, 1},
   {2, 3}};
#else  /* FTT_3D */
  static gint index[FTT_NEIGHBORS][FTT_CELLS/2] =
  {{1, 3, 5, 7},
   {0, 2, 4, 6},
   {0, 1, 4, 5},
   {2, 3, 6, 7},
   {0, 1, 2, 3},
   {4, 5, 6, 7}};
#endif /* FTT_3D */

  g_return_val_if_fail (cell != NULL, 0);
  g_return_val_if_fail (!FTT_CELL_IS_LEAF (cell), 0);
  g_return_val_if_fail (d < FTT_NEIGHBORS, 0);
  g_return_val_if_fail (children != NULL, 0);

  oct = cell->children;

#if FTT_2D3
  if (d >= FTT_NEIGHBORS_2D) {
    for (i = 0; i < FTT_CELLS; i++)
      children->c[i] = FTT_CELL_IS_DESTROYED (&(oct->cell[i])) ? NULL : &(oct->cell[i]);
    return FTT_CELLS;
  }
#endif /* 2D3 */

  for (i = 0; i < FTT_CELLS/2; i++)
    children->c[i] = FTT_CELL_IS_DESTROYED (&(oct->cell[index[d][i]])) ? 
      NULL : &(oct->cell[index[d][i]]);
  return FTT_CELLS/2;
}

G_INLINE_FUNC
FttCell * ftt_cell_child_corner (const FttCell * cell,
				 FttDirection d[FTT_DIMENSION])
{
#if (FTT_2D || FTT_2D3)
  static gint index[FTT_NEIGHBORS_2D][FTT_NEIGHBORS_2D] = {
    {-1,-1,1,3},
    {-1,-1,0,2},
    {1,0,-1,-1},
    {3,2,-1,-1}
  };
  gint i;

  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (!FTT_CELL_IS_LEAF (cell), NULL);

  g_return_val_if_fail (d[0] < FTT_NEIGHBORS, NULL);
  g_return_val_if_fail (d[1] < FTT_NEIGHBORS, NULL);

#  if FTT_2D3
  if (d[0] >= FTT_NEIGHBORS_2D)
    i = index[d[1]][d[2]];
  else if (d[1] >= FTT_NEIGHBORS_2D)
    i = index[d[0]][d[2]];
  else
#  endif
    i = index[d[0]][d[1]];
#else  /* FTT_3D */
  static gint index[FTT_NEIGHBORS][FTT_NEIGHBORS][FTT_NEIGHBORS] = {
    {{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},
     {-1,-1,-1,-1,1,5},{-1,-1,-1,-1,3,7},
     {-1,-1,1,3,-1,-1},{-1,-1,5,7,-1,-1}},
    {{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},
     {-1,-1,-1,-1,0,4},{-1,-1,-1,-1,2,6},
     {-1,-1,0,2,-1,-1},{-1,-1,4,6,-1,-1}},
    {{-1,-1,-1,-1,1,5},{-1,-1,-1,-1,0,4},
     {-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},
     {1,0,-1,-1,-1,-1},{5,4,-1,-1,-1,-1}},
    {{-1,-1,-1,-1,3,7},{-1,-1,-1,-1,2,6},
     {-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},
     {3,2,-1,-1,-1,-1},{7,6,-1,-1,-1,-1}},
    {{-1,-1,1,3,-1,-1},{-1,-1,0,2,-1,-1},
     {1,0,-1,-1,-1,-1},{3,2,-1,-1,-1,-1},
     {-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1}},
    {{-1,-1,5,7,-1,-1},{-1,-1,4,6,-1,-1},
     {5,4,-1,-1,-1,-1},{7,6,-1,-1,-1,-1},
     {-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1}},
  };
  gint i;

  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (!FTT_CELL_IS_LEAF (cell), NULL);
  g_return_val_if_fail (d[0] < FTT_NEIGHBORS, NULL);
  g_return_val_if_fail (d[1] < FTT_NEIGHBORS, NULL);
  g_return_val_if_fail (d[2] < FTT_NEIGHBORS, NULL);

  i = index[d[0]][d[1]][d[2]];
#endif /* FTT_3D */

  g_return_val_if_fail (i >= 0, NULL);

  return FTT_CELL_IS_DESTROYED (&(cell->children->cell[i])) ? NULL:
    &(cell->children->cell[i]);
}

G_INLINE_FUNC
void ftt_cell_neighbors_not_cached (const FttCell * cell,
				    FttCellNeighbors * neighbors)
{
  static gint neighbor_index[FTT_NEIGHBORS][FTT_CELLS]
#if FTT_2D
    = {{1,-1,3,-3},
       {-2,0,-4,2},
       {-3,-4,0,1},
       {2,3,-1,-2}};
#elif FTT_2D3
    = {{1,-1,3,-3},
       {-2,0,-4,2},
       {-3,-4,0,1},
       {2,3,-1,-2},
       {-1,-2,-3,-4},
       {-1,-2,-3,-4}};
#else  /* FTT_3D */
    = {{1,-1,3,-3,5,-5,7,-7},
       {-2,0,-4,2,-6,4,-8,6},
       {-3,-4,0,1,-7,-8,4,5},
       {2,3,-1,-2,6,7,-5,-6},
       {-5,-6,-7,-8,0,1,2,3},
       {4,5,6,7,-1,-2,-3,-4}};
#endif /* FTT_3D */
  guint n, d;
  struct _FttOct * parent;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (neighbors != NULL);

  if (FTT_CELL_IS_ROOT (cell)) {
    memcpy (neighbors, &((struct _FttRootCell *) cell)->neighbors,
	    sizeof (FttCellNeighbors));
    return;
  }

  parent = cell->parent;
  n = FTT_CELL_ID (cell);
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    gint nn = neighbor_index[d][n];
    FttCell * c;

    if (nn >= 0) /* neighbor belongs to same Oct */
      c = &(parent->cell[nn]);
    else {       /* neighbor belongs to neighboring Cell or Oct */
      c = parent->neighbors.c[d];
      if (c != NULL && c->children != NULL)
	c = &(c->children->cell[- nn - 1]);
    }
    if (c == NULL || FTT_CELL_IS_DESTROYED (c))
      neighbors->c[d] = NULL;
    else
      neighbors->c[d] = c;
  }
}

G_INLINE_FUNC
FttCell * ftt_cell_neighbor_not_cached (const FttCell * cell,
					FttDirection d)
{
  static gint neighbor_index[FTT_NEIGHBORS][FTT_CELLS]
#if FTT_2D
    = {{1,-1,3,-3},
       {-2,0,-4,2},
       {-3,-4,0,1},
       {2,3,-1,-2}};
#elif FTT_2D3
    = {{1,-1,3,-3},
       {-2,0,-4,2},
       {-3,-4,0,1},
       {2,3,-1,-2},
       {-1,-2,-3,-4},
       {-1,-2,-3,-4}};
#else  /* FTT_3D */
    = {{1,-1,3,-3,5,-5,7,-7},
       {-2,0,-4,2,-6,4,-8,6},
       {-3,-4,0,1,-7,-8,4,5},
       {2,3,-1,-2,6,7,-5,-6},
       {-5,-6,-7,-8,0,1,2,3},
       {4,5,6,7,-1,-2,-3,-4}};
#endif /* FTT_3D */
  gint n;
  FttCell * c;

  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (d < FTT_NEIGHBORS, NULL);

  if (FTT_CELL_IS_ROOT (cell))
    return ((struct _FttRootCell *) cell)->neighbors.c[d];

  n = neighbor_index[d][FTT_CELL_ID (cell)];
  if (n >= 0) /* neighbor belongs to same Oct */
    c = &(cell->parent->cell[n]);
  else {      /* neighbor belongs to neighboring Cell or Oct */
    c = cell->parent->neighbors.c[d];
    if (c != NULL && c->children != NULL)
      c = &(c->children->cell[- n - 1]);
  }
  if (c == NULL || FTT_CELL_IS_DESTROYED (c))
    return NULL;
  else
    return c;
}

G_INLINE_FUNC
void ftt_cell_neighbors (const FttCell * cell,
			 FttCellNeighbors * neighbors)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (neighbors != NULL);

  if (!FTT_CELL_IS_LEAF (cell) && neighbors != &cell->children->neighbors) {
    memcpy (neighbors, &cell->children->neighbors, 
	    sizeof (FttCellNeighbors));
    return;
  }

  ftt_cell_neighbors_not_cached (cell, neighbors);
}

G_INLINE_FUNC
FttCell * ftt_cell_neighbor (const FttCell * cell,
			     FttDirection d)
{
  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (d < FTT_NEIGHBORS, NULL);

  if (!FTT_CELL_IS_LEAF (cell))
    return cell->children->neighbors.c[d];

  return ftt_cell_neighbor_not_cached (cell, d);
}

G_INLINE_FUNC
FttCellFace ftt_cell_face (FttCell * cell,
			   FttDirection d)
{
  FttCellFace f;

  g_assert (cell != NULL);

  f.cell = cell;
  f.neighbor = ftt_cell_neighbor (cell, d);
  f.d = d;

  return f;
}

G_INLINE_FUNC
FttFaceType ftt_face_type (const FttCellFace * face)
{
  g_return_val_if_fail (face != NULL, 0);

  if (face->neighbor == NULL)
    return FTT_BOUNDARY;
  if (ftt_cell_level (face->cell) > ftt_cell_level (face->neighbor))
    return FTT_FINE_COARSE;
  g_assert (ftt_cell_level (face->cell) == ftt_cell_level (face->neighbor));
  return FTT_FINE_FINE;
}

G_INLINE_FUNC
gboolean ftt_cell_neighbor_is_brother (FttCell * cell, 
				       FttDirection d)
{
  static gboolean b[FTT_CELLS][FTT_NEIGHBORS] = {
#if FTT_2D
    {1,0,0,1}, {0,1,0,1}, {1,0,1,0}, {0,1,1,0}
#elif FTT_2D3
    {1,0,0,1,0,0}, {0,1,0,1,0,0}, {1,0,1,0,0,0}, {0,1,1,0,0,0}
#else  /* 3D */
    {1,0,0,1,0,1}, {0,1,0,1,0,1}, {1,0,1,0,0,1}, {0,1,1,0,0,1},
    {1,0,0,1,1,0}, {0,1,0,1,1,0}, {1,0,1,0,1,0}, {0,1,1,0,1,0}
#endif /* 3D */
  };

  g_return_val_if_fail (cell != NULL, FALSE);
  
  if (FTT_CELL_IS_ROOT (cell))
    return FALSE;
  return b[FTT_CELL_ID (cell)][d];
}
#endif	/* G_CAN_INLINE */

void                 ftt_cell_set_neighbor           (FttCell * root,
						      FttCell * neighbor,
						      FttDirection d,
						      FttCellInitFunc init,
						      gpointer init_data);
void                 ftt_cell_set_neighbor_match     (FttCell * root,
						      FttCell * neighbor,
						      FttDirection d,
						      FttCellInitFunc init,
						      gpointer init_data);
void                 ftt_cell_relative_pos           (const FttCell * cell,
						      FttVector * pos);
void                 ftt_cell_pos                    (const FttCell * cell,
						      FttVector * pos);
void                 ftt_corner_relative_pos         (const FttCell * cell,
						      FttDirection d[FTT_DIMENSION],
						      FttVector * pos);
void                 ftt_corner_pos                  (const FttCell * cell,
						      FttDirection d[FTT_DIMENSION],
						      FttVector * pos);
void                 ftt_face_pos                    (const FttCellFace * face,
						      FttVector * pos);
void                 ftt_cell_set_pos                (FttCell * root,
						      const FttVector * pos);
void                 ftt_cell_set_level              (FttCell * root,
						      guint level);
void                 ftt_cell_draw                   (const FttCell * cell,
						      FILE * fp);
void                 ftt_face_draw                   (const FttCellFace * face,
						      FILE * fp);
gboolean             ftt_cell_check                  (const FttCell * cell);
typedef gboolean  (* FttCellRefineFunc)              (FttCell * cell,
						      gpointer data);
void                 ftt_cell_refine                 (FttCell * root,
						      FttCellRefineFunc refine,
						      gpointer refine_data,
						      FttCellInitFunc init,
						      gpointer init_data);
void                 ftt_cell_refine_single          (FttCell * cell,
						      FttCellInitFunc init,
						      gpointer init_data);
gboolean             ftt_refine_corner               (const FttCell * cell);
void                 ftt_cell_traverse               (FttCell * root,
						      FttTraverseType order,
						      FttTraverseFlags flags,
						      gint max_depth,
						      FttCellTraverseFunc func,
						      gpointer data);
void                 ftt_cell_traverse_condition     (FttCell * root,
						      FttTraverseType order,
						      FttTraverseFlags flags,
						      gint max_depth,
						      FttCellTraverseFunc func,
						      gpointer data,
						      gboolean (* condition) (FttCell *, 
									      gpointer),
						      gpointer cdata);
void                 ftt_cell_traverse_box           (FttCell * root,
						      GtsBBox * box,
						      FttTraverseType order,
						      FttTraverseFlags flags,
						      gint max_depth,
						      FttCellTraverseFunc func,
						      gpointer data);
void                 ftt_cell_traverse_boundary      (FttCell * root,
						      FttDirection d,
						      FttTraverseType order,
						      FttTraverseFlags flags,
						      gint max_depth,
						      FttCellTraverseFunc func,
						      gpointer data);
typedef void      (* FttFaceTraverseFunc)            (FttCellFace * face, 
						      gpointer data);
void                 ftt_face_traverse               (FttCell * root,
						      FttComponent c,
						      FttTraverseType order,
						      FttTraverseFlags flags,
						      gint max_depth,
						      FttFaceTraverseFunc func,
						      gpointer data);
void                 ftt_face_traverse_boundary      (FttCell * root,
						      FttDirection d,
						      FttTraverseType order,
						      FttTraverseFlags flags,
						      gint max_depth,
						      FttFaceTraverseFunc func,
						      gpointer data);
FttCell *            ftt_cell_locate                 (FttCell * root,
						      FttVector target,
						      gint max_depth);
gdouble              ftt_cell_point_distance2_min    (FttCell * cell, 
						      GtsPoint * p);
void                 ftt_cell_point_distance2_internal (FttCell * root,
							GtsPoint * p,
							gdouble d,
							gdouble (* distance2) (FttCell *, 
									       GtsPoint *, 
									       gpointer),
							gpointer data,
							FttCell ** closest,
							gdouble * dmin);
gdouble              ftt_cell_point_distance2        (FttCell * root,
						      GtsPoint * p,
						      gdouble (* distance2) (FttCell *, 
									     GtsPoint *, 
									     gpointer),
						      gpointer data,
						      FttCell ** closest);
void                 ftt_cell_bbox                   (const FttCell * cell, 
						      GtsBBox * bb);
typedef void      (* FttCellCopyFunc)                (const FttCell * from,
						      FttCell * to,
						      gpointer data);
FttCell *            ftt_cell_copy                   (const FttCell * root,
						      FttCellCopyFunc copy,
						      gpointer data);
typedef void      (* FttCellWriteFunc)               (const FttCell * cell,
						      FILE * fp,
						      gpointer data);
void                 ftt_cell_write                  (const FttCell * root,
						      gint max_depth,
						      FILE * fp,
						      FttCellWriteFunc write,
						      gpointer data);
void                 ftt_cell_write_binary           (const FttCell * root,
						      gint max_depth,
						      FILE * fp,
						      FttCellWriteFunc write,
						      gpointer data);
typedef void      (* FttCellReadFunc)                (FttCell * cell,
						      GtsFile * fp,
						      gpointer data);
FttCell *            ftt_cell_read                   (GtsFile * fp,
						      FttCellReadFunc read,
						      gpointer data);
FttCell *            ftt_cell_read_binary            (GtsFile * fp,
						      FttCellReadFunc read,
						      gpointer data);
typedef void      (* FttCellCleanupFunc)             (FttCell * cell,
						      gpointer data);
void                 ftt_cell_destroy           (FttCell * cell,
						 FttCellCleanupFunc cleanup,
						 gpointer data);
void                 ftt_cell_destroy_root      (FttCell * root,
						 FttCellChildren * children,
						 FttCellCleanupFunc cleanup,
						 gpointer data);
void                 ftt_cell_flatten           (FttCell * root, 
						 FttDirection d, 
						 FttCellCleanupFunc cleanup,
						 gpointer data);
typedef gboolean  (* FttCellCoarsenFunc)        (FttCell * cell,
						 gpointer data);
gboolean             ftt_cell_coarsen           (FttCell * root,
						 FttCellCoarsenFunc coarsen,
						 gpointer coarsen_data,
						 FttCellCleanupFunc cleanup,
						 gpointer cleanup_data);
FttDirection         ftt_direction_from_name    (const gchar * name);

struct _FttCellTraverse {
  FttCell ** cells;
  FttCell ** current;
};

typedef struct _FttCellTraverse FttCellTraverse;

FttCellTraverse *    ftt_cell_traverse_new      (FttCell * root,
						 FttTraverseType order,
						 FttTraverseFlags flags,
						 gint max_depth);
void                 ftt_cell_traverse_rewind   (FttCellTraverse * t);
void                 ftt_cell_traverse_destroy  (FttCellTraverse * t);
G_INLINE_FUNC
FttCell *            ftt_cell_traverse_next     (FttCellTraverse * t);

#ifdef G_CAN_INLINE
G_INLINE_FUNC
FttCell * ftt_cell_traverse_next (FttCellTraverse * t)
{
  g_return_val_if_fail (t != NULL, NULL);

  return *(t->current++);
}
#endif /* G_CAN_INLINE */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FTT_H__ */
