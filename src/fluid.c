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

#include "fluid.h"
#include "variable.h"
#include "domain.h"
#include "solid.h"

/* Permanent variables */
GfsVariable * gfs_div, * gfs_dp, * gfs_res;
GfsVariable * gfs_gx, * gfs_gy;
#if (!FTT_2D)
GfsVariable * gfs_gz;
#endif /* FTT_3D */
GfsVariable * gfs_centered_variables;
GfsVariable * gfs_p;

/* Derived variables */
GfsVariable * gfs_derived_first, * gfs_derived_last;

typedef struct _Gradient Gradient;

/* grad(p) = -a*p(cell) + b*p(neighbor) + c */
struct _Gradient {
  gdouble a, b, c;
};

static gdouble average_neighbor_value (const FttCellFace * face,
				       guint v,
				       gdouble * x)
{
  /* check for corner refinement violation (topology.fig) */
  g_assert (ftt_cell_level (face->neighbor) == ftt_cell_level (face->cell));
  
  if (FTT_CELL_IS_LEAF (face->neighbor))
    return GFS_VARIABLE (face->neighbor, v);
  else {
    FttCellChildren children;
    gdouble av = 0., a = 0.;
    guint i, n;
    
    n = ftt_cell_children_direction (face->neighbor,
				     FTT_OPPOSITE_DIRECTION (face->d),
				     &children);
    for (i = 0; i < n; i++) 
      if (children.c[i]) {
	a += 1.;
	av += GFS_VARIABLE (children.c[i], v);
      }
    if (a > 0.) {
      *x = 3./4.;
      return av/a;
    }
    else
      return GFS_VARIABLE (face->cell, v);
  }
}

#if (FTT_2D || FTT_2D3)

/* v = a*v(cell) + b 
 * 
 * Second order 1D interpolation.
 */
static GfsGradient interpolate_1D2 (FttCell * cell,
				    FttDirection d,
				    gdouble x,
				    guint v)
{
  GfsGradient p;
  FttCellFace f1, f2;
  gdouble p1 = 0., p2 = 0.;
  gdouble x1 = 1., x2 = 1.;
  gdouble a1, a2;

  g_return_val_if_fail (cell != NULL, p);
  g_return_val_if_fail (!GFS_IS_MIXED (cell), p);
  
  f1 = ftt_cell_face (cell, FTT_OPPOSITE_DIRECTION (d));
  if (f1.neighbor)
    p1 = average_neighbor_value (&f1, v, &x1);
  f2 = ftt_cell_face (cell, d);
  if (f2.neighbor)
    p2 = average_neighbor_value (&f2, v, &x2);

  a1 = x*(x - x2)/(x1*(x1 + x2));
  a2 = x*(x + x1)/(x2*(x1 + x2));

  p.a = 1. - a1 - a2;
  p.b = 0.;
  if (f1.neighbor)
    p.b += a1*p1;
  else
    p.a += a1;
  if (f2.neighbor)
    p.b += a2*p2;
  else
    p.a += a2;

  return p;
}

/* v = a*v(cell) + b 
 * 
 * First order 1D interpolation.
 */
static GfsGradient interpolate_1D1 (FttCell * cell,
				    FttDirection d,
				    gdouble x,
				    guint v)
{
  GfsGradient p;
  FttCellFace f;

  g_return_val_if_fail (cell != NULL, p);

  f = ftt_cell_face (cell, d);
  if (f.neighbor) {
    gdouble p2;
    gdouble x2 = 1.;

    p2 = average_neighbor_value (&f, v, &x2);
    p.a = 1. - x/x2;
    p.b = p2*x/x2;
  }
  else {
    p.a = 1.;
    p.b = 0.;
  }

  return p;
}

#else /* not FTT_2D */

/* v = a*v(cell) + b 
 * 
 * First order 2D interpolation.
 */
static GfsGradient interpolate_2D1 (FttCell * cell,
				    FttDirection d1, FttDirection d2,
				    gdouble x, gdouble y,
				    guint v)
{
  GfsGradient p;
  gdouble y1 = 1.;
  gdouble x2 = 1.;
  gdouble p1 = 0., p2 = 0.;
  gdouble a1, a2;
  FttCellFace f1, f2;

  g_return_val_if_fail (cell != NULL, p);
  
  f1 = ftt_cell_face (cell, d1);
  if (f1.neighbor)
    p1 = average_neighbor_value (&f1, v, &y1);
  f2 = ftt_cell_face (cell, d2);
  if (f2.neighbor)
    p2 = average_neighbor_value (&f2, v, &x2);

  a1 = y/y1;
  a2 = x/x2;

  p.a = 1. - a1 - a2;
  p.b = 0.;
  if (f1.neighbor)
    p.b += a1*p1;
  else
    p.a += a1;
  if (f2.neighbor)
    p.b += a2*p2;
  else
    p.a += a2;
  
  return p;
}

#endif /* not FTT_2D */

#if (FTT_2D || FTT_2D3)
static gint perpendicular[FTT_NEIGHBORS][FTT_CELLS] = 
  {{-1,  2, -1,  3},
   { 2, -1,  3, -1},
   { 1,  0, -1, -1},
   {-1, -1,  1,  0}};
#else  /* FTT_3D */
static gint perpendicular[FTT_NEIGHBORS][FTT_CELLS][2] = 
  {{{-1,-1},{2,4},{-1,-1},{3,4},{-1,-1},{2,5},{-1,-1},{3,5}},
   {{2,4},{-1,-1},{3,4},{-1,-1},{2,5},{-1,-1},{3,5},{-1,-1}},
   {{4,1},{4,0},{-1,-1},{-1,-1},{5,1},{5,0},{-1,-1},{-1,-1}},
   {{-1,-1},{-1,-1},{4,1},{4,0},{-1,-1},{-1,-1},{5,1},{5,0}},
   {{1,2},{0,2},{1,3},{0,3},{-1,-1},{-1,-1},{-1,-1},{-1,-1}},
   {{-1,-1},{-1,-1},{-1,-1},{-1,-1},{1,2},{0,2},{1,3},{0,3}}};
#endif /* FTT_3D */

static Gradient gradient_fine_coarse (const FttCellFace * face,
				      guint v,
				      gint max_level)
{
  Gradient g;
  guint level;
  FttCell * n;
  GfsGradient p;
  GfsSolidVector * s;
#if (FTT_2D || FTT_2D3)
  gint dp;
#else  /* FTT_3D */
  gint * dp;
#endif /* FTT_3D */

  g_return_val_if_fail (face != NULL, g);
  g_return_val_if_fail (ftt_face_type (face) == FTT_FINE_COARSE, g);

  level = ftt_cell_level (face->cell);
  s = GFS_STATE (face->cell)->solid;
  n = (s != NULL && s->s[FTT_OPPOSITE_DIRECTION (face->d)] == 0.) ? NULL :
    ftt_cell_neighbor (face->cell, FTT_OPPOSITE_DIRECTION (face->d));
  g_assert (n == NULL || ftt_cell_level (n) == level);

  dp = perpendicular[face->d][FTT_CELL_ID (face->cell)];
#if (FTT_2D || FTT_2D3)
  g_assert (dp >= 0);
  p = interpolate_1D1 (face->neighbor, dp, 1./4., v);
#else  /* FTT_3D */
  g_assert (dp[0] >= 0 && dp[1] >= 0);
  p = interpolate_2D1 (face->neighbor, dp[0], dp[1], 1./4., 1./4., v);
#endif /* FTT_3D */

  if (n == NULL) {
    g.a = 2./3.;
    g.b = 2.*p.a/3.;
    g.c = 2.*p.b/3.;
  }
  else if (level == max_level || FTT_CELL_IS_LEAF (n)) {
    g.a = 1./3.;
    g.b = 8.*p.a/15.;
    g.c = 8.*p.b/15. - GFS_VARIABLE (n, v)/5.;
  }
  else {
    FttCellChildren children;
    guint i, j;
    gdouble sa = 0.;

    j = ftt_cell_children_direction (n, face->d, &children);
    g.a = 2./9.;
    g.b = 14.*p.a/27.;
    g.c = 0.;
    for (i = 0; i < j; i++) 
      if (children.c[i]) {
	g.c += GFS_VARIABLE (children.c[i], v);
	sa += 1.;
      }
    g.c = (14.*p.b - 8.*g.c/sa)/27.;
  }

  return g;
}

#define REFINE_CORNER(cell) {if (cell && FTT_CELL_IS_LEAF (cell) && \
                              ftt_cell_level (cell) < level - 1) \
	                    ftt_cell_refine_single (cell, init, data);}

void ftt_cell_refine_corners (FttCell * cell,
			      FttCellInitFunc init,
			      gpointer data)
{
  FttDirection d;
  FttCellNeighbors neighbor;
  guint level;

  g_return_if_fail (cell != NULL);

  level = ftt_cell_level (cell);
  ftt_cell_neighbors (cell, &neighbor);
#if FTT_2D3
  for (d = 0; d < FTT_NEIGHBORS_2D; d++)
#else /* 2D && 3D */
  for (d = 0; d < FTT_NEIGHBORS; d++)
#endif  /* 2D && 3D */
    if (neighbor.c[d] && ftt_cell_level (neighbor.c[d]) < level) {
      if (GFS_CELL_IS_BOUNDARY (neighbor.c[d]))
	ftt_cell_refine_single (neighbor.c[d], init, data);
      else {
	FttCell * n;
#if (FTT_2D || FTT_2D3)
	gint dp;
#else  /* FTT_3D */
	gint * dp;
#endif /* FTT_3D */
	
	dp = perpendicular[d][FTT_CELL_ID (cell)];
#if (FTT_2D || FTT_2D3)
	g_assert (dp >= 0);
	n = ftt_cell_neighbor (neighbor.c[d], dp);
	REFINE_CORNER (n)
#else  /* FTT_3D */
	g_assert (dp[0] >= 0 && dp[1] >= 0);
	n = ftt_cell_neighbor (neighbor.c[d], dp[0]);
	REFINE_CORNER (n)
	n = ftt_cell_neighbor (neighbor.c[d], dp[1]);
	REFINE_CORNER (n)
#endif /* FTT_3D */
      }
    }
}

static gdouble neighbor_value (const FttCellFace * face,
			       guint v,
			       gdouble * x)
{
  GfsGradient vc;
#if (FTT_2D || FTT_2D3)
  gint dp;
#else  /* FTT_3D */
  gint * dp;
#endif /* FTT_3D */

  if (ftt_cell_level (face->neighbor) == ftt_cell_level (face->cell))
    /* neighbor at same level */
    return average_neighbor_value (face, v, x);
  else {
    /* neighbor at coarser level */
    dp = perpendicular[face->d][FTT_CELL_ID (face->cell)];
#if (FTT_2D || FTT_2D3)
    g_assert (dp >= 0);
    vc = interpolate_1D1 (face->neighbor, dp, 1./4., v);
#else  /* FTT_3D */
    g_assert (dp[0] >= 0 && dp[1] >= 0);
    vc = interpolate_2D1 (face->neighbor, dp[0], dp[1], 1./4., 1./4., v);
#endif /* FTT_3D */
    *x = 3./2.;
    return vc.a*GFS_VARIABLE (face->neighbor, v) + vc.b;
  }
}

/**
 * gfs_center_gradient:
 * @cell: a #FttCell.
 * @c: a component.
 * @v: a #GfsVariable index.
 *
 * The gradient is normalized by the size of the cell.
 *
 * Returns: the value of the @c component of the gradient of variable @v
 * at the center of the cell.  
 */
gdouble gfs_center_gradient (FttCell * cell,
			     FttComponent c,
			     guint v)
{
  FttDirection d = 2*c;
  FttCellFace f1;
  gdouble v0;

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (c < FTT_DIMENSION, 0.);

  f1 = ftt_cell_face (cell, FTT_OPPOSITE_DIRECTION (d));
  if (f1.neighbor == cell) /* periodic */
    return 0.;
  v0 = GFS_VARIABLE (cell, v);
  if (f1.neighbor) {
    FttCellFace f2 = ftt_cell_face (cell, d);
    gdouble x1 = 1., v1;
    
    v1 = neighbor_value (&f1, v, &x1);
    if (f2.neighbor) {
      /* two neighbors: second-order differencing (parabola) */
      gdouble x2 = 1., v2;
      
      v2 = neighbor_value (&f2, v, &x2);
      return (x1*x1*(v2 - v0) + x2*x2*(v0 - v1))/(x1*x2*(x2 + x1));
    }
    else
      /* one neighbor: first-order differencing */
      return (v0 - v1)/x1;
  }
  else {
    FttCellFace f2 = ftt_cell_face (cell, d);

    if (f2.neighbor) {
      gdouble x2 = 1.;
      
      /* one neighbor: first-order differencing */
      return (neighbor_value (&f2, v, &x2) - v0)/x2;
    }
  }
  /* no neighbors */
  return 0.;
}

/**
 * gfs_center_van_leer_gradient:
 * @cell: a #FttCell.
 * @c: a component.
 * @v: a #GfsVariable index.
 *
 * The gradient is normalized by the size of the cell and is limited
 * using van Leer's limiter.
 *
 * Returns: the value of the @c component of the gradient of variable @v
 * at the center of the cell.  
 */
gdouble gfs_center_van_leer_gradient (FttCell * cell,
				      FttComponent c,
				      guint v)
{
  FttDirection d = 2*c;
  FttCellFace f1;
  
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (c < FTT_DIMENSION, 0.);

  f1 = ftt_cell_face (cell, FTT_OPPOSITE_DIRECTION (d));
  if (f1.neighbor == cell) /* periodic */
    return 0.;
  if (f1.neighbor) {
    FttCellFace f2 = ftt_cell_face (cell, d);
    
    if (f2.neighbor) {
      /* two neighbors: second-order differencing (parabola)
	 + van Leer limiter */
      gdouble x1 = 1., x2 = 1., v0, v1, v2;
      gdouble s0, s1, s2;
      
      v0 = GFS_VARIABLE (cell, v);
      v1 = neighbor_value (&f1, v, &x1);
      v2 = neighbor_value (&f2, v, &x2);      
      
      s1 = 2.*(v0 - v1);
      s2 = 2.*(v2 - v0);

      if (s1*s2 <= 0.)
	return 0.;
      s0 = (x1*x1*(v2 - v0) + x2*x2*(v0 - v1))/(x1*x2*(x2 + x1));
      if (ABS (s2) < ABS (s1))
	s1 = s2;
      if (ABS (s0) < ABS (s1))
	return s0;
      return s1;
    }
  }
  /* only one or no neighbors */
  return 0.;
}
					 
/**
 * gfs_face_gradient:
 * @face: a #FttCellFace.
 * @g: the #GfsGradient.
 * @v: a #GfsVariable index.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 *
 * Set the value of @g as the gradient of variable @v on the
 * @face. The value returned is second order accurate in space and
 * conservative, in the sense that values at a coarse/fine cell
 * boundary are consistent.  
 */
void gfs_face_gradient (const FttCellFace * face,
			GfsGradient * g,
			guint v,
			gint max_level)
{
  guint level;

  g_return_if_fail (face != NULL);

  g->a = g->b = 0.;
  if (face->neighbor == NULL || GFS_FACE_FRACTION (face) == 0.)
    return;

  level = ftt_cell_level (face->cell);
  if (ftt_cell_level (face->neighbor) < level) {
    /* neighbor is at a shallower level */
    Gradient gcf;

    gcf = gradient_fine_coarse (face, v, max_level);
    g->a = gcf.a;
    g->b = gcf.b*GFS_VARIABLE (face->neighbor, v) + gcf.c;
  }
  else {
    if (level == max_level || FTT_CELL_IS_LEAF (face->neighbor)) {
      /* neighbor is at the same level */
      g->a = 1.;
      g->b = GFS_VARIABLE (face->neighbor, v);
    }
    else {
      /* neighbor is at a deeper level */
      FttCellChildren children;
      FttCellFace f;
      guint i, n;
      gdouble s;
      
      f.d = FTT_OPPOSITE_DIRECTION (face->d);
      n = ftt_cell_children_direction (face->neighbor, f.d, &children);
      f.neighbor = face->cell;
      for (i = 0; i < n; i++) 
	if ((f.cell = children.c[i])) {
	  Gradient gcf;	  
	  gcf = gradient_fine_coarse (&f, v, max_level);
	  s = GFS_FACE_FRACTION (&f);
	  g->a += s*gcf.b;
	  g->b += s*(gcf.a*GFS_VARIABLE (f.cell, v) - gcf.c);
	}
      s = GFS_FACE_FRACTION (face);
#if (FTT_2D || FTT_2D3)
      g->a /= s;
      g->b /= s;
#else  /* 3D */
      g->a /= 2.*s;
      g->b /= 2.*s;
#endif /* 3D */
    }
  }
}

/**
 * gfs_face_weighted_gradient:
 * @face: a #FttCellFace.
 * @g: the #GfsGradient.
 * @v: a #GfsVariable index.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 *
 * Set the value of @g as the gradient of variable @v on the @face
 * weighted by the value of the @v field of the face state vector of the
 * corresponding cell. The value returned is second order accurate in
 * space and conservative, in the sense that values at a coarse/fine
 * cell boundary are consistent.  
 */
void gfs_face_weighted_gradient (const FttCellFace * face,
				 GfsGradient * g,
				 guint v,
				 gint max_level)
{
  guint level;

  g_return_if_fail (face != NULL);

  g->a = g->b = 0.;
  if (face->neighbor == NULL)
    return;

  level = ftt_cell_level (face->cell);
  if (ftt_cell_level (face->neighbor) < level) {
    /* neighbor is at a shallower level */
    Gradient gcf;
    gdouble w = GFS_STATE (face->cell)->f[face->d].v;

    gcf = gradient_fine_coarse (face, v, max_level);
    g->a = w*gcf.a;
    g->b = w*(gcf.b*GFS_VARIABLE (face->neighbor, v) + gcf.c);
  }
  else {
    if (level == max_level || FTT_CELL_IS_LEAF (face->neighbor)) {
      /* neighbor is at the same level */
      gdouble w = GFS_STATE (face->cell)->f[face->d].v;

      g->a = w;
      g->b = w*GFS_VARIABLE (face->neighbor, v);
    }
    else {
      /* neighbor is at a deeper level */
      FttCellChildren children;
      FttCellFace f;
      guint i, n;
      
      f.d = FTT_OPPOSITE_DIRECTION (face->d);
      n = ftt_cell_children_direction (face->neighbor, f.d, &children);
      f.neighbor = face->cell;
      for (i = 0; i < n; i++) 
	if ((f.cell = children.c[i])) {
	  Gradient gcf;
	  gdouble w = GFS_STATE (f.cell)->f[f.d].v;
	
	  gcf = gradient_fine_coarse (&f, v, max_level);
	  g->a += w*gcf.b;
	  g->b += w*(gcf.a*GFS_VARIABLE (f.cell, v) - gcf.c);
	}
#if (!FTT_2D && !FTT_2D3)
      g->a /= 2.;
      g->b /= 2.;
#endif /* not 2D and not 2D3 */
    }
  }
}

static void fullest_directions (const FttCellFace * face,
				FttDirection d[FTT_DIMENSION])
{
  FttComponent c = face->d/2, i;
  FttCell * mixed = GFS_IS_MIXED (face->cell) ? face->cell : face->neighbor;
  GfsSolidVector * s = GFS_STATE (mixed)->solid;

  d[0] = face->d;
  for (i = 1; i < FTT_DIMENSION; i++) {
    guint cp = (c + i) % FTT_DIMENSION;
    d[i] = s->s[2*cp] > s->s[2*cp + 1] ? 2*cp : 2*cp + 1;
  }
}

static FttCell * cell_corner_neighbor1 (FttCell * cell,
					FttDirection * d,
					gint max_level)
{
  FttCell * neighbor = ftt_cell_neighbor (cell, d[0]);
  if (!neighbor)
    return NULL;
  else {
    guint level = ftt_cell_level (cell);
    if (ftt_cell_level (neighbor) < level)
      /* neighbor is at a shallower level */
      return neighbor;
    else {
      if (level == max_level || FTT_CELL_IS_LEAF (neighbor))
	/* neighbor is at the same level */
	return neighbor;
      else {
	/* neighbor is at a deeper level */
	guint i;
	FttDirection d1[FTT_DIMENSION];
	d1[0] = FTT_OPPOSITE_DIRECTION (d[0]);
	for (i = 1; i < FTT_DIMENSION; i++)
	  d1[i] = d[i];
	return ftt_cell_child_corner (neighbor, d1);
      }
    }
  }
}

#if FTT_2D
# define N_CELLS 4
#else  /* 2D3 and 3D */
# define N_CELLS 8
#endif /* 2D3 and 3D */

static gboolean inverse (gdouble mi[N_CELLS - 1][N_CELLS - 1])
{
#if FTT_2D
  gdouble m[N_CELLS - 1][N_CELLS - 1], det;
  guint i, j;
  
  for (i = 0; i < N_CELLS - 1; i++)
    for (j = 0; j < N_CELLS - 1; j++)
      m[i][j] = mi[i][j];
  
  det = (m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2]) -
	 m[0][1]*(m[1][0]*m[2][2] - m[2][0]*m[1][2]) +
	 m[0][2]*(m[1][0]*m[2][1] - m[2][0]*m[1][1]));
  if (det == 0.)
    return FALSE;
  
  mi[0][0] = (m[1][1]*m[2][2] - m[1][2]*m[2][1])/det; 
  mi[0][1] = (m[2][1]*m[0][2] - m[0][1]*m[2][2])/det;
  mi[0][2] = (m[0][1]*m[1][2] - m[1][1]*m[0][2])/det; 
  mi[1][0] = (m[1][2]*m[2][0] - m[1][0]*m[2][2])/det; 
  mi[1][1] = (m[0][0]*m[2][2] - m[2][0]*m[0][2])/det; 
  mi[1][2] = (m[1][0]*m[0][2] - m[0][0]*m[1][2])/det; 
  mi[2][0] = (m[1][0]*m[2][1] - m[2][0]*m[1][1])/det; 
  mi[2][1] = (m[2][0]*m[0][1] - m[0][0]*m[2][1])/det; 
  mi[2][2] = (m[0][0]*m[1][1] - m[0][1]*m[1][0])/det; 
#else /* 3D */
  gint indxc[N_CELLS - 1], indxr[N_CELLS - 1], ipiv[N_CELLS - 1];
  gint i, icol = 0, irow = 0, j, k, l, ll;
  gdouble big, dum, pivinv, temp;
  
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

  for (j = 0; j < N_CELLS - 1; j++) 
    ipiv[j] = -1;

  for (i = 0; i < N_CELLS - 1; i++) {
    big = 0.0;
    for (j = 0; j < N_CELLS - 1; j++)
      if (ipiv[j] != 0)
	for (k = 0; k < N_CELLS - 1; k++) {
	  if (ipiv[k] == -1) {
	    if (fabs (mi[j][k]) >= big) {
	      big = fabs (mi[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }
	}
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < N_CELLS - 1; l++) 
	SWAP (mi[irow][l], mi[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (mi[icol][icol] == 0.)
      return FALSE;
    pivinv = 1.0/mi[icol][icol];
    mi[icol][icol] = 1.0;
    for (l = 0; l < N_CELLS - 1; l++) mi[icol][l] *= pivinv;
    for (ll = 0; ll < N_CELLS - 1; ll++)
      if (ll != icol) {
	dum = mi[ll][icol];
	mi[ll][icol] = 0.0;
	for (l = 0; l < N_CELLS - 1; l++)
	  mi[ll][l] -= mi[icol][l]*dum;
      }
  }
  for (l = N_CELLS - 1 - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < N_CELLS - 1; k++)
	SWAP (mi[k][indxr[l]], mi[k][indxc[l]]);
  }
#endif /* 3D */
  return TRUE;
}

#if (!FTT_2D)
static void draw_cell (FttCell * cell, gdouble r, gdouble g, gdouble b,
		       const gchar * name)
{
  FttVector p;
  gdouble size = ftt_cell_size (cell)/2.;

  ftt_cell_pos (cell, &p);
  fprintf (stderr,
	   "(geometry \"%s\" = OFF 8 6 12\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "4 3 2 1 0 %g %g %g\n"
	   "4 4 5 6 7 %g %g %g\n"
	   "4 2 3 7 6 %g %g %g\n"
	   "4 0 1 5 4 %g %g %g\n"
	   "4 0 4 7 3 %g %g %g\n"
	   "4 1 2 6 5 %g %g %g\n"
	   ")\n",
	   name,
	   p.x - size, p.y - size, p.z - size,
	   p.x + size, p.y - size, p.z - size,
	   p.x + size, p.y + size, p.z - size,
	   p.x - size, p.y + size, p.z - size,
	   p.x - size, p.y - size, p.z + size,
	   p.x + size, p.y - size, p.z + size,
	   p.x + size, p.y + size, p.z + size,
	   p.x - size, p.y + size, p.z + size,
	   r, g, b,
	   r, g, b,
	   r, g, b,
	   r, g, b,
	   r, g, b,
	   r, g, b);
  gfs_cell_cm (cell, &p);
  size /= 8.;
  fprintf (stderr,
	   "(geometry \"cm %s\" = OFF 8 6 12\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "4 3 2 1 0 %g %g %g\n"
	   "4 4 5 6 7 %g %g %g\n"
	   "4 2 3 7 6 %g %g %g\n"
	   "4 0 1 5 4 %g %g %g\n"
	   "4 0 4 7 3 %g %g %g\n"
	   "4 1 2 6 5 %g %g %g\n"
	   ")\n",
	   name,
	   p.x - size, p.y - size, p.z - size,
	   p.x + size, p.y - size, p.z - size,
	   p.x + size, p.y + size, p.z - size,
	   p.x - size, p.y + size, p.z - size,
	   p.x - size, p.y - size, p.z + size,
	   p.x + size, p.y - size, p.z + size,
	   p.x + size, p.y + size, p.z + size,
	   p.x - size, p.y + size, p.z + size,
	   r, g, b,
	   r, g, b,
	   r, g, b,
	   r, g, b,
	   r, g, b,
	   r, g, b);
}

static void output_error_mesh (FttCell ** n)
{
  draw_cell (n[0], 0., 0., 0., "n0");
  draw_cell (n[1], 0.5, 0., 0., "n1");
  draw_cell (n[2], 0., 0.5, 0., "n2");
  draw_cell (n[3], 0., 0., 0.5, "n3");
  draw_cell (n[4], 1., 1., 1., "n4");
  draw_cell (n[5], 1., 0., 0., "n5");
  draw_cell (n[6], 0., 1., 0., "n6");
  draw_cell (n[7], 0., 0., 1., "n7");
  g_assert_not_reached ();
}
#endif /* 3D */

static gboolean face_bilinear (const FttCellFace * face,
			       FttCell ** n,
			       FttVector * o,
			       void (*cell_pos) (const FttCell *, FttVector *),
			       gint max_level,
			       gdouble m[N_CELLS - 1][N_CELLS - 1])
{
  FttDirection d[3], d1[3];
  guint i;
  gdouble size = ftt_cell_size (face->cell);

  fullest_directions (face, d);
  n[0] = face->cell; n[1] = face->neighbor;
  d1[0] = d[1]; d1[1] = d[0]; d1[2] = d[2];
  if ((n[2] = cell_corner_neighbor1 (n[0], d1, max_level)) == NULL)
    return FALSE;
  d1[1] = FTT_OPPOSITE_DIRECTION (d[0]);
  if ((n[3] = cell_corner_neighbor1 (n[1], d1, max_level)) == NULL)
    return FALSE;
  if (n[3] == n[2]) {
    d1[0] = d[0]; d1[1] = FTT_OPPOSITE_DIRECTION (d[1]);
    if ((n[3] = cell_corner_neighbor1 (n[2], d1, max_level)) == NULL)
      return FALSE;
  }

#if FTT_2D
  for (i = 0; i < 3; i++) {
    FttVector cm;
    guint j;

    (*cell_pos) (n[i + 1], &cm);

    for (j = 0; j < FTT_DIMENSION; j++) {
      (&cm.x)[j] -= (&o->x)[j];
      /* fixme: this does not work for periodic boundaries */      
      (&cm.x)[j] /= size;
      g_assert (fabs ((&cm.x)[j]) <= 4.);
    }

    m[i][0] = cm.x;
    m[i][1] = cm.y; 
    m[i][2] = cm.x*cm.y;
  }
  g_assert (inverse (m));
#else /* 3D */
  d1[0] = d[2]; d1[1] = d[0]; d1[2] = d[1];
  if ((n[4] = cell_corner_neighbor1 (n[0], d1, max_level)) == NULL)
    return FALSE;
  d1[1] = FTT_OPPOSITE_DIRECTION (d[0]);
  if ((n[5] = cell_corner_neighbor1 (n[1], d1, max_level)) == NULL)
    return FALSE;
  if (n[5] == n[4]) {
    d1[0] = d[0]; d1[1] = d[1]; d1[2] = FTT_OPPOSITE_DIRECTION (d[2]);
    if ((n[5] = cell_corner_neighbor1 (n[4], d1, max_level)) == NULL)
      return FALSE;
  }
  d1[0] = d[2]; d1[1] = d[0]; d1[2] = FTT_OPPOSITE_DIRECTION (d[1]);
  if ((n[6] = cell_corner_neighbor1 (n[2], d1, max_level)) == NULL)
    return FALSE;
  if (n[6] == n[4]) {
    d1[0] = d[1]; d1[1] = d[0]; d1[2] = FTT_OPPOSITE_DIRECTION (d[2]);
    if ((n[6] = cell_corner_neighbor1 (n[4], d1, max_level)) == NULL)
      return FALSE;
  }
  d1[0] = d[2]; d1[1] = FTT_OPPOSITE_DIRECTION (d[0]); d1[2] = FTT_OPPOSITE_DIRECTION (d[1]);
  if ((n[7] = cell_corner_neighbor1 (n[3], d1, max_level)) == NULL)
    return FALSE;
  if (n[7] == n[4] || n[7] == n[5]) {
    d1[0] = d[1]; d1[2] = FTT_OPPOSITE_DIRECTION (d[2]);
    if ((n[7] = cell_corner_neighbor1 (n[5], d1, max_level)) == NULL)
      return FALSE;
  }
  if (n[7] == n[6]) {
    d1[0] = d[0]; d1[1] = FTT_OPPOSITE_DIRECTION (d[1]); d1[2] = FTT_OPPOSITE_DIRECTION (d[2]);
    if ((n[7] = cell_corner_neighbor1 (n[6], d1, max_level)) == NULL)
      return FALSE;
  }

  for (i = 0; i < 7; i++) {
    FttVector cm;
    guint j;

    for (j = i + 1; j < 7; j++)
      if (n[i + 1] == n[j + 1])
	output_error_mesh (n);

    (*cell_pos) (n[i + 1], &cm);

    for (j = 0; j < FTT_DIMENSION; j++) {
      (&cm.x)[j] -= (&o->x)[j];
      /* fixme: this does not work for periodic boundaries */      
      (&cm.x)[j] /= size;
      if (fabs ((&cm.x)[j]) > 4.)
	output_error_mesh (n);
    }

    m[i][0] = cm.x;
    m[i][1] = cm.y; 
    m[i][2] = cm.z;
    m[i][3] = cm.x*cm.y;
    m[i][4] = cm.x*cm.z;
    m[i][5] = cm.y*cm.z;
    m[i][6] = cm.x*cm.y*cm.z;
  }
  if (!inverse (m))
    output_error_mesh (n);
#endif /* 3D */

  return TRUE;
}

/* grad(v) = = -a*v(cell) + b*v(neighbor) + c */
static gboolean mixed_face_gradient (const FttCellFace * face,
				     Gradient * g,
				     guint v,
				     gint max_level)
{
  FttCell * n[N_CELLS];
  gdouble m[N_CELLS - 1][N_CELLS - 1];
  FttVector o, cm;
  FttComponent c = face->d/2;
  gdouble h = ftt_cell_size (face->cell);

  gfs_cell_cm (face->cell, &o);
  if (!face_bilinear (face, n, &o, gfs_cell_cm, max_level, m))
    return FALSE;

  gfs_face_ca (face, &cm);

#if FTT_2D
  {
    FttComponent cp = FTT_ORTHOGONAL_COMPONENT (c);
    gdouble vp;
    
    vp = ((&cm.x)[cp] - (&o.x)[cp])/h;
    g->a = ((m[c][0] + m[2][0]*vp) +
	    (m[c][1] + m[2][1]*vp) +
	    (m[c][2] + m[2][2]*vp));
    g->b = m[c][0] + m[2][0]*vp;
    g->c = ((m[c][1] + m[2][1]*vp)*GFS_VARIABLE (n[2], v) +
	    (m[c][2] + m[2][2]*vp)*GFS_VARIABLE (n[3], v));
  }
#else /* 3D */
  {
    guint j;

    cm.x = (cm.x - o.x)/h;
    cm.y = (cm.y - o.y)/h;
    cm.z = (cm.z - o.z)/h;
    g->c = 0.;
    
    switch (c) {
    case FTT_X:
      g->a = g->b = m[0][0] + cm.y*m[3][0] + cm.z*m[4][0] + cm.y*cm.z*m[6][0];
      for (j = 1; j < N_CELLS - 1; j++) {
	gdouble a = m[0][j] + cm.y*m[3][j] + cm.z*m[4][j] + cm.y*cm.z*m[6][j];
	g->a += a;
	g->c += a*GFS_VARIABLE (n[j+1], v);
      }
      break;
    case FTT_Y:
      g->a = g->b = m[1][0] + cm.x*m[3][0] + cm.z*m[5][0] + cm.x*cm.z*m[6][0];
      for (j = 1; j < N_CELLS - 1; j++) {
	gdouble a = m[1][j] + cm.x*m[3][j] + cm.z*m[5][j] + cm.x*cm.z*m[6][j];
	g->a += a;
	g->c += a*GFS_VARIABLE (n[j+1], v);
      }
      break;
    case FTT_Z:
      g->a = g->b = m[2][0] + cm.x*m[4][0] + cm.y*m[5][0] + cm.x*cm.y*m[6][0];
      for (j = 1; j < N_CELLS - 1; j++) {
	gdouble a = m[2][j] + cm.x*m[4][j] + cm.y*m[5][j] + cm.x*cm.y*m[6][j];
	g->a += a;
	g->c += a*GFS_VARIABLE (n[j+1], v);
      }
      break;
    default:
      g_assert_not_reached ();
    }
  }
#endif /* 3D */

  if (!FTT_FACE_DIRECT (face)) {
    g->a = - g->a;
    g->b = - g->b;
    g->c = - g->c;
  }

  return TRUE;
}

/**
 * gfs_face_gradient_flux:
 * @face: a #FttCellFace.
 * @g: the #GfsGradient.
 * @v: a #GfsVariable index.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 *
 * Set the value of @g as the gradient of variable @v on the @face
 * weighted by the value of the @v field of the face state vector of
 * the corresponding cell. Variable @v is defined at the center of
 * mass of its cell. Linear interpolation is used to evaluate the
 * gradient in the vicinity of cut cells.
 */
void gfs_face_gradient_flux (const FttCellFace * face,
			     GfsGradient * g,
			     guint v,
			     gint max_level)
{
  guint level;
  Gradient gcf;
  gdouble w;

  g_return_if_fail (face != NULL);
  g_return_if_fail (g != NULL);

  g->a = g->b = 0.;
  if (face->neighbor == NULL || (w = GFS_STATE (face->cell)->f[face->d].v) == 0.)
    return;

  level = ftt_cell_level (face->cell);
  if (ftt_cell_level (face->neighbor) < level) {
    /* neighbor is at a shallower level */
    if (GFS_IS_MIXED (face->cell) || GFS_IS_MIXED (face->neighbor)) {
      if (!mixed_face_gradient (face, &gcf, v, max_level))
	gcf = gradient_fine_coarse (face, v, max_level);
    }
    else
      gcf = gradient_fine_coarse (face, v, max_level);
    g->a = w*gcf.a;
    g->b = w*(gcf.b*GFS_VARIABLE (face->neighbor, v) + gcf.c);
  }
  else {
    if (level == max_level || FTT_CELL_IS_LEAF (face->neighbor)) {
      /* neighbor is at the same level */
      if (!GFS_IS_MIXED (face->cell) && !GFS_IS_MIXED (face->neighbor)) {
	g->a = w;
	g->b = w*GFS_VARIABLE (face->neighbor, v);
      }
      else if (mixed_face_gradient (face, &gcf, v, max_level)) {
	g->a = w*gcf.a;
	g->b = w*(gcf.b*GFS_VARIABLE (face->neighbor, v) + gcf.c);
      }
      else {
	g->a = w;
	g->b = w*GFS_VARIABLE (face->neighbor, v);
      }
    }
    else {
      /* neighbor is at a deeper level */
      FttCellChildren children;
      FttCellFace f;
      guint i, n;
      
      f.d = FTT_OPPOSITE_DIRECTION (face->d);
      n = ftt_cell_children_direction (face->neighbor, f.d, &children);
      f.neighbor = face->cell;
      for (i = 0; i < n; i++) 
	if ((f.cell = children.c[i])) {
	  w = GFS_STATE (f.cell)->f[f.d].v;
	  if (GFS_IS_MIXED (f.cell) || GFS_IS_MIXED (f.neighbor)) {
	    if (!mixed_face_gradient (&f, &gcf, v, max_level))
	      gcf = gradient_fine_coarse (&f, v, max_level);
	  }
	  else
	    gcf = gradient_fine_coarse (&f, v, max_level);
	  g->a += w*gcf.b;
	  g->b += w*(gcf.a*GFS_VARIABLE (f.cell, v) - gcf.c);
	}
    }
  }
}

static gboolean cell_bilinear (FttCell * cell,
			       FttCell ** n,
			       FttVector * o,
			       void (*cell_pos) (const FttCell *, FttVector *),
			       gint max_level,
			       gdouble m[N_CELLS - 1][N_CELLS - 1])
{
  GfsSolidVector * s = GFS_STATE (cell)->solid;  
  FttCellFace f;
  FttDirection d[FTT_DIMENSION];
  FttComponent c;

  if ((s->s[FTT_RIGHT] == 0. && s->s[FTT_LEFT] == 0.) ||
      (s->s[FTT_RIGHT] == 1. && s->s[FTT_LEFT] == 1.))
    return FALSE;

  for (c = 0; c < FTT_DIMENSION; c++)
    d[c] = s->s[2*c] > s->s[2*c + 1] ? 2*c : 2*c + 1;
  f.cell = cell;
  f.d = d[0];
  f.neighbor = cell_corner_neighbor1 (cell, d, max_level);

  return face_bilinear (&f, n, o, cell_pos, max_level, m);
}

/**
 * gfs_cell_dirichlet_gradient:
 * @cell: a #FttCell.
 * @v: a #GfsVariable index.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 * @v0: the Dirichlet value on the boundary.
 * @grad: a #FttVector.
 *
 * Fills @grad with components of the gradient of variable @v
 * interpolated at the center of area of the solid boundary contained
 * in @cell. The gradient is scaled by the size of the cell.
 */
void gfs_cell_dirichlet_gradient (FttCell * cell,
				  guint v,
				  gint max_level,
				  gdouble v0,
				  FttVector * grad)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (grad != NULL);

  if (!GFS_IS_MIXED (cell))
    return;
  else {
    FttCell * n[N_CELLS];
    gdouble m[N_CELLS - 1][N_CELLS - 1];
    guint i, c;

    grad->x = grad->y = grad->z = 0.;
    if (!cell_bilinear (cell, n, &GFS_STATE (cell)->solid->ca, 
			gfs_cell_cm, max_level, m))
      return;

    for (i = 0; i < N_CELLS - 1; i++)
      for (c = 0; c < FTT_DIMENSION; c++)
	(&grad->x)[c] += m[c][i]*(GFS_VARIABLE (n[i + 1], v) - v0);
  }
}

/**
 * gfs_cell_dirichlet_gradient_flux:
 * @cell: a #FttCell.
 * @v: a #GfsVariable index.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 * @v0: the Dirichlet value on the boundary.
 *
 * Returns: the flux of the gradient of variable @v through the solid
 * boundary contained in @cell.
 */
gdouble gfs_cell_dirichlet_gradient_flux (FttCell * cell,
					  guint v,
					  gint max_level,
					  gdouble v0)
{
  g_return_val_if_fail (cell != NULL, 0.);

  if (!GFS_IS_MIXED (cell))
    return 0.;
  else {
    GfsSolidVector * s = GFS_STATE (cell)->solid;
    FttVector g;
    
    gfs_cell_dirichlet_gradient (cell, v, max_level, v0, &g);

    return (g.x*(s->s[1] - s->s[0]) + g.y*(s->s[3] - s->s[2])
#if (!FTT_2D)
	    + g.z*(s->s[5] - s->s[4])
#endif
	    )*s->v;
  }
}

/**
 * gfs_cell_dirichlet_value:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 *
 * Returns: the value of variable @v interpolated at the center of
 * area of the solid boundary contained in @cell.
 */
gdouble gfs_cell_dirichlet_value (FttCell * cell,
				  GfsVariable * v,
				  gint max_level)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  if (!GFS_IS_MIXED (cell))
    return 0.;
  else {
    FttCell * n[N_CELLS];
    FttVector p;
    gdouble m[N_CELLS - 1][N_CELLS - 1], a[N_CELLS - 1];
    GfsSolidVector * s = GFS_STATE (cell)->solid;
    gdouble v0, size = ftt_cell_size (cell);
    void (* cell_pos) (const FttCell *, FttVector *) = v->centered ? ftt_cell_pos : gfs_cell_cm;
    guint i, j;

    (*cell_pos) (cell, &p);
    if (!cell_bilinear (cell, n, &p, cell_pos, max_level, m))
      return 0.;

    v0 = GFS_VARIABLE (cell, v->i);
    p.x = (s->ca.x - p.x)/size;
    p.y = (s->ca.y - p.y)/size;
    for (i = 0; i < N_CELLS - 1; i++) {
      a[i] = 0.;
      for (j = 0; j < N_CELLS - 1; j++)
	a[i] += m[i][j]*(GFS_VARIABLE (n[j + 1], v->i) - v0);
    }
#if FTT_2D
    return (a[0]*p.x + a[1]*p.y + a[2]*p.x*p.y + v0);
#else /* 3D */
    p.z = (s->ca.z - p.z)/size;
    return (a[0]*p.x + a[1]*p.y + a[2]*p.z + 
	    a[3]*p.x*p.y + a[4]*p.x*p.z + a[5]*p.y*p.z + 
	    a[6]*p.x*p.y*p.z + v0);
#endif /* 3D */
  }
}

/**
 * gfs_get_from_above:
 * @cell: a #FttCell.
 * @v: a #GfsVariable to "get from above".
 *
 * Sets the value of the variable @v of @cell to the value of this
 * variable in its parent cell.
 *
 * This function fails if @cell is the root of the cell tree.
 */
void gfs_get_from_above (FttCell * cell, const GfsVariable * v)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (!FTT_CELL_IS_ROOT (cell));
  g_return_if_fail (v != NULL);

  GFS_VARIABLE (cell, v->i) = GFS_VARIABLE (ftt_cell_parent (cell), v->i);
}

/**
 * gfs_get_from_below_intensive:
 * @cell: a #FttCell.
 * @v: a #GfsVariable to "get from below".
 *
 * Sets the value of the "intensive" variable @v of @cell by taking
 * the volume weighted average of the values of its children cells.
 *
 * This functions fails if @cell is a leaf of the cell tree.
 */
void gfs_get_from_below_intensive (FttCell * cell, const GfsVariable * v)
{
  gdouble val = 0., sa = 0.;
  guint i;
  FttCellChildren child;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (!FTT_CELL_IS_LEAF (cell));
  g_return_if_fail (v != NULL);

  ftt_cell_children (cell, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
      gdouble a = GFS_IS_MIXED (child.c[i]) ? 
	GFS_STATE (child.c[i])->solid->a : 1.;

      val += GFS_VARIABLE (child.c[i], v->i)*a;
      sa += a;
    }
  if (sa > 0.)
    GFS_VARIABLE (cell, v->i) = val/sa;
}

/**
 * gfs_get_from_below_extensive:
 * @cell: a #FttCell.
 * @v: a #GfsVariable to "get from below".
 *
 * Sets the value of the "extensive" variable @v of @cell as (half in
 * 3D) the sum of the values of its children cells.
 *
 * This functions fails if @cell is a leaf of the cell tree.  
 */
void gfs_get_from_below_extensive (FttCell * cell, const GfsVariable * v)
{
  gdouble val = 0.;
  guint i;
  FttCellChildren child;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (!FTT_CELL_IS_LEAF (cell));
  g_return_if_fail (v != NULL);

  ftt_cell_children (cell, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i])
      val += GFS_VARIABLE (child.c[i], v->i);
  GFS_VARIABLE (cell, v->i) = val
#if (!FTT_2D && !FTT_2D3)
    /2.
#endif /* not FTT_2D */
    ;
}

/**
 * gfs_cell_cleanup:
 * @cell: a #FttCell.
 *
 * Frees the memory allocated for extra data associated with @cell.
 *
 * This function must be used as "cleanup function" when using
 * ftt_cell_destroy().
 */
void gfs_cell_cleanup (FttCell * cell)
{
  g_return_if_fail (cell != NULL);
  
  if (cell->data && GFS_STATE (cell)->solid) {
    g_free (GFS_STATE (cell)->solid);
    GFS_STATE (cell)->solid = NULL;
  }
  g_free (cell->data);
}

/**
 * gfs_cell_reset:
 * @cell: a #FttCell.
 * @v: a #GfsVariable to reset.
 *
 * Sets the value of the variable @v of @cell to zero.
 */
void gfs_cell_reset (FttCell * cell, GfsVariable * v)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);

  GFS_VARIABLE (cell, v->i) = 0.;
}

static void add_stats (const FttCell * cell, gpointer * data)
{
  GtsRange * s = data[0];
  GfsVariable * v = data[1];

  gts_range_add_value (s, GFS_VARIABLE (cell, v->i));
}

/**
 * gfs_stats_variable:
 * @root: the root #FttCell of the tree to obtain statistics from.
 * @v: the variable to consider for statistics.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 *
 * Traverses the cell tree defined by @root using ftt_cell_traverse()
 * and gathers statistics about variable @v.
 *
 * Returns: a #GtsRange containing the statistics about @v.
 */
GtsRange gfs_stats_variable (FttCell * root,
			     GfsVariable * v,
			     FttTraverseFlags flags,
			     gint max_depth)
{
  GtsRange s;
  gpointer data[2];

  g_return_val_if_fail (root != NULL, s);
  g_return_val_if_fail (v != NULL, s);
  
  gts_range_init (&s);
  data[0] = &s;
  data[1] = v;
  ftt_cell_traverse (root, FTT_PRE_ORDER, flags, max_depth, 
		     (FttCellTraverseFunc) add_stats, data);
  gts_range_update (&s);

  return s;
}

static void add_norm (const FttCell * cell, gpointer * data)
{
  GfsNorm * n = data[0];
  GfsVariable * v = data[1];

  gfs_norm_add (n, GFS_VARIABLE (cell, v->i), 
		ftt_cell_volume (cell)*(GFS_IS_MIXED (cell) ? 
				    GFS_STATE (cell)->solid->a : 1.));
}

/**
 * gfs_norm_variable:
 * @root: the root #FttCell of the tree to obtain norm from.
 * @v: the variable to consider for norm statistics.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 *
 * Traverses the cell tree defined by @root using ftt_cell_traverse()
 * and gathers norm statistics about variable @v.
 *
 * Returns: a #GfsNorm containing the norm statistics about @v.
 */
GfsNorm gfs_norm_variable (FttCell * root,
			   GfsVariable * v,
			   FttTraverseFlags flags,
			   gint max_depth)
{
  GfsNorm n;
  gpointer data[2];

  g_return_val_if_fail (root != NULL, n);
  g_return_val_if_fail (v != NULL, n);
  
  gfs_norm_init (&n);
  data[0] = &n;
  data[1] = v;
  ftt_cell_traverse (root, FTT_PRE_ORDER, flags, max_depth, 
		     (FttCellTraverseFunc) add_norm, data);
  gfs_norm_update (&n);

  return n;
}

/**
 * gfs_norm_init:
 * @n: a #GfsNorm.
 *
 * Initializes a #GfsNorm.
 */
void gfs_norm_init (GfsNorm * n)
{
  g_return_if_fail (n != NULL);

  n->bias = n->first = n->second = 0.;
  n->infty = - G_MAXDOUBLE;
  n->w = 0.;
}

/**
 * gfs_norm_reset:
 * @n: a #GfsNorm.
 *
 * Sets all the fields of @n to 0.
 */
void gfs_norm_reset (GfsNorm * n)
{
  g_return_if_fail (n != NULL);

  n->bias = n->first = n->second = 0.;
  n->infty = 0.;
  n->w = 0.;
}

/**
 * gfs_norm_add:
 * @n: a #GfsNorm.
 * @val: a value to add to @n.
 * @weight: weight of @val.
 *
 * Adds @val to @n.
 */
void gfs_norm_add (GfsNorm * n, gdouble val, gdouble weight)
{
  g_return_if_fail (n != NULL);

  n->bias += weight*val;
  val = fabs (val);
  if (val > n->infty)
    n->infty = val;
  n->first += weight*val;
  n->second += weight*val*val;
  n->w += weight;
}

/**
 * gfs_norm_update:
 * @n: a #GfsNorm.
 * 
 * Updates the fields of @n.
 */
void gfs_norm_update (GfsNorm * n)
{
  g_return_if_fail (n != NULL);

  if (n->w > 0.0) {
    n->bias /= n->w;
    n->first /= n->w;
    n->second = sqrt (n->second/n->w);
  }
  else
    n->infty = 0.0;
}

/**
 * gfs_face_interpolated_value:
 * @face: a #FttFace.
 * @v: a #GfsVariable index.
 *
 * Computes the value of variable @v on the @face using second-order
 * interpolation from the cell-centered values.
 *
 * Returns: the value of variable @v on the face.  
 */
gdouble gfs_face_interpolated_value (const FttCellFace * face,
				     guint v)
{
  gdouble x1 = 1., v1;
#if 1
  g_return_val_if_fail (face != NULL, 0.);

  if (face->neighbor) {
    v1 = neighbor_value (face, v, &x1);
    return ((x1 - 0.5)*GFS_VARIABLE (face->cell, v) + 0.5*v1)/x1;
  }
  else
    return GFS_VARIABLE (face->cell, v);
#else
  gdouble v0;
  FttCellFace f2;

  g_return_val_if_fail (face != NULL, 0.);

  v0 = GFS_VARIABLE (face->cell, v);
  v1 = neighbor_value (face, v, &x1);
  f2 = ftt_cell_face (face->cell, FTT_OPPOSITE_DIRECTION (face->d));
  if (f2.neighbor) {
    gdouble x2 = 1.;
    gdouble v2 = neighbor_value (&f2, v, &x2);

    return v0 + (x2*(v1 - v0)*(1. + 2.*x2) - x1*(v0 - v2)*(1. - 2.*x1))
      /(4.*x1*x2*(x1 + x2));
  }
  else
    return ((x1 - 0.5)*v0 + 0.5*v1)/x1;
#endif
}

/**
 * gfs_normal_divergence:
 * @cell: a #FttCell.
 *
 * Fills variable %GFS_DIV of @cell with the integral of the divergence
 * of the (MAC) velocity field in this cell.  
 */
void gfs_normal_divergence (FttCell * cell)
{
  FttComponent c;
  gdouble div = 0.;

  g_return_if_fail (cell != NULL);

  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * solid = GFS_STATE (cell)->solid;
    
    for (c = 0; c < FTT_DIMENSION; c++) {
      FttDirection d = 2*c;
      
      div += (solid->s[d]*GFS_STATE (cell)->f[d].un - 
	      solid->s[d + 1]*GFS_STATE (cell)->f[d + 1].un);
    }
  }
  else
    for (c = 0; c < FTT_DIMENSION; c++) {
      FttDirection d = 2*c;
      
      div += (GFS_STATE (cell)->f[d].un - 
	      GFS_STATE (cell)->f[d + 1].un);
    }
  GFS_STATE (cell)->div = div*ftt_cell_size (cell);
}

/**
 * gfs_normal_divergence_2D:
 * @cell: a #FttCell.
 *
 * Fills variable %GFS_DIV of @cell with the integral of the 2D
 * divergence of the (MAC) velocity field in this cell.
 */
void gfs_normal_divergence_2D (FttCell * cell)
{
  FttComponent c;
  gdouble div = 0.;

  g_return_if_fail (cell != NULL);

  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * solid = GFS_STATE (cell)->solid;
    
    for (c = 0; c < 2; c++) {
      FttDirection d = 2*c;
      
      div += (solid->s[d]*GFS_STATE (cell)->f[d].un - 
	      solid->s[d + 1]*GFS_STATE (cell)->f[d + 1].un);
    }
  }
  else
    for (c = 0; c < 2; c++) {
      FttDirection d = 2*c;
      
      div += (GFS_STATE (cell)->f[d].un - 
	      GFS_STATE (cell)->f[d + 1].un);
    }
  GFS_STATE (cell)->div = div*ftt_cell_size (cell);
}

/**
 * gfs_divergence:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 *
 * Fills variable @v of @cell with the divergence of the
 * (centered) velocity field in this cell.  
 */
void gfs_divergence (FttCell * cell, GfsVariable * v)
{
  FttComponent c;
  gdouble div = 0.;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);

  for (c = 0; c < FTT_DIMENSION; c++)
    div += gfs_center_gradient (cell, c, GFS_VELOCITY_INDEX (c));
  GFS_VARIABLE (cell, v->i) = div/ftt_cell_size (cell);
}

/**
 * gfs_vorticity_value:
 * @cell: a #FttCell.
 * @lambda: the dimensions of the domain containing @cell.
 *
 * Returns: the vorticity (norm of the vorticity vector in 3D) of the
 * velocity field in this cell.
 */
gdouble gfs_vorticity_value (FttCell * cell, FttVector * lambda)
{
  gdouble size;
#if (!FTT_2D)
  FttVector vort;
#endif /* FTT_3D */

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (lambda != NULL, 0.);

  size = ftt_cell_size (cell);
#if FTT_2D
  return (lambda->x*gfs_center_gradient (cell, FTT_X, GFS_V)/lambda->y -
	  lambda->y*gfs_center_gradient (cell, FTT_Y, GFS_U)/lambda->x)/size;
#else  /* FTT_3D */
  vort.x = (lambda->y*gfs_center_gradient (cell, FTT_Y, GFS_W)/lambda->z -
	    lambda->z*gfs_center_gradient (cell, FTT_Z, GFS_V)/lambda->y)/size;
  vort.y = (lambda->z*gfs_center_gradient (cell, FTT_Z, GFS_U)/lambda->x -
	    lambda->x*gfs_center_gradient (cell, FTT_X, GFS_W)/lambda->z)/size;
  vort.z = (lambda->x*gfs_center_gradient (cell, FTT_X, GFS_V)/lambda->y -
	    lambda->y*gfs_center_gradient (cell, FTT_Y, GFS_U)/lambda->x)/size;
  return sqrt (vort.x*vort.x + vort.y*vort.y + vort.z*vort.z);
#endif /* FTT_3D */
}

/**
 * gfs_vorticity:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 *
 * Fills variable @v of @cell with the vorticity (norm of the
 * vorticity vector in 3D) of the velocity field in this cell.  
 */
void gfs_vorticity (FttCell * cell,
		    GfsVariable * v)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);

  GFS_VARIABLE (cell, v->i) = gfs_vorticity_value (cell, 
			&GFS_DOMAIN (gfs_variable_parent (v))->lambda);
}

/**
 * gfs_velocity_norm:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 *
 * Fills variable @v of @cell with the norm of the velocity field in
 * this cell.
 */
void gfs_velocity_norm (FttCell * cell,
			GfsVariable * v)
{
  GfsStateVector * s;
  
  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);

  s = GFS_STATE (cell);
#if FTT_2D
  GFS_VARIABLE (cell, v->i) = sqrt (s->u*s->u + s->v*s->v);
#else  /* FTT_3D */
  GFS_VARIABLE (cell, v->i) = sqrt (s->u*s->u + s->v*s->v + s->w*s->w);
#endif /* FTT_3D */
}

/**
 * gfs_velocity_norm2:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 *
 * Fills variable @v of @cell with the squared norm of the velocity field in
 * this cell.
 */
void gfs_velocity_norm2 (FttCell * cell,
			 GfsVariable * v)
{
  GfsStateVector * s;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);

  s = GFS_STATE (cell);
#if FTT_2D
  GFS_VARIABLE (cell, v->i) = s->u*s->u + s->v*s->v;
#else  /* FTT_3D */
  GFS_VARIABLE (cell, v->i) = s->u*s->u + s->v*s->v + s->w*s->w;
#endif /* FTT_3D */
}

/**
 * gfs_velocity_lambda2:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 *
 * Fills variable @v of @cell with the lambda2 eigenvalue used by
 * Jeong and Hussain as vortex criterion (JFM 285, 69-94, 1995).
 *
 * The value is normalized by the square of the size of the cell.
 */
void gfs_velocity_lambda2 (FttCell * cell,
			   GfsVariable * v)
{
  gdouble J[FTT_DIMENSION][FTT_DIMENSION];
  gdouble S2O2[FTT_DIMENSION][FTT_DIMENSION];
  gdouble lambda[FTT_DIMENSION], ev[FTT_DIMENSION][FTT_DIMENSION];
  guint i, j, k;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);

  for (i = 0; i < FTT_DIMENSION; i++)
    for (j = 0; j < FTT_DIMENSION; j++)
      J[i][j] = gfs_center_gradient (cell, j, GFS_VELOCITY_INDEX (i));
  for (i = 0; i < FTT_DIMENSION; i++)
    for (j = 0; j < FTT_DIMENSION; j++) {
      S2O2[i][j] = 0.;
      for (k = 0; k < FTT_DIMENSION; k++)
	S2O2[i][j] += J[i][k]*J[k][j] + J[k][i]*J[j][k];
    }
  gfs_eigenvalues (S2O2, lambda, ev);
  GFS_VARIABLE (cell, v->i) = lambda[1]/2.;
}

static void cell_traverse_mixed (FttCell * cell,
				 FttTraverseType order,
				 FttTraverseFlags flags,
				 FttCellTraverseFunc func,
				 gpointer data)
{
  if (!GFS_IS_MIXED (cell))
    return;
  if (order == FTT_PRE_ORDER &&
      (flags == FTT_TRAVERSE_ALL ||
       ((flags & FTT_TRAVERSE_LEAFS) != 0 && FTT_CELL_IS_LEAF (cell)) ||
       ((flags & FTT_TRAVERSE_NON_LEAFS) != 0 && !FTT_CELL_IS_LEAF (cell))))
    (* func) (cell, data);
  if (!FTT_CELL_IS_LEAF (cell)) {
    struct _FttOct * children = cell->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);

      if (!FTT_CELL_IS_DESTROYED (c))
	cell_traverse_mixed (c, order, flags, func, data);
    }
  }
  if (order == FTT_POST_ORDER &&
      (flags == FTT_TRAVERSE_ALL ||
       ((flags & FTT_TRAVERSE_LEAFS) != 0 && FTT_CELL_IS_LEAF (cell)) ||
       ((flags & FTT_TRAVERSE_NON_LEAFS) != 0 && !FTT_CELL_IS_LEAF (cell))))
    (* func) (cell, data);
}

/**
 * gfs_cell_traverse_mixed:
 * @root: the root #FttCell of the tree to traverse.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 * 
 * Traverses a cell tree starting at the given root #FttCell. Calls
 * the given function for each mixed cell.
 */
void gfs_cell_traverse_mixed (FttCell * root,
			      FttTraverseType order,
			      FttTraverseFlags flags,
			      FttCellTraverseFunc func,
			      gpointer data)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (func != NULL);

  cell_traverse_mixed (root, order, flags, func, data);
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
 * @s: a #GtsSurface.
 *
 * Returns: a new #GtsSurface containing the faces of @s which may
 * intersect @cell or %NULL if no faces of @s intersects @cell.
 */
GtsSurface * gfs_cell_is_cut (FttCell * cell, GtsSurface * s)
{
  GtsSurface * s1 = NULL;
  gpointer data[2];
  GtsBBox bb;
  FttVector p;
  gdouble h;

  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (s != NULL, NULL);

  h = ftt_cell_size (cell)/1.99999;
  ftt_cell_pos (cell, &p);
  bb.x1 = p.x - h; bb.y1 = p.y - h;
  bb.x2 = p.x + h; bb.y2 = p.y + h; 
#if FTT_2D
  bb.z1 = bb.z2 = 0.;
#elif FTT_2D3
  bb.z1 = p.z - 1./1.99999; bb.z2 = p.z + 1./1.99999;
#else  /* 3D */
  bb.z1 = p.z - h; bb.z2 = p.z + h;
#endif /* 3D */
  data[0] = &bb;
  data[1] = &s1;
  gts_surface_foreach_face (s, (GtsFunc) face_overlaps_box, data);
  return s1;
}

static void cell_traverse_cut (FttCell * cell,
			       GtsSurface * s,
			       FttTraverseType order,
			       FttTraverseFlags flags,
			       FttCellTraverseCutFunc func,
			       gpointer data)
{
  GtsSurface * s1 = gfs_cell_is_cut (cell, s);

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
	cell_traverse_cut (c, s1, order, flags, func, data);
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
 * @s: a #GtsSurface.
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
			    GtsSurface * s,
			    FttTraverseType order,
			    FttTraverseFlags flags,
			    FttCellTraverseCutFunc func,
			    gpointer data)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (func != NULL);

  cell_traverse_cut (root, s, order, flags, func, data);
}

/**
 * gfs_interpolate:
 * @cell: a #FttCell containing location @p.
 * @p: the location at which to interpolate.
 * @v: a #GfsVariable.
 *
 * Interpolates the @v variable of @cell, at location @p. Linear
 * interpolation is used and the boundaries of the domain are treated
 * as planes of symmetry for all variables.
 *
 * Returns: the interpolated value of variable @v at location @p.
 */
gdouble gfs_interpolate (FttCell * cell,
			 FttVector p,
			 GfsVariable * v)
{
  FttVector o;
  FttDirection e[FTT_DIMENSION];
  gdouble size;

  g_return_val_if_fail (cell != NULL, 0.);

  ftt_cell_pos (cell, &o);
  size = ftt_cell_size (cell)/2.;
  p.x = (p.x - o.x)/size;
  p.y = (p.y - o.y)/size;
#if FTT_2D
  {
    gdouble f[4], a, b, c, d;
    
    e[0] = FTT_LEFT; e[1] = FTT_BOTTOM;
    f[0] = gfs_cell_corner_value (cell, e, v, -1);
    e[0] = FTT_RIGHT; e[1] = FTT_BOTTOM;
    f[1] = gfs_cell_corner_value (cell, e, v, -1);
    e[0] = FTT_RIGHT; e[1] = FTT_TOP;
    f[2] = gfs_cell_corner_value (cell, e, v, -1);
    e[0] = FTT_LEFT; e[1] = FTT_TOP;
    f[3] = gfs_cell_corner_value (cell, e, v, -1);
    a = f[1] + f[2] - f[0] - f[3];
    b = f[2] + f[3] - f[0] - f[1];
    c = f[0] - f[1] + f[2] - f[3];
    d = f[0] + f[1] + f[2] + f[3];
    return (a*p.x + b*p.y + c*p.x*p.y + d)/4.;
  }
#else  /* 3D */
  {
    gdouble f[8], c[8];
    
    p.z = (p.z - o.z)/size;
    e[0] = FTT_LEFT; e[1] = FTT_BOTTOM; e[2] = FTT_FRONT;
    f[0] = gfs_cell_corner_value (cell, e, v, -1);
    e[0] = FTT_RIGHT; e[1] = FTT_BOTTOM;
    f[1] = gfs_cell_corner_value (cell, e, v, -1);
    e[0] = FTT_RIGHT; e[1] = FTT_TOP;
    f[2] = gfs_cell_corner_value (cell, e, v, -1);
    e[0] = FTT_LEFT; e[1] = FTT_TOP;
    f[3] = gfs_cell_corner_value (cell, e, v, -1);
    
    e[0] = FTT_LEFT; e[1] = FTT_BOTTOM; e[2] = FTT_BACK;
    f[4] = gfs_cell_corner_value (cell, e, v, -1);
    e[0] = FTT_RIGHT; e[1] = FTT_BOTTOM;
    f[5] = gfs_cell_corner_value (cell, e, v, -1);
    e[0] = FTT_RIGHT; e[1] = FTT_TOP;
    f[6] = gfs_cell_corner_value (cell, e, v, -1);
    e[0] = FTT_LEFT; e[1] = FTT_TOP;
    f[7] = gfs_cell_corner_value (cell, e, v, -1);

    c[0] = - f[0] + f[1] + f[2] - f[3] - f[4] + f[5] + f[6] - f[7];
    c[1] = - f[0] - f[1] + f[2] + f[3] - f[4] - f[5] + f[6] + f[7];
    c[2] =   f[0] + f[1] + f[2] + f[3] - f[4] - f[5] - f[6] - f[7];
    c[3] =   f[0] - f[1] + f[2] - f[3] + f[4] - f[5] + f[6] - f[7];
    c[4] = - f[0] + f[1] + f[2] - f[3] + f[4] - f[5] - f[6] + f[7];
    c[5] = - f[0] - f[1] + f[2] + f[3] + f[4] + f[5] - f[6] - f[7];
    c[6] =   f[0] - f[1] + f[2] - f[3] - f[4] + f[5] - f[6] + f[7];
    c[7] =   f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7];

    return (c[0]*p.x + c[1]*p.y + c[2]*p.z + 
	    c[3]*p.x*p.y + c[4]*p.x*p.z + c[5]*p.y*p.z + 
	    c[6]*p.x*p.y*p.z + 
	    c[7])/8.;
  }
#endif /* 3D */
}

/**
 * gfs_center_curvature:
 * @cell: a #FttCell.
 * @c: a component.
 * @v: a #GfsVariable index.
 *
 * The curvature is normalized by the square of the size of the cell.
 *
 * Returns: the value of the @c component of the curvature of variable @v
 * at the center of the cell.  
 */
gdouble gfs_center_curvature (FttCell * cell,
			      FttComponent c,
			      guint v)
{
  FttCellFace f;
  GfsGradient g = { 0., 0. };

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (c < FTT_DIMENSION, 0.);

  if (GFS_IS_MIXED (cell))
    return 0.;

  f.cell = cell;
  for (f.d = 2*c; f.d <= 2*c + 1; f.d++)
    if ((f.neighbor = ftt_cell_neighbor (cell, f.d))) {
      GfsGradient e;

      gfs_face_gradient (&f, &e, v, -1);
      g.a += e.a;
      g.b += e.b;
    }

  return g.b - g.a*GFS_VARIABLE (cell, v);
}

/**
 * gfs_streamline_curvature:
 * @cell: a #FttCell.
 *
 * The curvature is normalized by the size of the cell.
 *
 * Returns: the value of the curvature of the streamline passing
 * through the center of the cell.
 */
gdouble gfs_streamline_curvature (FttCell * cell)
{
  GfsStateVector * s;
  gdouble u2;

  g_return_val_if_fail (cell != NULL, 0.);

  s = GFS_STATE (cell);
#if FTT_2D
  u2 = s->u*s->u + s->v*s->v;
#else  /* FTT_3D */
  u2 = s->u*s->u + s->v*s->v + s->w*s->w;
#endif /* FTT_3D */

  if (u2 == 0.)
    return 0.;
  else {
    FttComponent i;
    gdouble ugu = 0.;

    for (i = 0; i < FTT_DIMENSION; i++) {
      FttComponent j;
      gdouble ugui = 0.;

      for (j = 0; j < FTT_DIMENSION; j++)
	ugui += GFS_VARIABLE (cell, GFS_VELOCITY_INDEX (j))*
	  gfs_center_gradient (cell, j, GFS_VELOCITY_INDEX (i));
      ugu += ugui*ugui;
    }
    return sqrt (ugu)/u2;
  }
}

/**
 * gfs_cell_laplacian:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 *
 * Returns: an evaluation of the Laplacian of @v at the center of
 * @cell normalized by h^2, where h is the cell size.
 */
gdouble gfs_cell_laplacian (FttCell * cell, GfsVariable * v)
{
  FttCellFace f;
  FttCellNeighbors n;
  GfsGradient g = { 0., 0. };
  FttComponent c;
  gdouble v0;

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  if (GFS_IS_MIXED (cell))
    return 0.;

  c = GFS_VELOCITY_COMPONENT (v->i);

  v0 = GFS_VARIABLE (cell, v->i);
  f.cell = cell;
  ftt_cell_neighbors (cell, &n);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++)
    if (n.c[f.d]) {
      GfsGradient e;

      f.neighbor = n.c[f.d];
      gfs_face_gradient (&f, &e, v->i, -1);
      g.a += e.a;
      g.b += e.b;
    }
    else if (f.d/2 == c) {
      g.a += 1.;
      g.b -= v0;
    }

  return g.b - g.a*v0;
}

static FttCell * cell_corner_neighbor (FttCell * cell,
				       FttDirection * d,
				       gint max_level,
				       gboolean * t_junction)
{
  FttCell * neighbor = ftt_cell_neighbor (cell, d[0]);
  if (!neighbor)
    return NULL;
  else {
    guint level = ftt_cell_level (cell);
    if (ftt_cell_level (neighbor) < level) {
      /* neighbor is at a shallower level */
      if (ftt_cell_child_corner (ftt_cell_parent (cell), d) != cell)
	*t_junction = TRUE;
      return neighbor;
    }
    else {
      if (level == max_level || FTT_CELL_IS_LEAF (neighbor))
	/* neighbor is at the same level */
	return neighbor;
      else {
	/* neighbor is at a deeper level */
	FttCell * n;
	guint i;
	FttDirection d1[FTT_DIMENSION];
	d1[0] = FTT_OPPOSITE_DIRECTION (d[0]);
	for (i = 1; i < FTT_DIMENSION; i++)
	  d1[i] = d[i];
	n = ftt_cell_child_corner (neighbor, d1);
	return n ? n : neighbor;
      }
    }
  }
}

static void interpolator_merge (GfsInterpolator * a, GfsInterpolator * b)
{
  guint i;
  for (i = 0; i < b->n; i++) {
    FttCell * c = b->c[i];
    guint j;

    for (j = 0; j < a->n && c != a->c[j]; j++)
      ;
    if (j < a->n)
      a->w[j] += b->w[i];
    else {
#if FTT_2D
      g_assert (j < 7);
#else
      g_assert (j < 29);
#endif
      a->c[j] = c;
      a->w[j] = b->w[i];
      a->n++;
    }
  }
}

static void interpolator_scale (GfsInterpolator * a, gdouble b)
{
  guint i;
  for (i = 0; i < a->n; i++)
    a->w[i] *= b;
}

static void t_junction_interpolator (FttCell * cell,
				     FttDirection * d,
				     FttCell * n,
				     gint max_level,
				     gboolean centered,
				     GfsInterpolator * inter)
{
  FttDirection d1[FTT_DIMENSION];
  GfsInterpolator a;

  d1[0] = FTT_OPPOSITE_DIRECTION (d[0]);
#if FTT_2D
  d1[1] = d[1];
  gfs_cell_corner_interpolator (n, d1, max_level, centered, inter);
  d1[1] = FTT_OPPOSITE_DIRECTION (d[1]);
  gfs_cell_corner_interpolator (n, d1, max_level, centered, &a);
  interpolator_merge (inter, &a);
  interpolator_scale (inter, 0.5);
#else /* 3D */
  d1[1] = d[1]; d1[2] = d[2];
  gfs_cell_corner_interpolator (n, d1, max_level, centered, inter);
  if (ftt_cell_neighbor_is_brother (cell, d[1])) {
    d1[1] = FTT_OPPOSITE_DIRECTION (d[1]);
    gfs_cell_corner_interpolator (n, d1, max_level, centered, &a);
    interpolator_merge (inter, &a);
    if (ftt_cell_neighbor_is_brother (cell, d[2])) {
      d1[2] = FTT_OPPOSITE_DIRECTION (d[2]);
      gfs_cell_corner_interpolator (n, d1, max_level, centered, &a);
      interpolator_merge (inter, &a);
      d1[1] = d[1];
      gfs_cell_corner_interpolator (n, d1, max_level, centered, &a);
      interpolator_merge (inter, &a);
      interpolator_scale (inter, 0.25);
    }
    else
      interpolator_scale (inter, 0.5);
  }
  else {
    d1[2] = FTT_OPPOSITE_DIRECTION (d[2]);
    gfs_cell_corner_interpolator (n, d1, max_level, centered, &a);
    interpolator_merge (inter, &a);
    interpolator_scale (inter, 0.5);
  }
#endif /* 3D */
}

static gboolean do_path (FttCell * cell, gint i,
			 FttCell * n[N_CELLS],
			 FttDirection * d,
			 gint max_level,
			 gboolean centered,
			 GfsInterpolator * inter)
{
  /* paths from each cell to neighbors. Cell indices are as in
     doc/figures/indices.fig
     first index: cell index
     second index: path index
     third index: < FTT_DIMENSION: directions.
                  = FTT_DIMENSION: neighboring cell index */
  static gint path[N_CELLS][FTT_DIMENSION][FTT_DIMENSION + 1] = {
#if FTT_2D
    {{1,2,1},   {2,1,2}},
    {{2,-1,3},  {-1,2,0}},
    {{1,-2,3},  {-2,1,0}},
    {{-1,-2,2}, {-2,-1,1}}
#else /* 3D */
    {{1,2,3,1},    {2,1,3,2},    {3,1,2,4}},
    {{2,3,-1,3},   {3,2,-1,5},   {-1,2,3,0}},
    {{1,-2,3,3},   {3,-2,1,6},   {-2,3,1,0}},
    {{3,-1,-2,7},  {-2,-1,3,1},  {-1,-2,3,2}},
    {{2,1,-3,6},   {1,2,-3,5},   {-3,1,2,0}},
    {{2,-1,-3,7},  {-1,2,-3,4},  {-3,2,-1,1}},
    {{1,-2,-3,7},  {-2,1,-3,4},  {-3,1,-2,2}},
    {{-1,-2,-3,6}, {-2,-1,-3,5}, {-3,-1,-2,3}}
#endif /* 3D */
  };
  guint j;

  for (j = 0; j < FTT_DIMENSION; j++) {
    guint k = path[i][j][FTT_DIMENSION];

    if (n[k] == NULL) {
      gboolean t_junction = FALSE;
      FttDirection d1[FTT_DIMENSION];
      guint l;

      for (l = 0; l < FTT_DIMENSION; l++)
	d1[l] = path[i][j][l] < 0 ? FTT_OPPOSITE_DIRECTION (d[- path[i][j][l] - 1]) : 
	  d[path[i][j][l] - 1];
      n[k] = cell_corner_neighbor (cell, d1, max_level, &t_junction);
      if (t_junction) {
	t_junction_interpolator (cell, d1, n[k], max_level, centered, inter);
	return TRUE;
      }
      if (n[k]) {
	t_junction = do_path (n[k], k, n, d, max_level, centered, inter);
	if (t_junction)
	  return TRUE;
      }
    }
  }
  return FALSE;
}

static gdouble distance (FttVector * c, FttCell * cell, gboolean centered)
{
  if (centered || !GFS_IS_MIXED (cell))
    return ftt_cell_size (cell)*
#if FTT_2D
      0.707106781185
#else  /* 3D */
      0.866025403785
#endif /* 3D */
      ;
  else {
    FttVector cm;
    gfs_cell_cm (cell, &cm);
    /* fixme: what about periodic boundaries? */
    return sqrt ((cm.x - c->x)*(cm.x - c->x) + (cm.y - c->y)*(cm.y - c->y)
#if (!FTT_2D)
      + (cm.z - c->z)*(cm.z - c->z)
#endif /* 3D */
		 );
  }
}

/**
 * gfs_cell_corner_interpolator:
 * @cell: a #FttCell.
 * @d: a set of perpendicular directions.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 * @centered: %TRUE if the interpolator is cell-centered.
 * @inter: a #GfsInterpolator.
 *
 * Fills @inter with the interpolator for the corner of @cell defined by @d.
 */
void gfs_cell_corner_interpolator (FttCell * cell,
				   FttDirection d[FTT_DIMENSION],
				   gint max_level,
				   gboolean centered,
				   GfsInterpolator * inter)
{
  FttCell * n[N_CELLS];
  guint i;
  gboolean t_junction;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (inter != NULL);

  while (!FTT_CELL_IS_LEAF (cell) && 
	 ftt_cell_level (cell) != max_level &&
	 (n[0] = ftt_cell_child_corner (cell, d)))
    cell = n[0];
  n[0] = cell;
  for (i = 1; i < N_CELLS; i++)
    n[i] = NULL;
  t_junction = do_path (cell, 0, n, d, max_level, centered, inter);
  if (t_junction)
    return;

  {
    FttVector c;
    gdouble w = 0.;

    inter->n = 0;
    ftt_corner_pos (cell, d, &c);
    for (i = 0; i < N_CELLS; i++)
      if (n[i]) {
	gdouble a;
	a = 1./(distance (&c, n[i], centered) + 1e-12);
	inter->c[inter->n] = n[i];
	inter->w[inter->n++] = a;
	w += a;
      }
    g_assert (w > 0.);
    interpolator_scale (inter, 1./w);
  }
}

/**
 * gfs_cell_corner_value:
 * @cell: a #FttCell.
 * @d: a set of perpendicular directions.
 * @v: a #GfsVariable.
 * @max_level: the maximum cell level to consider (-1 means no restriction).
 *
 * Returns: the value of variable @v interpolated at the corner of
 * @cell defined by @d.
 */
gdouble gfs_cell_corner_value (FttCell * cell,
			       FttDirection d[FTT_DIMENSION],
			       GfsVariable * v,
			       gint max_level)
{
  GfsInterpolator inter;
  gdouble val = 0.;
  guint i;

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  gfs_cell_corner_interpolator (cell, d, max_level, v->centered, &inter);
  for (i = 0; i < inter.n; i++)
    val += inter.w[i]*GFS_VARIABLE (inter.c[i], v->i);
  return val;
}
