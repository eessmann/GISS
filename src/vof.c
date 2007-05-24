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
#include "vof.h"
#include "variable.h"
#include "adaptive.h"
#include "graphic.h"

#define THRESHOLD(c) {if ((c) < 0.) c = 0.; else if ((c) > 1.) c = 1.;}

/**
 * gfs_line_area:
 * @m: normal to the line.
 * @alpha: line constant.
 *
 * Returns: the area of the fraction of a cell lying under the line
 * (@m,@alpha).
 */
gdouble gfs_line_area (FttVector * m, gdouble alpha)
{
  FttVector n;
  gdouble alpha1, a, v;

  g_return_val_if_fail (m != NULL, 0.);

  n = *m;
  alpha1 = alpha;
  if (n.x < 0.) {
    alpha1 -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha1 -= n.y;
    n.y = - n.y;
  }

  if (alpha1 <= 0.)
    return 0.;

  if (alpha1 >= n.x + n.y)
    return 1.;

  if (n.x == 0.)
    return alpha1/n.y;
  else if (n.y == 0.)
    return alpha1/n.x;
  else {
    v = alpha1*alpha1;

    a = alpha1 - n.x;
    if (a > 0.)
      v -= a*a;
    
    a = alpha1 - n.y;
    if (a > 0.)
      v -= a*a;

    return v/(2.*n.x*n.y);
  }
}

/**
 * gfs_line_alpha:
 * @m: a #FttVector.
 * @c: a volume fraction.
 *
 * Returns: the value @alpha such that the area of a square cell
 * lying under the line defined by @m.@x = @alpha is equal to @c. 
 */
gdouble gfs_line_alpha (FttVector * m, gdouble c)
{
  gdouble alpha, m1, m2, v1;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (c >= 0. && c <= 1., 0.);
  
  m1 = fabs (m->x); m2 = fabs (m->y);
  if (m1 > m2) {
    v1 = m1; m1 = m2; m2 = v1;
  }
  
  v1 = m1/2.;
  if (c <= v1/m2)
    alpha = sqrt (2.*c*m1*m2);
  else if (c <= 1. - v1/m2)
    alpha = c*m2 + v1;
  else
    alpha = m1 + m2 - sqrt (2.*m1*m2*(1. - c));

  if (m->x < 0.)
    alpha += m->x;
  if (m->y < 0.)
    alpha += m->y;

  return alpha;
}

/**
 * gfs_line_center:
 * @m: normal to the line.
 * @alpha: line constant.
 * @a: area of cell fraction.
 * @p: a #FttVector.
 *
 * Fills @p with the position of the center of mass of the fraction of
 * a square cell lying under the line (@m,@alpha).
 */
void gfs_line_center (FttVector * m, gdouble alpha, gdouble a, FttVector * p)
{
  FttVector n;
  gdouble b;

  g_return_if_fail (m != NULL);
  g_return_if_fail (p != NULL);

  n = *m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }

  p->z = 0.;
  if (alpha <= 0.) {
    p->x = p->y = 0.;
    return;
  }

  if (alpha >= n.x + n.y) {
    p->x = p->y = 0.5;
    return;
  }

  g_return_if_fail (a > 0. && a < 1.);

  n.x += 1e-4; n.y += 1e-4;

  p->x = p->y = alpha*alpha*alpha;

  b = alpha - n.x;
  if (b > 0.) {
    p->x -= b*b*(alpha + 2.*n.x);
    p->y -= b*b*b;
  }

  b = alpha - n.y;
  if (b > 0.) {
    p->y -= b*b*(alpha + 2.*n.y);
    p->x -= b*b*b;
  }
  
  p->x /= 6.*n.x*n.x*n.y*a;
  p->y /= 6.*n.x*n.y*n.y*a;

  if (m->x < 0.)
    p->x = 1. - p->x;
  if (m->y < 0.)
    p->y = 1. - p->y;
}

/**
 * gfs_line_area_center:
 * @m: normal to the line.
 * @alpha: line constant.
 * @p: a #FttVector.
 *
 * Fills @p with the position of the center of area of the fraction of
 * a square cell lying under the line (@m,@alpha).
 *
 * Returns: the length of the facet.
 */
gdouble gfs_line_area_center (FttVector * m, gdouble alpha, FttVector * p)
{
  FttVector n;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (p != NULL, 0.);

  n = *m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }

  p->z = 0.;
  if (alpha <= 0. || alpha >= n.x + n.y) {
    p->x = p->y = 0.;
    return 0.;
  }

  n.x += 1e-4; n.y += 1e-4;

  p->x = p->y = 0.;

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  gdouble ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

  p->x /= 2.;
  p->y /= 2.;

  THRESHOLD (p->x);
  THRESHOLD (p->y);

  if (m->x < 0.)
    p->x = 1. - p->x;
  if (m->y < 0.)
    p->y = 1. - p->y;

  return sqrt (ax*ax + ay*ay);
}

#if (!FTT_2D)
/**
 * gfs_plane_volume:
 * @m: normal to the plane.
 * @alpha: plane constant.
 *
 * Returns: the volume of a cell lying under the plane (@m,@alpha).
 */
gdouble gfs_plane_volume (FttVector * m, gdouble alpha)
{
  g_return_val_if_fail (m != NULL, 0.);

  gdouble al = alpha + MAX(0., -m->x) + MAX(0., -m->y) + MAX(0., -m->z);
  if (al <= 0.)
    return 0.;
  gdouble tmp = fabs(m->x) + fabs(m->y) + fabs(m->z);
  if (al >= tmp)
    return 1.;
  g_assert (tmp > 0.);
  gdouble n1 = fabs(m->x)/tmp;
  gdouble n2 = fabs(m->y)/tmp;
  gdouble n3 = fabs(m->z)/tmp;
  al = MAX(0., MIN(1., al/tmp));
  gdouble al0 = MIN(al, 1. - al);
  gdouble b1 = MIN(n1*1, n2);
  gdouble b3 = MAX(n1*1, n2);
  gdouble b2 = n3;
  if (b2 < b1) {
    tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  else if (b2 > b3) {
    tmp = b3;
    b3 = b2;
    b2 = tmp;
  }
  gdouble b12 = b1 + b2;
  gdouble bm = MIN(b12, b3);
  gdouble pr = MAX(6.*b1*b2*b3, 1e-50);
  if (al0 < b1)
    tmp = al0*al0*al0/pr;
  else if (al0 < b2)
    tmp = 0.5*al0*(al0 - b1)/(b2*b3) +  b1*b1*b1/pr;
  else if (al0 < bm)
    tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) + b2*b2*(b2 - 3.*al0))/pr;
  else if (b12 < b3)
    tmp = (al0 - 0.5*bm)/b3;
  else
    tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) + 
	   b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;

  return al <= 0.5 ? tmp : 1. - tmp;
}

/**
 * gfs_plane_alpha:
 * @m: a #FttVector.
 * @c: a volume fraction.
 *
 * Returns: the value @alpha such that the volume of a cubic cell
 * lying under the plane defined by @m.@x = @alpha is equal to @c. 
 */
gdouble gfs_plane_alpha (FttVector * m, gdouble c)
{
  gdouble alpha;
  FttVector n;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (c >= 0. && c <= 1., 0.);

  n.x = fabs (m->x); n.y = fabs (m->y); n.z = fabs (m->z);

  gdouble m1, m2, m3;
  m1 = MIN(n.x, n.y);
  m3 = MAX(n.x, n.y);
  m2 = n.z;
  if (m2 < m1) {
    gdouble tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  else if (m2 > m3) {
    gdouble tmp = m3;
    m3 = m2;
    m2 = tmp;
  }
  gdouble m12 = m1 + m2;
  gdouble pr = MAX(6.*m1*m2*m3, 1e-50);
  gdouble V1 = m1*m1*m1/pr;
  gdouble V2 = V1 + (m2 - m1)/(2.*m3), V3;
  gdouble mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  gdouble ch = MIN(c, 1. - c);
  if (ch < V1)
    alpha = pow (pr*ch, 1./3.);
  else if (ch < V2)
    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    gdouble p = 2.*m1*m2;
    gdouble q = 3.*m1*m2*(m12 - 2.*m3*ch)/2.;
    gdouble p12 = sqrt (p);
    gdouble teta = acos(q/(p*p12))/3.;
    gdouble cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 < m3)
    alpha = m3*ch + mm/2.;
  else {
    gdouble p = m1*(m2 + m3) + m2*m3 - 1./4.;
    gdouble q = 3.*m1*m2*m3*(1./2. - ch)/2.;
    gdouble p12 = sqrt(p);
    gdouble teta = acos(q/(p*p12))/3.;
    gdouble cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (m->x < 0.)
    alpha += m->x;
  if (m->y < 0.)
    alpha += m->y;
  if (m->z < 0.)
    alpha += m->z;

  return alpha;
}

/**
 * gfs_plane_center:
 * @m: normal to the plane.
 * @alpha: plane constant.
 * @a: volume of cell fraction.
 * @p: a #FttVector.
 *
 * Fills @p with the position of the center of mass of the fraction of
 * a cubic cell lying under the plane (@m,@alpha).
 */
void gfs_plane_center (FttVector * m, gdouble alpha, gdouble a, FttVector * p)
{
  FttVector n;
  gdouble b, amax;

  g_return_if_fail (m != NULL);
  g_return_if_fail (p != NULL);

  n = *m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }  
  if (n.z < 0.) {
    alpha -= n.z;
    n.z = - n.z;
  }  

  if (alpha <= 0.) {
    p->x = p->y = p->z = 0.;
    return;
  }

  if (alpha >= n.x + n.y + n.z) {
    p->x = p->y = p->z = 0.5;
    return;
  }

  g_return_if_fail (a > 0. && a < 1.);

  n.x += 1e-4; n.y += 1e-4; n.z += 1e-4;

  amax = n.x + n.y + n.z;
  p->x = p->y = p->z = alpha*alpha*alpha*alpha;

  b = alpha - n.x;
  if (b > 0.) {
    p->x -= b*b*b*(3.*n.x + alpha);
    p->y -= b*b*b*b;
    p->z -= b*b*b*b;
  }
  b = alpha - n.y;
  if (b > 0.) {
    p->y -= b*b*b*(3.*n.y + alpha);
    p->x -= b*b*b*b;
    p->z -= b*b*b*b;
  }
  b = alpha - n.z;
  if (b > 0.) {
    p->z -= b*b*b*(3.*n.z + alpha);
    p->x -= b*b*b*b;
    p->y -= b*b*b*b;
  }

  amax = alpha - amax;
  b = amax + n.x;
  if (b > 0.) {
    p->y += b*b*b*(3.*n.y + alpha - n.z);
    p->z += b*b*b*(3.*n.z + alpha - n.y);
    p->x += b*b*b*b;
  }
  b = amax + n.y;
  if (b > 0.) {
    p->x += b*b*b*(3.*n.x + alpha - n.z);
    p->z += b*b*b*(3.*n.z + alpha - n.x);
    p->y += b*b*b*b;
  }
  b = amax + n.z;
  if (b > 0.) {
    p->x += b*b*b*(3.*n.x + alpha - n.y);
    p->y += b*b*b*(3.*n.y + alpha - n.x);
    p->z += b*b*b*b;
  }

  b = 24.*n.x*n.y*n.z*a;
  p->x /= b*n.x; p->y /= b*n.y; p->z /= b*n.z;

  if (m->x < 0.) p->x = 1. - p->x;
  if (m->y < 0.) p->y = 1. - p->y;
  if (m->z < 0.) p->z = 1. - p->z;
}

/**
 * gfs_plane_area_center:
 * @m: normal to the plane.
 * @alpha: plane constant.
 * @p: a #FttVector.
 *
 * Fills @p with the position of the center of area of the fraction of
 * a cubic cell lying under the plane (@m,@alpha).
 *
 * Returns: the area of the facet.
 */
gdouble gfs_plane_area_center (FttVector * m, gdouble alpha, FttVector * p)
{
  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (p != NULL, 0.);

  if (fabs (m->x) < 1e-4) {
    FttVector n, q;
    n.x = m->y;
    n.y = m->z;
    gdouble area = gfs_line_area_center (&n, alpha, &q);
    p->x = 0.5;
    p->y = q.x;
    p->z = q.y;
    return area;
  }
  if (fabs (m->y) < 1e-4) {
    FttVector n, q;
    n.x = m->z;
    n.y = m->x;
    gdouble area = gfs_line_area_center (&n, alpha, &q);
    p->x = q.y;
    p->y = 0.5;
    p->z = q.x;
    return area;
  }
  if (fabs (m->z) < 1e-4) {
    gdouble area = gfs_line_area_center (m, alpha, p);
    p->z = 0.5;
    return area;
  }

  FttVector n = *m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }  
  if (n.z < 0.) {
    alpha -= n.z;
    n.z = - n.z;
  }

  gdouble amax = n.x + n.y + n.z;
  if (alpha <= 0. || alpha >= amax) {
    p->x = p->y = p->z = 0.;
    return 0.;
  }

  gdouble area = alpha*alpha;
  p->x = p->y = p->z = area*alpha;

  gdouble b = alpha - n.x;
  if (b > 0.) {
    area -= b*b;
    p->x -= b*b*(2.*n.x + alpha);
    p->y -= b*b*b;
    p->z -= b*b*b;
  }
  b = alpha - n.y;
  if (b > 0.) {
    area -= b*b;
    p->y -= b*b*(2.*n.y + alpha);
    p->x -= b*b*b;
    p->z -= b*b*b;
  }
  b = alpha - n.z;
  if (b > 0.) {
    area -= b*b;
    p->z -= b*b*(2.*n.z + alpha);
    p->x -= b*b*b;
    p->y -= b*b*b;
  }

  amax = alpha - amax;
  b = amax + n.x;
  if (b > 0.) {
    area += b*b;
    p->y += b*b*(2.*n.y + alpha - n.z);
    p->z += b*b*(2.*n.z + alpha - n.y);
    p->x += b*b*b;
  }
  b = amax + n.y;
  if (b > 0.) {
    area += b*b;
    p->x += b*b*(2.*n.x + alpha - n.z);
    p->z += b*b*(2.*n.z + alpha - n.x);
    p->y += b*b*b;
  }
  b = amax + n.z;
  if (b > 0.) {
    area += b*b;
    p->x += b*b*(2.*n.x + alpha - n.y);
    p->y += b*b*(2.*n.y + alpha - n.x);
    p->z += b*b*b;
  }

  area *= 3.;
  p->x /= area*n.x;
  p->y /= area*n.y;
  p->z /= area*n.z;

  THRESHOLD (p->x);
  THRESHOLD (p->y);
  THRESHOLD (p->z);

  if (m->x < 0.) p->x = 1. - p->x;
  if (m->y < 0.) p->y = 1. - p->y;
  if (m->z < 0.) p->z = 1. - p->z;

  return area*sqrt (1./(n.x*n.x*n.y*n.y) + 1./(n.x*n.x*n.z*n.z) + 1./(n.z*n.z*n.y*n.y))/6.;
}
#endif /* 3D */

/**
 * gfs_youngs_gradient:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 * @g: a #FttVector.
 *
 * Fills @g with the Youngs-averaged gradients of @v 
 * normalised by the size of @cell.
 */
void gfs_youngs_gradient (FttCell * cell, GfsVariable * v, FttVector * g)
{
  static FttDirection d[(FTT_DIMENSION - 1)*4][FTT_DIMENSION] = {
#if FTT_2D
    {FTT_RIGHT, FTT_TOP}, {FTT_LEFT, FTT_TOP}, {FTT_LEFT, FTT_BOTTOM}, {FTT_RIGHT, FTT_BOTTOM}
#else  /* 3D */
    {FTT_RIGHT, FTT_TOP, FTT_FRONT}, {FTT_LEFT, FTT_TOP, FTT_FRONT}, 
    {FTT_LEFT, FTT_BOTTOM, FTT_FRONT}, {FTT_RIGHT, FTT_BOTTOM, FTT_FRONT},
    {FTT_RIGHT, FTT_TOP, FTT_BACK}, {FTT_LEFT, FTT_TOP, FTT_BACK}, 
    {FTT_LEFT, FTT_BOTTOM, FTT_BACK}, {FTT_RIGHT, FTT_BOTTOM, FTT_BACK},
#endif /* 3D */
  };
  gdouble u[(FTT_DIMENSION - 1)*4];
  guint i;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);
  g_return_if_fail (g != NULL);

  for (i = 0; i < (FTT_DIMENSION - 1)*4; i++)
    u[i] = gfs_cell_corner_value (cell, d[i], v, -1);

#if FTT_2D
  g->x = (u[0] + u[3] - u[1] - u[2])/2.;
  g->y = (u[0] + u[1] - u[2] - u[3])/2.;
#else  /* 3D */
  g->x = (u[0] + u[3] + u[4] + u[7] - u[1] - u[2] - u[5] - u[6])/4.;
  g->y = (u[0] + u[1] + u[4] + u[5] - u[2] - u[3] - u[6] - u[7])/4.;
  g->z = (u[0] + u[1] + u[2] + u[3] - u[4] - u[5] - u[6] - u[7])/4.;
#endif /* 3D */
}

/**
 * gfs_vof_plane_interpolate:
 * @cell: a #FttCell containing location @p.
 * @p: the center of the virtual cell.
 * @level: the level of the virtual cell.
 * @t: a #GfsVariableTracerVOF.
 * @m: a #FttVector.
 *
 * Computes the equation @m.x = alpha of the volume fraction plane of
 * a virtual cell at @level centered on @p.
 *
 * Returns: alpha for the virtual cell.
 */
gdouble gfs_vof_plane_interpolate (FttCell * cell,
				   FttVector * p,
				   guint level,
				   GfsVariableTracerVOF * t,
				   FttVector * m)
{
  guint l = ftt_cell_level (cell);

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (l <= level, 0.);
  g_return_val_if_fail (t != NULL, 0.);
  g_return_val_if_fail (m != NULL, 0.);

  GfsVariable * v = GFS_VARIABLE1 (t);
  gdouble f = GFS_VARIABLE (cell, v->i);
  g_return_val_if_fail (!GFS_IS_FULL (f), 0.);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&m->x)[c] = GFS_VARIABLE (cell, t->m[c]->i);

  gdouble alpha = GFS_VARIABLE (cell, t->alpha->i);
  if (l < level) {
    gdouble h = ftt_level_size (level);
    gdouble H = ftt_cell_size (cell);
    FttVector q;
    
    ftt_cell_pos (cell, &q);
    alpha *= H;
    for (c = 0; c < FTT_DIMENSION; c++)
      alpha -= (&m->x)[c]*((&p->x)[c] - h/2. - (&q.x)[c] + H/2);
    alpha /= h;
  }
  return alpha;
}

/**
 * gfs_vof_interpolate:
 * @cell: a #FttCell containing location @p.
 * @p: the center of the virtual cell.
 * @level: the level of the virtual cell.
 * @t: a #GfsVariableTracerVOF.
 *
 * Computes the volume fraction of a virtual cell at @level centered
 * on @p.
 *
 * Returns: the volume fraction of the virtual cell.
 */
gdouble gfs_vof_interpolate (FttCell * cell,
			     FttVector * p,
			     guint level,
			     GfsVariableTracerVOF * t)
{
  guint l = ftt_cell_level (cell);

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (l <= level, 0.);
  g_return_val_if_fail (t != NULL, 0.);

  GfsVariable * v = GFS_VARIABLE1 (t);
  gdouble f = GFS_VARIABLE (cell, v->i);
  if (l == level || GFS_IS_FULL (f))
    return f;
  else {
    FttVector m;
    gdouble alpha = gfs_vof_plane_interpolate (cell, p, level, t, &m);
    return gfs_plane_volume (&m, alpha);
  }
}

/* GfsVariableTracerVOF: object */

static FttCell * domain_and_boundary_locate (GfsDomain * domain, FttVector p, guint level)
{
  FttCell * cell = gfs_domain_locate (domain, p, level);
  if (cell)
    return cell;
  return gfs_domain_boundary_locate (domain, p, level);
}

#if FTT_2D
# define F(x,y,z) f[x][y]
#else
# define F(x,y,z) f[x][y][z]
#endif

static void stencil (FttCell * cell, GfsVariable * v, gdouble F(3,3,3))
{
  gdouble h = ftt_cell_size (cell);
  guint level = ftt_cell_level (cell);
  FttVector p;
  gint x, y, z = 0;
  
  F(1,1,1) = GFS_VARIABLE (cell, v->i);
  ftt_cell_pos (cell, &p);
#if !FTT_2D
  for (z = -1; z <= 1; z++)
#endif
    for (x = -1; x <= 1; x++)
      for (y = -1; y <= 1; y++)
	if (x != 0 || y != 0 || z != 0) {
	  FttVector o;
	  o.x = p.x + h*x; o.y = p.y + h*y; o.z = p.z + h*z;
	  FttCell * neighbor = domain_and_boundary_locate (v->domain, o, level);
	  if (neighbor)
	    F(x + 1, y + 1, z + 1) =
	      gfs_vof_interpolate (neighbor, &o, level, GFS_VARIABLE_TRACER_VOF (v));
	  else
	    F(x + 1, y + 1, z + 1) = -1.;
	}
  /* boundary conditions (symmetry) */
#if FTT_2D
  for (x = 0; x <= 2; x++) {
    if (f[x][0] < 0.) f[x][0] = f[x][1];
    if (f[x][2] < 0.) f[x][2] = f[x][1];
  }
  for (y = 0; y <= 2; y++) {
    if (f[0][y] < 0.) f[0][y] = f[1][y];
    if (f[2][y] < 0.) f[2][y] = f[1][y];
  }
#else /* 3D */
  for (x = 0; x <= 2; x++)
    for (y = 0; y <= 2; y++) {
      if (f[x][y][0] < 0.) f[x][y][0] = f[x][y][1];
      if (f[x][y][2] < 0.) f[x][y][2] = f[x][y][1];
    }
  for (x = 0; x <= 2; x++)
    for (z = 0; z <= 2; z++) {
      if (f[x][0][z] < 0.) f[x][0][z] = f[x][1][z];
      if (f[x][2][z] < 0.) f[x][2][z] = f[x][1][z];
    }
  for (z = 0; z <= 2; z++)
    for (y = 0; y <= 2; y++) {
      if (f[0][y][z] < 0.) f[0][y][z] = f[1][y][z];
      if (f[2][y][z] < 0.) f[2][y][z] = f[1][y][z];
    }
#endif /* 3D */
}

static void youngs_normal (FttCell * cell, GfsVariable * v, FttVector * n)
{
  gdouble F(3,3,3);

  stencil (cell, v, f);
#if FTT_2D
  n->x = (f[0][2] + 2.*f[0][1] + f[0][0] - 2.*f[2][1] - f[2][2] - f[2][0])/8.;
  n->y = (f[2][0] + 2.*f[1][0] + f[0][0] - 2.*f[1][2] - f[2][2] - f[0][2])/8.;
  n->z = 0.;
#else  /* 3D */
  gdouble mm1 = f[0][0][0] + f[0][0][2] + f[0][2][0] + f[0][2][2] +
    2.*(f[0][0][1] + f[0][2][1] + f[0][1][0] + f[0][1][2]) + 
    4.*f[0][1][1];
  gdouble mm2 = f[2][0][0] + f[2][0][2] + f[2][2][0] + f[2][2][2] + 
    2.*(f[2][0][1] + f[2][2][1] + f[2][1][0] + f[2][1][2]) + 
    4.*f[2][1][1];
  n->x = (mm1 - mm2)/32.;
    
  mm1 = f[0][0][0] + f[0][0][2] + f[2][0][0] + f[2][0][2] + 
    2.*(f[0][0][1] + f[2][0][1] + f[1][0][0] + f[1][0][2]) + 
    4.*f[1][0][1];
  mm2 = f[0][2][0] + f[0][2][2] + f[2][2][0] + f[2][2][2] + 
    2.*(f[0][2][1] + f[2][2][1] + f[1][2][0] + f[1][2][2]) + 
    4.*f[1][2][1];
  n->y = (mm1 - mm2)/32.;
                  
  mm1 = f[0][0][0] + f[0][2][0] + f[2][0][0] + f[2][2][0] +
    2.*(f[0][1][0] + f[2][1][0] + f[1][0][0] + f[1][2][0]) + 
    4.*f[1][1][0];
  mm2 = f[0][0][2] + f[0][2][2] + f[2][0][2] + f[2][2][2] + 
    2.*(f[0][1][2] + f[2][1][2] + f[1][0][2] + f[1][2][2]) + 
    4.*f[1][1][2];
  n->z = (mm1 - mm2)/32.;
#endif /* 3D */
}

#if FTT_2D
# include "myc2d.h"
#else
# include "myc.h"
#endif

static void myc_normal (FttCell * cell, GfsVariable * v, FttVector * n)
{
  gdouble F(3,3,3);  

  stencil (cell, v, f);
  mycs (f, &n->x);
#if FTT_2D
  n->z = 0.;
#endif
}

static void vof_plane (FttCell * cell, GfsVariable * v)
{
  if (FTT_CELL_IS_LEAF (cell)) {
    GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
    gdouble f = GFS_VARIABLE (cell, v->i);
    FttComponent c;

    THRESHOLD (f);
    if (GFS_IS_FULL (f)) {
      for (c = 1; c < FTT_DIMENSION; c++)
	GFS_VARIABLE (cell, t->m[c]->i) = 0.;
      GFS_VARIABLE (cell, t->m[0]->i) = 1.;
      GFS_VARIABLE (cell, t->alpha->i) = f;
    }
    else {
      FttVector m;
      gdouble n = 0.;

      myc_normal (cell, v, &m);
      for (c = 0; c < FTT_DIMENSION; c++)
	n += fabs ((&m.x)[c]);
      if (n > 0.)
	for (c = 0; c < FTT_DIMENSION; c++)
	  (&m.x)[c] /= n;
      else /* fixme: this is a small fragment */
	m.x = 1.;
      for (c = 0; c < FTT_DIMENSION; c++)
	GFS_VARIABLE (cell, t->m[c]->i) = (&m.x)[c];
      GFS_VARIABLE (cell, t->alpha->i) = gfs_plane_alpha (&m, f);
    }
  }
}

static void variable_tracer_vof_update (GfsVariable * v, GfsDomain * domain)
{
  GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
  guint l, depth = gfs_domain_depth (domain);
  FttComponent c;
  for (l = 0; l <= depth; l++) {
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, l,
			      (FttCellTraverseFunc) vof_plane, v);
    for (c = 0; c < FTT_DIMENSION; c++)
      gfs_domain_bc (domain, FTT_TRAVERSE_LEVEL, l, t->m[c]);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEVEL, l, t->alpha);
  }
}

static gboolean variable_tracer_vof_event (GfsEvent * event, 
					   GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_tracer_vof_class ())->parent_class)->event)
      (event, sim)) {
    GfsVariable * v = GFS_VARIABLE1 (event);
    GfsDomain * domain = GFS_DOMAIN (sim);
    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_get_from_below_intensive, v);
    gfs_domain_bc (GFS_DOMAIN (sim), FTT_TRAVERSE_ALL, -1, v);

    GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
    if (!t->alpha) {
      FttComponent c;
      for (c = 0; c < FTT_DIMENSION; c++)
	t->m[c] = gfs_temporary_variable (domain);
      t->alpha = gfs_temporary_variable (domain);
    }
    variable_tracer_vof_update (v, domain);
    return TRUE;
  }
  return FALSE;
}

static void variable_tracer_vof_destroy (GtsObject * o)
{
  GfsVariableTracerVOF * v = GFS_VARIABLE_TRACER_VOF (o);

  if (v->alpha) {
    FttComponent c;
    for (c = 0; c < FTT_DIMENSION; c++)
      gts_object_destroy (GTS_OBJECT (v->m[c]));
    gts_object_destroy (GTS_OBJECT (v->alpha));
  }

  (* GTS_OBJECT_CLASS (gfs_variable_tracer_vof_class ())->parent_class->destroy) (o);
}

static void variable_tracer_vof_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_tracer_vof_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (GFS_VARIABLE_TRACER (*o)->advection.cfl > 0.5) {
    gts_file_error (fp, "cfl `%g' is out of range `]0,0.5]'", 
		    GFS_VARIABLE_TRACER (*o)->advection.cfl);
    return;
  }  

  GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (*o);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    t->m[c] = gfs_temporary_variable (GFS_VARIABLE1 (*o)->domain);
  t->alpha = gfs_temporary_variable (GFS_VARIABLE1 (*o)->domain);
}

static void variable_tracer_vof_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = variable_tracer_vof_event;
  klass->destroy = variable_tracer_vof_destroy;
  klass->read = variable_tracer_vof_read;
}

static void vof_coarse_fine (FttCell * parent, GfsVariable * v)
{
  GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
  gdouble f = GFS_VARIABLE (parent, v->i);
  FttCellChildren child;
  guint i;
  
  ftt_cell_children (parent, &child);
  if (GFS_IS_FULL (f))
    for (i = 0; i < FTT_CELLS; i++) {
      FttComponent c;
      if (!child.c[i])
	g_assert_not_implemented ();
      GFS_VARIABLE (child.c[i], v->i) = f;
      for (c = 1; c < FTT_DIMENSION; c++)
	GFS_VARIABLE (child.c[i], t->m[c]->i) = 0.;
      GFS_VARIABLE (child.c[i], t->m[0]->i) = 1.;
      GFS_VARIABLE (child.c[i], t->alpha->i) = f;
    }
  else {
    gdouble alpha = GFS_VARIABLE (parent, t->alpha->i);
    FttVector m;
    
    for (i = 0; i < FTT_DIMENSION; i++)
      (&m.x)[i] = GFS_VARIABLE (parent, t->m[i]->i);
    for (i = 0; i < FTT_CELLS; i++) {
      gdouble alpha1 = alpha;
      FttComponent c;
      FttVector p;
      
      if (!child.c[i])
	g_assert_not_implemented ();
      ftt_cell_relative_pos (child.c[i], &p);
      for (c = 0; c < FTT_DIMENSION; c++) {
	alpha1 -= (&m.x)[c]*(0.25 + (&p.x)[c]);
	GFS_VARIABLE (child.c[i], t->m[c]->i) = (&m.x)[c];
      }
      GFS_VARIABLE (child.c[i], v->i) = gfs_plane_volume (&m, 2.*alpha1);
      GFS_VARIABLE (child.c[i], t->alpha->i) = 2.*alpha1;
    }
  }
}

static void vof_fine_coarse (FttCell * parent, GfsVariable * v)
{
  GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
  gfs_get_from_below_intensive (parent, v);
  gdouble f = GFS_VARIABLE (parent, v->i);
  FttComponent c;

  if (GFS_IS_FULL (f)) {
    for (c = 1; c < FTT_DIMENSION; c++)
      GFS_VARIABLE (parent, t->m[c]->i) = 0.;
    GFS_VARIABLE (parent, t->m[0]->i) = 1.;
    GFS_VARIABLE (parent, t->alpha->i) = f;
  }
  else {
    FttCellChildren child;
    FttVector m = {0., 0., 0.};
    guint i;
    
    ftt_cell_children (parent, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i]) {
	gdouble f = GFS_VARIABLE (child.c[i], v->i);
	gdouble a = f*(1. - f);

	for (c = 0; c < FTT_DIMENSION; c++)
	  (&m.x)[c] += a*GFS_VARIABLE (child.c[i], t->m[c]->i);
      }
    
    gdouble n = 0.;
    for (c = 0; c < FTT_DIMENSION; c++)
      n += fabs ((&m.x)[c]);
    if (n > 0.)
      for (c = 0; c < FTT_DIMENSION; c++)
	(&m.x)[c] /= n;
    else /* fixme: this is a small fragment */
      m.x = 1.;
    for (c = 0; c < FTT_DIMENSION; c++)
      GFS_VARIABLE (parent, t->m[c]->i) = (&m.x)[c];
    GFS_VARIABLE (parent, t->alpha->i) = gfs_plane_alpha (&m, f);
  }
}

static void variable_tracer_vof_init (GfsVariable * v)
{
  GFS_EVENT (v)->start = -1;
  GFS_EVENT (v)->istep = G_MAXINT/2;
  v->coarse_fine = vof_coarse_fine;
  v->fine_coarse = vof_fine_coarse;
  //  v->face_value = gfs_vof_face_value;
  GFS_VARIABLE_TRACER (v)->advection.cfl = 0.45;
}

GfsVariableClass * gfs_variable_tracer_vof_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_tracer_vof_info = {
      "GfsVariableTracerVOF",
      sizeof (GfsVariableTracerVOF),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_tracer_vof_class_init,
      (GtsObjectInitFunc) variable_tracer_vof_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_tracer_class ()), 
				  &gfs_variable_tracer_vof_info);
  }

  return klass;
}

typedef struct {
  GfsAdvectionParams * par;
  FttComponent c;
  GfsDomain * domain;
  guint depth, too_coarse;
} VofParms;

static gdouble plane_volume_shifted (FttVector m, gdouble alpha, FttVector p[2])
{
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++) {
    alpha -= (&m.x)[c]*(&p[0].x)[c];
    (&m.x)[c] *= (&p[1].x)[c] - (&p[0].x)[c];
  }
  return gfs_plane_volume (&m, alpha);
}

static gdouble fine_fraction (FttCellFace * face, VofParms * p, gdouble un)
{
  gdouble f = GFS_VARIABLE (face->cell, p->par->v->i);
  if (f == 0. || f == 1.)
    return f;
  else if (GFS_CELL_IS_BOUNDARY (face->cell))
    return GFS_STATE (face->cell)->f[face->d].v;
  else {
    FttComponent c;
    FttVector m;
    gdouble alpha = GFS_VARIABLE (face->cell, GFS_VARIABLE_TRACER_VOF (p->par->v)->alpha->i);

    for (c = 0; c < FTT_DIMENSION; c++)
      (&m.x)[c] = GFS_VARIABLE (face->cell, GFS_VARIABLE_TRACER_VOF (p->par->v)->m[c]->i);
    if (face->d % 2 != 0) {
      (&m.x)[face->d/2] = - (&m.x)[face->d/2];
      alpha += (&m.x)[face->d/2];
    }

    FttVector q[2] = {{0., 0., 0.},{1., 1., 1.}};

    (&q[0].x)[face->d/2] = 1. - un; (&q[1].x)[face->d/2] = 1.;
    return plane_volume_shifted (m, alpha, q);
  }
}

static gdouble coarse_fraction (FttCellFace * face, VofParms * p, gdouble un)
{
  gdouble f = GFS_VARIABLE (face->neighbor, p->par->v->i);
  if (f == 0. || f == 1.)
    return f;
  else {
    FttVector q[2] = {{0., 0., 0.},{1., 1., 1.}};
    FttComponent c;
    FttVector m, o;
    gdouble alpha = GFS_VARIABLE (face->neighbor, GFS_VARIABLE_TRACER_VOF (p->par->v)->alpha->i);
    
    for (c = 0; c < FTT_DIMENSION; c++)
      (&m.x)[c] = GFS_VARIABLE (face->neighbor, GFS_VARIABLE_TRACER_VOF (p->par->v)->m[c]->i);
    if (!FTT_FACE_DIRECT (face)) {
      (&m.x)[face->d/2] = - (&m.x)[face->d/2];
      alpha += (&m.x)[face->d/2];
    }
    
    /* shift interface perpendicularly */
    ftt_cell_relative_pos (face->cell, &o);
    for (c = 0; c < FTT_DIMENSION; c++)
      if (c != face->d/2) {
	(&q[0].x)[c] = (&o.x)[c] + 0.25;
	(&q[1].x)[c] = (&o.x)[c] + 0.75;
      }
    (&q[1].x)[face->d/2] = un;
    return plane_volume_shifted (m, alpha, q);
  }
}

#define TOO_COARSE(cell) (GFS_VARIABLE (cell, p->par->fv->i))

/* Marks coarse cells which should be refined because an interface in
   a neighboring finer cell will be advected into them */
static void face_too_coarse (FttCellFace * face, VofParms * p)
{
  if (ftt_face_type (face) == FTT_FINE_COARSE) {
    gdouble un = GFS_FACE_NORMAL_VELOCITY (face);
    if (!FTT_FACE_DIRECT (face))
      un = - un;
    if (un > 0.) {
      gdouble f = GFS_VARIABLE (face->neighbor, p->par->v->i);
      if (GFS_IS_FULL (f) &&
	  fine_fraction (face, p, un*p->par->dt/ftt_cell_size (face->cell)) != f) {
	p->too_coarse++;
	TOO_COARSE (face->neighbor) = TRUE;
      }
    }
  }
}

static void refine_too_coarse (FttCell * cell, VofParms * p)
{
  if (TOO_COARSE (cell)) {
    guint level = ftt_cell_level (cell);

    TOO_COARSE (cell) = FALSE;
    ftt_cell_refine_corners (cell, (FttCellInitFunc) gfs_cell_fine_init, p->domain);
    ftt_cell_refine_single (cell, (FttCellInitFunc) gfs_cell_fine_init, p->domain);
    if (level + 1 > p->depth)
      p->depth = level + 1;
  }
}

static void vof_face_value (FttCellFace * face, VofParms * p)
{
  gdouble un = GFS_FACE_NORMAL_VELOCITY (face)*p->par->dt/ftt_cell_size (face->cell);
  if (!FTT_FACE_DIRECT (face))
    un = - un;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE: {
    if (un < 0.) {
      FttCell * tmp = face->cell;
      face->cell = face->neighbor;
      face->neighbor = tmp;
      face->d = FTT_OPPOSITE_DIRECTION (face->d);
      un = - un;
    }
    GFS_STATE (face->cell)->f[face->d].v = fine_fraction (face, p, un);
    break;
  }
  case FTT_FINE_COARSE: {
    GFS_STATE (face->cell)->f[face->d].v = 
      un > 0. ? fine_fraction (face, p, un) : coarse_fraction (face, p, -un/2.);
    break;
  }
  default:
    g_assert_not_reached ();
  }
}
  
static void vof_flux (FttCellFace * face, VofParms * p)
{
  gdouble un = GFS_FACE_NORMAL_VELOCITY (face)*p->par->dt/ftt_cell_size (face->cell);
  if (!FTT_FACE_DIRECT (face))
    un = - un;
  if (fabs (un) > 0.5) {
    FttVector p;
    ftt_face_pos (face, &p);
    g_warning ("CFL (%g) at (%g,%g,%g) is larger than 0.5!", un, p.x, p.y, p.z);
  }
  un *= GFS_FACE_FRACTION (face);

  GFS_VARIABLE (face->cell, p->par->fv->i) += GFS_VARIABLE (face->cell, p->par->v->i)*un;
  gdouble flux = GFS_STATE (face->cell)->f[face->d].v*un;
  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE: {
    if (un < 0.)
      flux = GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v*un;
    GFS_VARIABLE (face->neighbor, p->par->fv->i) +=
      flux - GFS_VARIABLE (face->neighbor, p->par->v->i)*un;
    break;
  }
  case FTT_FINE_COARSE: {
    GFS_VARIABLE (face->neighbor, p->par->fv->i) += 
      (flux - coarse_fraction (face, p, 1.)*un)/FTT_CELLS;
    break;
  }
  default:
    g_assert_not_reached ();
  }
  GFS_VARIABLE (face->cell, p->par->fv->i) -= flux;
}

static void clamp (FttCell * cell, GfsVariable * v)
{
  gdouble f = GFS_VARIABLE (cell, v->i);
  GFS_VARIABLE (cell, v->i) = f < 1e-10 ? 0. : f > 1. - 1e-10 ? 1. : f;
}

/**
 * gfs_tracer_vof_advection:
 * @domain: a #GfsDomain.
 * @par: the advection parameters.
 *
 * Advects the @v field of @par using the current face-centered (MAC)
 * velocity field.
 */
void gfs_tracer_vof_advection (GfsDomain * domain,
			       GfsAdvectionParams * par)
{
  VofParms p;
  static FttComponent cstart = 0;
  FttComponent c;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (GFS_IS_VARIABLE_TRACER_VOF (par->v));
  g_return_if_fail (par->cfl <= 0.5);

  gfs_domain_timer_start (domain, "tracer_vof_advection");

  p.par = par;
  par->fv = gfs_temporary_variable (domain);
  for (c = 0; c < FTT_DIMENSION; c++) {
    p.c = (cstart + c) % FTT_DIMENSION;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, par->fv);
    p.too_coarse = 0;
    gfs_domain_face_traverse (domain, p.c,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) face_too_coarse, &p);
    if (p.too_coarse > 0) {
      p.depth = 0;
      p.domain = domain;
      gfs_domain_cell_traverse (domain,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) refine_too_coarse, &p);
      gfs_domain_reshape (domain, p.depth);
    }
    gfs_domain_face_traverse (domain, p.c,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) vof_face_value, &p);
    gfs_domain_face_bc (domain, p.c, par->v);
    gfs_domain_face_traverse (domain, p.c,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) vof_flux, &p);
    gfs_domain_traverse_merged (domain, (GfsMergedTraverseFunc) gfs_advection_update, par);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
    			      (FttCellTraverseFunc) clamp, par->v);
    gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) par->v->fine_coarse, par->v);
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, par->v);

    variable_tracer_vof_update (p.par->v, domain);
  }
  cstart = (cstart + 1) % FTT_DIMENSION;
  gts_object_destroy (GTS_OBJECT (par->fv));
  par->fv = NULL;

  gfs_domain_timer_stop (domain, "tracer_vof_advection");
}

static gdouble face_value (FttCell * cell, FttDirection d, GfsVariable * v)
{
  gdouble f = GFS_VARIABLE (cell, v->i);

  if (GFS_IS_FULL (f))
    return f;
  else {
    GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
    gdouble alpha = GFS_VARIABLE (cell, t->alpha->i);
    FttComponent c;
    FttVector m;
    
    for (c = 0; c < FTT_DIMENSION; c++)
      (&m.x)[c] = GFS_VARIABLE (cell, t->m[c]->i);
    (&m.x)[d/2] /= 2.;
    if (d % 2)
      alpha -= (&m.x)[d/2];
    return gfs_plane_volume (&m, alpha);
  }
}

/**
 * gfs_vof_face_value:
 * @face: a #FttCellFace.
 * @t: a #GfsVariableTracerVOF.
 *
 * Returns: the value of the VOF fraction defined by @t, interpolated
 * on @face.
 */
gdouble gfs_vof_face_value (const FttCellFace * face, GfsVariableTracerVOF * t)
{
  g_return_val_if_fail (face != NULL, 0.);
  g_return_val_if_fail (t != NULL, 0.);

  GfsVariable * v = GFS_VARIABLE1 (t);
  gdouble vright, vleft = GFS_VARIABLE (face->cell, v->i); //face_value (face->cell, face->d, v);
  if (ftt_face_type (face) == FTT_FINE_COARSE) {
    gdouble f = GFS_VARIABLE (face->neighbor, v->i);

    if (GFS_IS_FULL (f))
      vright = f;
    else {
      gdouble alpha = GFS_VARIABLE (face->neighbor, t->alpha->i);
      FttComponent c;
      FttVector m;

      for (c = 0; c < FTT_DIMENSION; c++)
	(&m.x)[c] = GFS_VARIABLE (face->neighbor, t->m[c]->i);

      FttVector p, o;
      ftt_face_pos (face, &p);
      ftt_cell_pos (face->neighbor, &o);
      gdouble h = ftt_cell_size (face->neighbor);

      (&p.x)[face->d/2] += face->d % 2 ? -h/4. : h/4.;
      for (c = 0; c < FTT_DIMENSION; c++)
	alpha -= (&m.x)[c]*(0.25 - ((&p.x)[c] - (&o.x)[c])/h);
      //      for (c = 0; c < FTT_DIMENSION; c++)
      //	(&m.x)[c] /= 2.;
      //      (&m.x)[face->d/2] /= 2.;
      //      if (!(face->d % 2))
      //	alpha -= (&m.x)[face->d/2];
      vright = gfs_plane_volume (&m, 2.*alpha);
#if 0
      if (vright > 0.2 && vright < 0.8) {
	fprintf (stderr, "%d (%g,%g) (%g,%g) %g\n", face->d, p.x, p.y, o.x, o.y, vright);
	g_assert_not_reached ();
      }
#endif
    }
  }
  else
    vright = GFS_VARIABLE (face->neighbor, v->i); //face_value (face->neighbor, FTT_OPPOSITE_DIRECTION (face->d), v);
  return (vright + vleft)/2.;
}

/**
 * gfs_vof_facet:
 * @cell: a #FttCell.
 * @t: a #GfsVariableTracerVOF.
 * @p: a #FttVector array (of size 2 in 2D and 6 in 3D)
 * @m: a #FttVector.
 *
 * Fills @p with the coordinates of points defining the
 * VOF-reconstructed interface facet defined by @t.
 *
 * Fills @m with the normal to the interface.
 *
 * Returns: the number of points defining the facet.
 */
guint gfs_vof_facet (FttCell * cell,
		     GfsVariableTracerVOF * t,
		     FttVector * p,
		     FttVector * m)
{
  g_return_val_if_fail (cell != NULL, 0);
  g_return_val_if_fail (t != NULL, 0);
  g_return_val_if_fail (p != NULL, 0);
  g_return_val_if_fail (m != NULL, 0);

  if (GFS_IS_FULL (GFS_VARIABLE (cell, GFS_VARIABLE1 (t)->i)))
    return 0;

  guint n = 0;
  FttVector q;
  ftt_cell_pos (cell, &q);
  gdouble h = ftt_cell_size (cell);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&m->x)[c] = GFS_VARIABLE (cell, t->m[c]->i);
  gdouble alpha = GFS_VARIABLE (cell, t->alpha->i);

#if FTT_2D
  gdouble x, y;

  if (fabs (m->y) > 1e-4) {
    y = (alpha - m->x)/m->y;
    if (y >= 0. && y <= 1.) {
      p[n].x = q.x + h/2.; p[n].y = q.y + h*(y - 0.5); p[n++].z = 0.;
    }
  }
  if (fabs (m->x) > 1e-4) {
    x = (alpha - m->y)/m->x;
    if (x >= 0. && x <= 1.) {
      p[n].x = q.x + h*(x - 0.5); p[n].y = q.y + h/2.; p[n++].z = 0.;
    }
  }
  if (fabs (m->y) > 1e-4) {
    y = alpha/m->y;
    if (y >= 0. && y <= 1.) {
      p[n].x = q.x - h/2.; p[n].y = q.y + h*(y - 0.5); p[n++].z = 0.;
    }
  }
  if (fabs (m->x) > 1e-4) {
    x = alpha/m->x;
    if (x >= 0. && x <= 1.) {
      p[n].x = q.x + h*(x - 0.5); p[n].y = q.y - h/2.; p[n++].z = 0.;
    }
  }
  g_assert (n <= 2);
#else /* 3D */
  gdouble max = fabs (m->x);
  c = FTT_X;
  if (fabs (m->y) > max) {
    max = fabs (m->y);
    c = FTT_Y;
  }
  if (fabs (m->z) > max)
    c = FTT_Z;
  q.x -= h/2.; q.y -= h/2.; q.z -= h/2.;
  (&q.x)[c] += h*alpha/(&m->x)[c];
  FttVector m1 = *m;
  gts_vector_normalize (&m1.x);

  FttDirection d[12];
  n = gfs_cut_cube_vertices (cell, -1, &q, &m1, p, d, NULL, NULL);
  g_assert (n <= 6);
#endif /* 3D */
  return n;
}

/**
 * gfs_vof_center:
 * @cell: a #FttCell.
 * @t: a #GfsVariableTracerVOF.
 * @p: a #FttVector.
 *
 * Fills @p with the coordinates of the center of mass of the
 * VOF-reconstructed interface facet defined by @t.
 *
 * Returns: the area (length in 2D) of the VOF-reconstructed facet or 0. if the
 * cell is not cut by the interface.
 */
gdouble gfs_vof_center (FttCell * cell, GfsVariableTracerVOF * t, FttVector * p)
{
  g_return_val_if_fail (cell != NULL, FALSE);
  g_return_val_if_fail (t != NULL, FALSE);
  g_return_val_if_fail (p != NULL, 0);

  if (GFS_IS_FULL (GFS_VARIABLE (cell, GFS_VARIABLE1 (t)->i)))
    return 0.;

  FttVector m, o;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&m.x)[c] = GFS_VARIABLE (cell, t->m[c]->i);
  gdouble area = gfs_plane_area_center (&m, GFS_VARIABLE (cell, t->alpha->i), p);
  ftt_cell_pos (cell, &o);
  gdouble h = ftt_cell_size (cell);
  for (c = 0; c < FTT_DIMENSION; c++)
    (&p->x)[c] = (&o.x)[c] + h*((&p->x)[c] - 0.5);
  return area;
}

static gdouble fraction (FttVector * p,
			 guint level,
			 GfsVariable * v)
{
  FttCell * cell = domain_and_boundary_locate (v->domain, *p, level);
  if (cell)
    return gfs_vof_interpolate (cell, p, level, GFS_VARIABLE_TRACER_VOF (v));
  else /* fixme: boundary conditions? */
    return 2.;
}

#define NMAX 10

#define ADD_H(f) { H += f; n++; }
#define SIGN(x) ((x) > 0. ? 1. : -1.)

static gdouble local_height (FttVector * p,
			     FttVector * origin,
			     guint level,
			     GfsVariable * v,
			     FttDirection d,
			     GtsVector interface)
{
  gdouble h = ftt_level_size (level), h1 = d % 2 ? - h : h, H = 0.;
  gdouble right = fraction (p, level, v), left = right;
  FttVector pright = *p, pleft = pright;
  FttComponent c = d/2;
  guint n = 0;

  ADD_H (right);
  gboolean found_interface = (right > 0.);
  while (n < NMAX && (!found_interface || !GFS_IS_FULL (right))) {
    (&pright.x)[c] += h1;
    right = fraction (&pright, level, v);
    if (right > 1.)
      return G_MAXDOUBLE;
    ADD_H (right);
    if (!GFS_IS_FULL (right))
      found_interface = TRUE;
  }
  if (right != 1.)
    return G_MAXDOUBLE;

  found_interface = (left < 1.);
  while (n < NMAX && (!found_interface || !GFS_IS_FULL (left))) {
    (&pleft.x)[c] -= h1;
    left = fraction (&pleft, level, v);
    if (left > 1.)
      return G_MAXDOUBLE;
    ADD_H (left);
    if (!GFS_IS_FULL (left))
      found_interface = TRUE;
  }
  if (left != 0.)
    return G_MAXDOUBLE;

  H -= ((&pright.x)[c] - (&origin->x)[c])/h1 + 0.5;
  interface[0] = (p->x - origin->x)/h;
  interface[1] = (p->y - origin->y)/h;
  interface[2] = (p->z - origin->z)/h;
  interface[c] = - SIGN (h1)*H;
  return H;
}

/* fixme: does not work for periodic boundary conditions along direction c
 * for cells close to the boundary */
static gboolean curvature_along_direction (FttCell * cell, 
					   GfsVariableTracerVOF * t,
					   FttComponent c,
					   gdouble * kappa,
					   gdouble * kmax,
					   GtsVector * interface,
					   guint * n)
{
  GfsVariable * v = GFS_VARIABLE1 (t);

  FttVector m;
  FttComponent i;
  for (i = 0; i < FTT_DIMENSION; i++)
    (&m.x)[i] = GFS_VARIABLE (cell, t->m[i]->i);
  FttDirection d = 2*c + ((&m.x)[c] > 0.);

  FttVector p;
  ftt_cell_pos (cell, &p);
  guint level = ftt_cell_level (cell);
  gdouble size = ftt_level_size (level), H;

  gboolean found_all_heights = TRUE;
  H = local_height (&p, &p, level, v, d, interface[*n]);
  if (H == G_MAXDOUBLE)
    found_all_heights = FALSE;
  else
    (*n)++;
#if 0
  if (H < -0.5 || H > 0.5)
    found_all_heights = FALSE;
#endif

#ifdef FTT_2D
  FttComponent cp = FTT_ORTHOGONAL_COMPONENT (c);
  FttVector q = p;
  gdouble h[2];
  (&q.x)[cp] += size;
  h[0] = local_height (&q, &p, level, v, d, interface[*n]);
  if (h[0] == G_MAXDOUBLE)
    found_all_heights = FALSE;
  else
    (*n)++;

  q = p;
  (&q.x)[cp] -= size;
  h[1] = local_height (&q, &p, level, v, d, interface[*n]);
  if (h[1] == G_MAXDOUBLE)
    found_all_heights = FALSE;
  else
    (*n)++;

  if (found_all_heights) {
    gdouble hxx = h[0] - 2.*H + h[1];
    gdouble hx = (h[0] - h[1])/2.;
    gdouble dnm = 1. + hx*hx;
    *kappa = hxx/(size*sqrt (dnm*dnm*dnm));
    if (kmax)
      *kmax = fabs (*kappa);
  }
#else  /* 3D */  
  static FttComponent or[3][2] = { { FTT_Y, FTT_Z }, { FTT_X, FTT_Z }, { FTT_X, FTT_Y } };
  gdouble h[3][3];
  gint x, y;

  for (x = -1; x <= 1; x++)
    for (y = -1; y <= 1; y++)
      if (x != 0 || y != 0) {
	FttVector q = p;
	(&q.x)[or[c][0]] += size*x;
	(&q.x)[or[c][1]] += size*y;
	h[x + 1][y + 1] = local_height (&q, &p, level, v, d, interface[*n]);
	if (h[x + 1][y + 1] == G_MAXDOUBLE)
	  found_all_heights = FALSE;
	else
	  (*n)++;
      }

  if (found_all_heights) {
    gdouble hxx = h[2][1] - 2.*H + h[0][1];
    gdouble hyy = h[1][2] - 2.*H + h[1][0];
    gdouble hx = (h[2][1] - h[0][1])/2.;
    gdouble hy = (h[1][2] - h[1][0])/2.;
    gdouble hxy = (h[2][2] + h[0][0] - h[2][0] - h[0][2])/4.;
    gdouble dnm = 1. + hx*hx + hy*hy; 
    *kappa = (hxx + hyy + hxx*hy*hy + hyy*hx*hx - 2.*hxy*hx*hy)/(size*sqrt (dnm*dnm*dnm));  
    if (kmax) {
      gdouble km = *kappa/2.;
      /* Gaussian curvature */
      gdouble kg = (hxx*hyy - hxy*hxy)/(size*size*dnm*dnm);
      gdouble a = km*km - kg;
      g_assert (a >= 0.);
      *kmax = fabs (km) + sqrt (a);
    }
  }
#endif /* 3D */
  return found_all_heights;
}

#define PARABOLA_FIT_CENTER_WEIGHT .1
#define PARABOLA_SIMPLER 0

typedef struct {
  GtsVector o;
#if FTT_2D /* y = a[0]*x^2 + a[0]*x + a[1] */
  GtsVector m;
  GtsMatrix * M;
  GtsVector rhs, a;
#else /* 3D */
# if PARABOLA_SIMPLER /* z = a[0]*x^2 + a[1]*y^2 + a[2]*x*y */
  GtsMatrix * M;
  GtsVector rhs, a;
# else /* z = a[0]*x^2 + a[1]*y^2 + a[2]*x*y + a[3]*x + a[4]*y + a[5] */
  gdouble ** M, rhs[6], a[6];
# endif
  GtsVector t[3];
#endif /* 3D */
} ParabolaFit;

static void parabola_fit_init (ParabolaFit * p, FttVector * o, FttVector * m)
{
  p->o[0] = o->x; p->o[1] = o->y; p->o[2] = o->z;
#if FTT_2D
  p->m[0] = m->x; p->m[1] = m->y; p->m[2] = m->z;
  gts_vector_normalize (p->m);
  p->M = gts_matrix_zero (NULL);
  p->rhs[0] = p->rhs[1] = p->rhs[2] = 0.;
#else /* 3D */
  gdouble max;
  GtsVector nx = {0., 0., 0.}, ny, nz;
  guint d = 0;

  nz[0] = m->x; nz[1] = m->y; nz[2] = m->z;
  gts_vector_normalize (nz);
  max = nz[0]*nz[0];
  /* build a vector orthogonal to nz */
  if (nz[1]*nz[1] > max) { max = nz[1]*nz[1]; d = 1; }
  if (nz[2]*nz[2] > max) d = 2;
  switch (d) {
  case 0: nx[0] = - nz[2]/nz[0]; nx[2] = 1.0; break;
  case 1: nx[1] = - nz[2]/nz[1]; nx[2] = 1.0; break;
  case 2: nx[2] = - nz[0]/nz[2]; nx[0] = 1.0; break;
  }
  gts_vector_normalize (nx);

  /* build a second vector orthogonal to nx and nz */
  gts_vector_cross (ny, nz, nx);

  /* transformation matrix from (i,j,k) to (nx, ny, nz) */
  p->t[0][0] = nx[0]; p->t[0][1] = nx[1]; p->t[0][2] = nx[2];
  p->t[1][0] = ny[0]; p->t[1][1] = ny[1]; p->t[1][2] = ny[2];
  p->t[2][0] = nz[0]; p->t[2][1] = nz[1]; p->t[2][2] = nz[2];

# if PARABOLA_SIMPLER
  p->M = gts_matrix_zero (NULL);
  p->rhs[0] = p->rhs[1] = p->rhs[2] = 0.;
# else
  p->M = gfs_matrix_new (6, sizeof (gdouble));
  p->rhs[0] = p->rhs[1] = p->rhs[2] = p->rhs[3] = p->rhs[4] = p->rhs[5] = 0.;
# endif
#endif /* 3D */
}

static void parabola_fit_add (ParabolaFit * p, GtsVector m, gdouble w)
{
#if FTT_2D
  gdouble x1 = m[0] - p->o[0];
  gdouble y1 = m[1] - p->o[1];
  gdouble x = p->m[1]*x1 - p->m[0]*y1;
  gdouble y = p->m[0]*x1 + p->m[1]*y1;
  gdouble x2 = w*x*x, x3 = x2*x, x4 = x3*x;
  p->M[0][0] += x4;
  p->M[1][0] += x3; p->M[1][1] += x2;
  p->M[2][1] += w*x; p->M[2][2] += w;
  p->rhs[0] += x2*y;
  p->rhs[1] += w*x*y;
  p->rhs[2] += w*y;
#else /* 3D */
  gdouble x1 = m[0] - p->o[0];
  gdouble y1 = m[1] - p->o[1];
  gdouble z1 = m[2] - p->o[2];
  gdouble x = p->t[0][0]*x1 + p->t[0][1]*y1 + p->t[0][2]*z1;
  gdouble y = p->t[1][0]*x1 + p->t[1][1]*y1 + p->t[1][2]*z1;
  gdouble z = p->t[2][0]*x1 + p->t[2][1]*y1 + p->t[2][2]*z1;
  gdouble x2 = x*x, x3 = x2*x, x4 = x3*x;
  gdouble y2 = y*y, y3 = y2*y, y4 = y3*y;
# if PARABOLA_SIMPLER
  p->M[0][0] += w*x4;
  p->M[1][0] += w*x2*y2; p->M[1][1] += w*y4;
  p->M[2][0] += w*x3*y;  p->M[2][1] += w*x*y3;
  p->rhs[0] += w*z*x2;   p->rhs[1] += w*z*y2;  p->rhs[2] += w*z*x*y;
# else
  p->M[0][0] += w*x4; p->M[1][1] += w*y4; p->M[2][2] += w*x2*y2; 
  p->M[3][3] += w*x2; p->M[4][4] += w*y2; p->M[5][5] += w;
  p->M[0][2] += w*x3*y; p->M[0][3] += w*x3; p->M[0][4] += w*x2*y;
  p->M[1][2] += w*x*y3; p->M[1][3] += w*x*y2; p->M[1][4] += w*y3;
  p->M[2][5] += w*x*y;
  p->M[3][5] += w*x;
  p->M[4][5] += w*y;
  p->rhs[0] += w*x2*z; p->rhs[1] += w*y2*z; p->rhs[2] += w*x*y*z;
  p->rhs[3] += w*x*z; p->rhs[4] += w*y*z; p->rhs[5] += w*z;
# endif
#endif /* 3D */
}

static void parabola_fit_solve (ParabolaFit * p)
{
#if FTT_2D
  p->M[0][1] = p->M[1][0];
  p->M[0][2] = p->M[2][0] = p->M[1][1];
  p->M[1][2] = p->M[2][1];
  GtsMatrix * M = gts_matrix3_inverse ((GtsMatrix *) p->M);
  if (M) {
    p->a[0] = M[0][0]*p->rhs[0] + M[0][1]*p->rhs[1] + M[0][2]*p->rhs[2];
    p->a[1] = M[1][0]*p->rhs[0] + M[1][1]*p->rhs[1] + M[1][2]*p->rhs[2];
    gts_matrix_destroy (M);
  }
  else /* this may be a degenerate/isolated interface fragment */
    p->a[0] = p->a[1] = 0.;
#else /* 3D */
# if PARABOLA_SIMPLER
  p->M[0][1] = p->M[1][0]; p->M[0][2] = p->M[2][0];
  p->M[1][2] = p->M[2][1]; p->M[2][2] = p->M[1][0];
  GtsMatrix * M = gts_matrix3_inverse ((GtsMatrix *) p->M);
  if (M) {
    p->a[0] = M[0][0]*p->rhs[0] + M[0][1]*p->rhs[1] + M[0][2]*p->rhs[2];
    p->a[1] = M[1][0]*p->rhs[0] + M[1][1]*p->rhs[1] + M[1][2]*p->rhs[2];
    p->a[2] = M[2][0]*p->rhs[0] + M[2][1]*p->rhs[1] + M[2][2]*p->rhs[2];
    gts_matrix_destroy (M);
  }
  else /* this may be a degenerate/isolated interface fragment */
    p->a[0] = p->a[1] = p->a[2] = 0.;
# else
  p->M[0][1] = p->M[2][2]; p->M[0][5] = p->M[3][3];
  p->M[1][5] = p->M[4][4];
  p->M[2][3] = p->M[0][4]; p->M[2][4] = p->M[1][3];
  p->M[3][4] = p->M[2][5];
  guint i, j;
  for (i = 1; i < 6; i++)
    for (j = 0; j < i; j++)
      p->M[i][j] = p->M[j][i];
  if (gfs_matrix_inverse (p->M, 6, 1e-10)) {
    for (i = 0; i < 6; i++) {
      p->a[i] = 0.;
      for (j = 0; j < 6; j++)
	p->a[i] += p->M[i][j]*p->rhs[j];
    }
  }
  else { /* this may be a degenerate/isolated interface fragment */
    g_warning ("singular matrix");
    p->a[0] = p->a[1] = p->a[2] = 0.;
  }
# endif
#endif /* 3D */
}

static gdouble parabola_fit_curvature (ParabolaFit * p, gdouble kappamax,
				       gdouble * kmax)
{
  gdouble kappa;
#if FTT_2D
  gdouble dnm = 1. + p->a[1]*p->a[1];
  kappa = 2.*p->a[0]/sqrt (dnm*dnm*dnm);
  if (kmax)
    *kmax = fabs (kappa);
#else /* 3D */
  gdouble hxx = 2.*p->a[0];
  gdouble hyy = 2.*p->a[1];
  gdouble hxy = p->a[2];
  gdouble hx, hy;
# if PARABOLA_SIMPLER
  hx = hy = 0.;
# else
  hx = p->a[3];
  hy = p->a[4];
# endif
  gdouble dnm = 1. + hx*hx + hy*hy;
  kappa = (hxx + hyy + hxx*hy*hy + hyy*hx*hx - 2.*hxy*hx*hy)/sqrt (dnm*dnm*dnm);
  if (kmax) {
    gdouble kg = (hxx*hyy - hxy*hxy)/(dnm*dnm);
    gdouble a = kappa*kappa/4. - kg;
    g_assert (a >= 0.);
    *kmax = fabs (kappa/2.) + sqrt (a);
  }
#endif /* 3D */
  if (fabs (kappa) > kappamax) {
    if (kmax)
      *kmax = kappamax;
    return kappa > 0. ? kappamax : - kappamax;
  }
  return kappa;
}

static void parabola_fit_destroy (ParabolaFit * p)
{
#if (FTT_2D || PARABOLA_SIMPLER)
  gts_matrix_destroy (p->M);
#else
  gfs_matrix_free (p->M);
#endif
}

static void add_vof_center (FttCell * cell, FttVector * p, guint level,
			    FttVector * origin,
			    GfsVariableTracerVOF * t,
			    ParabolaFit * fit, gdouble w)
{
  gdouble f = GFS_VARIABLE (cell, GFS_VARIABLE1 (t)->i);
  if (!GFS_IS_FULL (f)) {
    FttVector m, c;
    gdouble alpha = gfs_vof_plane_interpolate (cell, p, level, t, &m);
    gdouble area = gfs_plane_area_center (&m, alpha, &c);
    gdouble h = ftt_level_size (level);
    FttComponent i;
    for (i = 0; i < FTT_DIMENSION; i++)
      (&c.x)[i] = ((&p->x)[i] - (&origin->x)[i])/h + (&c.x)[i] - 0.5;
    parabola_fit_add (fit, &c.x, w*area);
  }
}

static void fit_from_fractions (FttCell * cell, GfsVariable * v, ParabolaFit * fit)
{
  gdouble h = ftt_cell_size (cell);
  guint level = ftt_cell_level (cell);
  gint x, y, z = 0;
  FttVector p;
  
  ftt_cell_pos (cell, &p);
#if !FTT_2D
  for (z = -1; z <= 1; z++)
#endif
    for (x = -1; x <= 1; x++)
      for (y = -1; y <= 1; y++)
	if (x != 0 || y != 0 || z != 0) {
	  FttVector o;
	  o.x = p.x + h*x; o.y = p.y + h*y; o.z = p.z + h*z;
	  FttCell * neighbor = domain_and_boundary_locate (v->domain, o, level);
	  if (neighbor)
	    add_vof_center (neighbor, &o, level, &p, GFS_VARIABLE_TRACER_VOF (v),
			    fit, 1.);
	}
}

/**
 * gfs_fit_curvature:
 * @cell: a #FttCell containing an interface.
 * @v: a #GfsVariableTracerVOF.
 * @kmax: a pointer or %NULL. 
 *
 * Computes an approximation of the curvature of the interface
 * contained in @cell using paraboloid fitting of the centroids of the
 * reconstructed interface segments.
 *
 * If @kmax is not %NULL, it is filled with the absolute value of the
 * maximum surface curvature (note that in 2D this is just the absolute value of
 * the mean curvature).
 *
 * Returns: (double in 3D) the mean curvature of the interface contained in @cell.
 */
gdouble gfs_fit_curvature (FttCell * cell, GfsVariableTracerVOF * t, gdouble * kmax)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (t != NULL, 0.);

  GfsVariable * v = GFS_VARIABLE1 (t);
  g_return_val_if_fail (!GFS_IS_FULL (GFS_VARIABLE (cell,  v->i)), 0.);

  FttVector m;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&m.x)[c] = GFS_VARIABLE (cell, t->m[c]->i);

  ParabolaFit fit;
  FttVector p, fc;
  ftt_cell_pos (cell, &p);
  gdouble area = gfs_vof_center (cell, t, &fc);
  gdouble h = ftt_cell_size (cell);
  fc.x = (fc.x - p.x)/h;
  fc.y = (fc.y - p.y)/h;
  fc.z = (fc.z - p.z)/h;
  parabola_fit_init (&fit, &fc, &m);
  parabola_fit_add (&fit, &fc.x, area);
  fit_from_fractions (cell, GFS_VARIABLE1 (t), &fit);
  parabola_fit_solve (&fit);
  gdouble kappa = parabola_fit_curvature (&fit, 2., kmax)/h;
  if (kmax)
    *kmax /= h;
  parabola_fit_destroy (&fit);
  return kappa;
}

#if FTT_2D
# define NI 3
#else
# define NI 9
#endif

static void orientation (FttVector * m, FttComponent * c)
{
  FttComponent i, j;
  for (i = 0; i < FTT_DIMENSION; i++)
    c[i] = i;
  for (i = 0; i < FTT_DIMENSION - 1; i++)
    for (j = 0; j < FTT_DIMENSION - 1 - i; j++)
      if (fabs ((&m->x)[c[j + 1]]) > fabs ((&m->x)[c[j]])) {
	FttComponent tmp = c[j];
	c[j] = c[j + 1];
	c[j + 1] = tmp;
      }
}

static guint independent_positions (GtsVector * interface, guint n)
{
  if (n < 2)
    return n;

  guint j, ni = 1;
  for (j = 1; j < n; j++) {
    guint i;
    gboolean depends = FALSE;
    for (i = 0; i < j && !depends; i++) {
      gdouble d2 = 0.;
      FttComponent c;
      for (c = 0; c < FTT_DIMENSION; c++)
	d2 += (interface[i][c] - interface[j][c])*(interface[i][c] - interface[j][c]);
      depends = d2 < 0.5*0.5;
    }
    ni += !depends;
  }
  return ni;
}

/**
 * gfs_height_curvature:
 * @cell: a #FttCell containing an interface.
 * @v: a #GfsVariableTracerVOF.
 * @kmax: a pointer or %NULL.
 *
 * An implementation of the Height-Function (HF) method generalised to
 * adaptive meshes.
 *
 * If @kmax is not %NULL, it is filled with the absolute value of the
 * maximum surface curvature (note that in 2D this is just the absolute value of
 * the mean curvature).
 *
 * Returns: (double in 3D) the mean curvature of the interface contained in @cell.
 */
gdouble gfs_height_curvature (FttCell * cell, GfsVariableTracerVOF * t, gdouble * kmax)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (t != NULL, 0.);

  GfsVariable * v = GFS_VARIABLE1 (t);
  gdouble f = GFS_VARIABLE (cell,  v->i);
  g_return_val_if_fail (!GFS_IS_FULL (f), 0.);

  FttVector m;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&m.x)[c] = GFS_VARIABLE (cell, t->m[c]->i);

  FttComponent try[FTT_DIMENSION];
  orientation (&m, try); /* sort directions according to normal */

  gdouble kappa = 0.;
  GtsVector interface[FTT_DIMENSION*NI];
  guint n = 0;
  for (c = 0; c < FTT_DIMENSION; c++) /* try each direction */
    if (curvature_along_direction (cell, t, try[c], &kappa, kmax, interface, &n))
      return kappa;

  /* Could not compute curvature from the simple algorithm along any direction:
   * Try circle fitting of the collected interface positions */

  if (independent_positions (interface, n) < 3*(FTT_DIMENSION - 1))
    return G_MAXDOUBLE;

  gdouble h = ftt_cell_size (cell);
  ParabolaFit fit;
  guint j;
  
  FttVector p, fc;
  ftt_cell_pos (cell, &p);
  gdouble area = gfs_vof_center (cell, t, &fc);
  fc.x = (fc.x - p.x)/h;
  fc.y = (fc.y - p.y)/h;
  fc.z = (fc.z - p.z)/h;
  parabola_fit_init (&fit, &fc, &m);
#if FTT_2D
  parabola_fit_add (&fit, &fc.x, PARABOLA_FIT_CENTER_WEIGHT);
#elif !PARABOLA_SIMPLER
  parabola_fit_add (&fit, &fc.x, area*100.);
#endif
  for (j = 0; j < n; j++)
    parabola_fit_add (&fit, interface[j], 1.);
  parabola_fit_solve (&fit);
  kappa = parabola_fit_curvature (&fit, 2., kmax)/h;
  if (kmax)
    *kmax /= h;
  parabola_fit_destroy (&fit);
  return kappa;
}
