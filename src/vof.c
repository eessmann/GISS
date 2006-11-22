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
#include "vof.h"
#include "variable.h"

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
  g_return_if_fail (m->x >= 0. && m->y >= 0.);

  if (alpha <= 0.) {
    p->x = p->y = 0.;
    return;
  }

  if (alpha >= m->x + m->y) {
    p->x = p->y = 0.5;
    return;
  }

  g_return_if_fail (a > 0. && a < 1.);

  n = *m; n.x += 1e-4; n.y += 1e-4;

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
  g_return_if_fail (m->x >= 0. && m->y >= 0. && m->z >= 0.);

  if (alpha <= 0.) {
    p->x = p->y = p->z = 0.;
    return;
  }

  if (alpha >= m->x + m->y + m->z) {
    p->x = p->y = p->z = 0.5;
    return;
  }

  g_return_if_fail (a > 0. && a < 1.);

  n = *m; n.x += 1e-4; n.y += 1e-4; n.z += 1e-4;

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

static FttCell * domain_and_boundary_locate (GfsDomain * domain, FttVector p, guint level)
{
  FttCell * cell = gfs_domain_locate (domain, p, level);
  if (cell)
    return cell;
  return gfs_domain_boundary_locate (domain, p, level);
}

static void stencil (FttCell * cell, GfsVariable * v, gdouble f[3][3][3])
{
  gdouble h = ftt_cell_size (cell);
  guint level = ftt_cell_level (cell);
  FttVector p;
  gint x, y, z = 0;
  
  f[1][1][1] = GFS_VARIABLE (cell, v->i);
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
	  if (!neighbor) /* fixme: boundary conditions */
	    f[x + 1][y + 1][z + 1] = f[1][1][1];
	  else {
	    guint l = ftt_cell_level (neighbor);
	    FttVector m;
	    gdouble alpha;
	    if (l == level || !gfs_vof_plane (neighbor, v, &m, &alpha))
	      f[x + 1][y + 1][z + 1] = GFS_VARIABLE (neighbor, v->i);
	    else {
	      FttComponent c;
	      FttVector q;
	      
	      g_assert (l == level - 1);
	      ftt_cell_pos (neighbor, &q);
	      for (c = 0; c < FTT_DIMENSION; c++) {
		gdouble a = ((&o.x)[c] - (&q.x)[c])/h;
		g_assert (fabs (a) == 0.5);
		alpha -= (&m.x)[c]*(0.25 - a/2.);
	      }
	      f[x + 1][y + 1][z + 1] = gfs_plane_volume (&m, 2.*alpha);
	    }
	  }
	}
}

/**
 * gfs_youngs_normal:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 * @n: a #FttVector.
 *
 * Fills @n with the components of the normal to the interface defined
 * by volume fraction @v. Youngs' method is used.
 */
void gfs_youngs_normal (FttCell * cell, GfsVariable * v, FttVector * n)
{
  gdouble f[3][3][3];

  stencil (cell, v, f);
#if FTT_2D
  n->x = (2.*f[2][1][1] + f[2][2][1] + f[2][0][1] - f[0][2][1] - 2.*f[0][1][1] - f[0][0][1])/8.;
  n->y = (2.*f[1][2][1] + f[2][2][1] + f[0][2][1] - f[2][0][1] - 2.*f[1][0][1] - f[0][0][1])/8.;
  n->z = 0.;
#else  /* 3D */
  gdouble mm1 = f[0][0][0] + f[0][0][2] + f[0][2][0] + f[0][2][2] +
    2.*(f[0][0][1] + f[0][2][1] + f[0][1][0] + f[0][1][2]) + 
    4.*f[0][1][1];
  gdouble mm2 = f[2][0][0] + f[2][0][2] + f[2][2][0] + f[2][2][2] + 
    2.*(f[2][0][1] + f[2][2][1] + f[2][1][0] + f[2][1][2]) + 
    4.*f[2][1][1];
  n->x = (mm2 - mm1)/32.;
    
  mm1 = f[0][0][0] + f[0][0][2] + f[2][0][0] + f[2][0][2] + 
    2.*(f[0][0][1] + f[2][0][1] + f[1][0][0] + f[1][0][2]) + 
    4.*f[1][0][1];
  mm2 = f[0][2][0] + f[0][2][2] + f[2][2][0] + f[2][2][2] + 
    2.*(f[0][2][1] + f[2][2][1] + f[1][2][0] + f[1][2][2]) + 
    4.*f[1][2][1];
  n->y = (mm2 - mm1)/32.;
                  
  mm1 = f[0][0][0] + f[0][2][0] + f[2][0][0] + f[2][2][0] +
    2.*(f[0][1][0] + f[2][1][0] + f[1][0][0] + f[1][2][0]) + 
    4.*f[1][1][0];
  mm2 = f[0][0][2] + f[0][2][2] + f[2][0][2] + f[2][2][2] + 
    2.*(f[0][1][2] + f[2][1][2] + f[1][0][2] + f[1][2][2]) + 
    4.*f[1][1][2];
  n->z = (mm2 - mm1)/32.;
#endif /* 3D */
}

static gdouble plane_volume_shifted (FttVector m, gdouble alpha, FttVector p[2])
{
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++) {
    alpha -= (&m.x)[c]*(&p[0].x)[c];
    (&m.x)[c] *= (&p[1].x)[c] - (&p[0].x)[c];
  }
  return gfs_plane_volume (&m, alpha);
}

typedef struct {
  GfsVariable * m[FTT_DIMENSION], * alpha;
  GfsAdvectionParams * par;
  FttComponent c;
} VofParms;

static void vof_plane (FttCell * cell, VofParms * p)
{
  FttVector m;
  gdouble alpha;

  if (GFS_IS_MIXED (cell))
    g_assert_not_implemented ();

  if (gfs_vof_plane (cell, p->par->v, &m, &alpha)) {
    FttComponent c;

    for (c = 0; c < FTT_DIMENSION; c++) {
      GFS_VARIABLE (cell, p->m[c]->i) = - (&m.x)[c];
      alpha -= (&m.x)[c];
    }
    GFS_VARIABLE (cell, p->alpha->i) = alpha;
  }
}

static gdouble fine_fraction (FttCellFace * face, VofParms * p, gdouble un)
{
  gdouble f = GFS_VARIABLE (face->cell, p->par->v->i);
  if (f == 0. || f == 1.)
    return f;
  else {
    FttComponent c;
    FttVector m;
    gdouble alpha = GFS_VARIABLE (face->cell, p->alpha->i);

    for (c = 0; c < FTT_DIMENSION; c++)
      (&m.x)[c] = GFS_VARIABLE (face->cell, p->m[c]->i);
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
    gdouble alpha = GFS_VARIABLE (face->neighbor, p->alpha->i);
    
    for (c = 0; c < FTT_DIMENSION; c++)
      (&m.x)[c] = GFS_VARIABLE (face->neighbor, p->m[c]->i);
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
  
static void vof_flux (FttCellFace * face, VofParms * p)
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
    gdouble f = fine_fraction (face, p, un);
    GFS_VARIABLE (face->cell, p->par->fv->i) += 
      (GFS_VARIABLE (face->cell, p->par->v->i) - f)*un;
    GFS_VARIABLE (face->neighbor, p->par->fv->i) -= 
      (GFS_VARIABLE (face->neighbor, p->par->v->i) - f)*un;
    break;
  }
  case FTT_FINE_COARSE: {
    GFS_VARIABLE (face->cell, p->par->fv->i) += GFS_VARIABLE (face->cell, p->par->v->i)*un;
    GFS_VARIABLE (face->neighbor, p->par->fv->i) -= coarse_fraction (face, p, 1.)*un/FTT_CELLS;

    gdouble flux = un > 0. ? fine_fraction (face, p, un)*un : coarse_fraction (face, p, -un/2.)*un;
    GFS_VARIABLE (face->cell, p->par->fv->i) -= flux;
    GFS_VARIABLE (face->neighbor, p->par->fv->i) += flux/FTT_CELLS;
    break;
  }
  default:
    g_assert_not_reached ();
  }
}

static void vof_update (FttCell * cell, GfsAdvectionParams * par)
{
  gdouble f = GFS_VARIABLE (cell, par->v->i) + GFS_VARIABLE (cell, par->fv->i);
  GFS_VARIABLE (cell, par->v->i) = f < 0. ? 0. : f > 1. ? 1. : f;
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
  g_return_if_fail (par->cfl <= 0.5);

  gfs_domain_timer_start (domain, "tracer_vof_advection");

  p.par = par;
  for (c = 0; c < FTT_DIMENSION; c++)
    p.m[c] = gfs_temporary_variable (domain);
  p.alpha = gfs_temporary_variable (domain);
  par->fv = gfs_temporary_variable (domain);
  for (c = 0; c < FTT_DIMENSION; c++) {
    p.c = (cstart + c) % FTT_DIMENSION;

    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) vof_plane, &p);
    /* fixme: this does not work for all types of boundary conditions */
    FttComponent i;
    for (i = 0; i < FTT_DIMENSION; i++)
      gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, p.m[i]);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, p.alpha);

    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, par->fv);
    gfs_domain_face_traverse (domain, p.c,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) vof_flux, &p);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) vof_update, par);
    gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) par->v->fine_coarse, par->v);
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, par->v);
  }
  cstart = (cstart + 1) % FTT_DIMENSION;
  for (c = 0; c < FTT_DIMENSION; c++)
    gts_object_destroy (GTS_OBJECT (p.m[c]));
  gts_object_destroy (GTS_OBJECT (p.alpha));
  gts_object_destroy (GTS_OBJECT (par->fv));
  par->fv = NULL;

  gfs_domain_timer_stop (domain, "tracer_vof_advection");
}

/**
 * gfs_vof_coarse_fine:
 * @parent: a #FttCell.
 * @v: a #GfsVariable.
 *
 * Fills the @v variable of the children of @parent with
 * VOF-interpolated values.
 */
void gfs_vof_coarse_fine (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  FttVector m;
  gdouble alpha;
  guint i;

  g_return_if_fail (parent != NULL);
  g_return_if_fail (!FTT_CELL_IS_LEAF (parent));
  g_return_if_fail (v != NULL);

  ftt_cell_children (parent, &child);
  if (gfs_vof_plane (parent, v, &m, &alpha))
    for (i = 0; i < FTT_CELLS; i++) {
      gdouble alpha1 = alpha;
      FttComponent c;
      FttVector p;

      ftt_cell_relative_pos (child.c[i], &p);
      for (c = 0; c < FTT_DIMENSION; c++)
	alpha1 -= (&m.x)[c]*(0.25 - (&p.x)[c]);
      GFS_VARIABLE (child.c[i], v->i) = gfs_plane_volume (&m, 2.*alpha1);
    }
  else {
    gdouble f = GFS_VARIABLE (parent, v->i);
    for (i = 0; i < FTT_CELLS; i++)
      GFS_VARIABLE (child.c[i], v->i) = f;
  }
}

/**
 * gfs_vof_plane:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 * @m: a #FttVector.
 * @alpha: a double.
 *
 * If @cell is cut by the interface defined by @v, fills @m and @alpha
 * with the equation of the VOF-reconstructed interface.
 *
 * Returns: %TRUE if @cell is cut by the interface, %FALSE otherwise.
 */
gboolean gfs_vof_plane (FttCell * cell, GfsVariable * v,
			FttVector * m, gdouble * alpha)
{
  gdouble f;

  g_return_val_if_fail (cell != NULL, FALSE);
  g_return_val_if_fail (v != NULL, FALSE);
  g_return_val_if_fail (m != NULL, FALSE);
  g_return_val_if_fail (alpha != NULL, FALSE);

  f = GFS_VARIABLE (cell, v->i);
  THRESHOLD (f);

  if (GFS_IS_FULL (f))
    return FALSE;
  else {
    FttComponent c;
    gdouble n = 0.;

    gfs_youngs_normal (cell, v, m);
    for (c = 0; c < FTT_DIMENSION; c++)
      n += fabs ((&m->x)[c]);
    if (n > 0.)
      for (c = 0; c < FTT_DIMENSION; c++)
	(&m->x)[c] /= n;
    else /* fixme: this is a small fragment */
      m->x = 1.;
    *alpha = gfs_plane_alpha (m, f);
    return TRUE;
  }
}

/**
 * gfs_vof_facet:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 *
 * Returns: a list of newly allocated #GtsPoint defining the
 * VOF-reconstructed interface facet defined by @v, or %NULL if @cell
 * is not cut by the interface.
 */
GSList * gfs_vof_facet (FttCell * cell, GfsVariable * v)
{
  FttVector m;
  gdouble alpha;

  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (v != NULL, NULL);

#if FTT_3D
  g_assert_not_implemented ();
#endif  

  if (!gfs_vof_plane (cell, v, &m, &alpha))
    return NULL;
  else {
    GSList * l = NULL;
    gdouble x, y;

    if (m.y != 0.) {
      y = alpha/m.y;
      if (y >= 0. && y <= 1.)
	l = g_slist_prepend (l, gts_point_new (gts_point_class (), 0.5, 0.5 - y, 0.));
    }
    if (m.x != 0.) {
      x = alpha/m.x;
      if (x >= 0. && x <= 1.)
	l = g_slist_prepend (l, gts_point_new (gts_point_class (), 0.5 - x, 0.5, 0.));
    }
    if (m.y != 0.) {
      y = (alpha - m.x)/m.y;
      if (y >= 0. && y <= 1.)
	l = g_slist_prepend (l, gts_point_new (gts_point_class (), -0.5, 0.5 - y, 0.));
    }
    if (m.x != 0.) {
      x = (alpha - m.y)/m.x;
      if (x >= 0. && x <= 1.)
	l = g_slist_prepend (l, gts_point_new (gts_point_class (), 0.5 - x, -0.5, 0.));
    }
    g_assert (g_slist_length (l) == 2);
    return l;
  }
}

/**
 * gfs_vof_interpolate:
 * @cell: a #FttCell containing location @p.
 * @p: the center of the virtual cell.
 * @level: the level of the virtual cell.
 * @v: a #GfsVariable (volume fraction).
 *
 * Computes the volume fraction of a virtual cell at @level centered
 * on @p.
 *
 * Returns: the volume fraction of the virtual cell.
 */
gdouble gfs_vof_interpolate (FttCell * cell,
			     FttVector * p,
			     guint level,
			     GfsVariable * v)
{
  FttVector m;
  gdouble alpha;
  guint l = ftt_cell_level (cell);

  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (l <= level, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  if (l == level || !gfs_vof_plane (cell, v, &m, &alpha))
    return GFS_VARIABLE (cell, v->i);
  else {
    gdouble h = ftt_level_size (level);
    gdouble H = ftt_cell_size (cell);
    FttComponent c;
    FttVector q;

    ftt_cell_pos (cell, &q);
    alpha *= H;
    for (c = 0; c < FTT_DIMENSION; c++)
      alpha -= (&m.x)[c]*((&q.x)[c] + H/2 - (&p->x)[c] - h/2.);
    return gfs_plane_volume (&m, alpha/h);
  }
}

static gdouble fraction (FttVector * p,
			 guint level,
			 GfsVariable * v)
{
  FttCell * cell = domain_and_boundary_locate (v->domain, *p, level);
  g_assert (cell); /* fixme: boundary conditions? */
  return gfs_vof_interpolate (cell, p, level, v);
}

#define NMAX 10

static guint local_height (FttVector * p,
			   FttVector * origin,
			   guint level,
			   GfsVariable * v,
			   FttComponent c,
			   gdouble * H)
{
  gdouble h = ftt_level_size (level);
  gdouble right = fraction (p, level, v), left = right;
  FttVector pright = *p, pleft = pright;
  guint n = 1;

  *H = right;
  while (fabs (right - left) < 1. && n < NMAX)
    if (GFS_IS_FULL (left)) {
      (&pright.x)[c] += h;
      right = fraction (&pright, level, v); n++;
      *H += right;
      if (fabs (right - left) < 1.) {
	(&pleft.x)[c] -= h;
	left = fraction (&pleft, level, v); n++;
	*H += left;
      }
    }
    else {
      (&pleft.x)[c] -= h;
      left = fraction (&pleft, level, v); n++;
      *H += left;
      if (fabs (right - left) < 1.) {
	(&pright.x)[c] += h;
	right = fraction (&pright, level, v); n++;
	*H += right;
      }
    }
  if (fabs (right - left) < 1.)
    return 0;
  if (left > 0.5)
    *H += ((&pleft.x)[c] - (&origin->x)[c])/h - 0.5;
  else {
    g_assert (right > 0.5);
    *H -= ((&pright.x)[c] - (&origin->x)[c])/h + 0.5;
  }
  return n;
}

static FttComponent orientation (FttVector * m)
{
  gdouble max = fabs (m->x);
  FttComponent c = FTT_X, i;
  for (i = 1; i < FTT_DIMENSION; i++)
    if (fabs ((&m->x)[i]) > max) {
      max = fabs ((&m->x)[i]);
      c = i;
    }
  return c;
}

/**
 * gfs_height_curvature:
 * @cell: a #FttCell containing an interface.
 * @v: a #GfsVariable.
 *
 * An implementation of the Height-Function (HF) method generalised to
 * adaptive meshes.
 *
 * Returns: the curvature of the interface contained in @cell or
 * G_MAXDOUBLE if the curvature cannot be computed using the HF
 * method.
 */
gdouble gfs_height_curvature (FttCell * cell, GfsVariable * v)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);
  g_return_val_if_fail (!GFS_IS_FULL (GFS_VARIABLE (cell,  v->i)), 0.);

  FttVector m;
  gfs_youngs_normal (cell, v, &m);
  FttComponent c = orientation (&m);

  FttVector p;
  ftt_cell_pos (cell, &p);
  guint level = ftt_cell_level (cell);
  gdouble H;
  if (!local_height (&p, &p, level, v, c, &H))
    return G_MAXDOUBLE;
  if (H < -0.5 || H > 0.5)
    return G_MAXDOUBLE;
      
  gdouble size = ftt_level_size (level);
#ifdef FTT_2D
  FttComponent cp = FTT_ORTHOGONAL_COMPONENT (c);
  gdouble h[2];
  FttVector q = p;
  gdouble slope = rint ((&m.x)[cp]/(&m.x)[c]);

  (&q.x)[cp] += size;
  (&q.x)[c] -= slope*size;
  //  fprintf (stderr, "\n  (%g %g) ", q.x, q.y);
  if (!local_height (&q, &p, level, v, c, &h[0])) {
    g_warning ("Failed to compute local height at (%g,%g)", q.x, q.y);
    return G_MAXDOUBLE;
  }

  q = p;
  (&q.x)[cp] -= size;
  (&q.x)[c] += slope*size;
  //  fprintf (stderr, "\n  (%g %g) ", q.x, q.y);
  if (!local_height (&q, &p, level, v, c, &h[1])) {
    g_warning ("Failed to compute local height at (%g,%g)", q.x, q.y);
    return G_MAXDOUBLE;
  }

  gdouble hxx = h[0] - 2*H + h[1];
  gdouble hx = (h[0] - h[1])/2.;
  gdouble dnm = 1. + hx*hx;
  //  fprintf (stderr, " %g\n", hxx/(size*sqrt (dnm*dnm*dnm)));
  return hxx/(size*sqrt (dnm*dnm*dnm));
#else  /* 3D */  
  static FttComponent or[3][2] = { { FTT_Y, FTT_Z }, { FTT_X, FTT_Z }, { FTT_X, FTT_Y } };
  gdouble h[3][3];
  gint x, y;

  for (x = -1; x <= 1; x++)
    for (y = -1; y <= 1; y++)
      if (x != 0 || y != 0) {
	FttVector q = p;
	gdouble slope = rint ((&m.x)[or[c][0]]/(&m.x)[c]*x + (&m.x)[or[c][1]]/(&m.x)[c]*y);

	(&q.x)[or[c][0]] += size*x;
	(&q.x)[or[c][1]] += size*y;
	(&q.x)[c] -= slope*size;
	guint n = local_height (&q, &p, level, v, c, &h[x + 1][y + 1]);
	if (!n) {
	  g_warning ("Failed to compute local height at (%g,%g,%g)", q.x, q.y, q.z);
	  return G_MAXDOUBLE;
	}
      }
  
  gdouble hxx = h[2][1] - 2.*H + h[0][1];
  gdouble hyy = h[1][2] - 2.*H + h[1][0];
  gdouble hx = (h[2][1] - h[0][1])/2.;
  gdouble hy = (h[1][2] - h[1][0])/2.;
  gdouble hxy = (h[2][2] + h[0][0] - h[2][0] - h[0][2])/2.;
  gdouble dnm = 1. + hx*hx + hy*hy; 
  return (hxx + hyy + hxx*hy*hy + hyy*hx*hx - hxy*hx*hy)/(size*sqrt (dnm*dnm*dnm));  
#endif /* 3D */
}
