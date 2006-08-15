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
  FttVector n;
  gdouble alpha1, a, amax, v;
  gdouble * md;
  guint j;

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
  if (n.z < 0.) {
    alpha1 -= n.z;
    n.z = - n.z;
  }

  if (alpha1 <= 0.)
    return 0.;

  if (alpha1 >= n.x + n.y + n.z)
    return 1.;

  n.x += 1e-4; n.y += 1e-4; n.z += 1e-4;
  amax = n.x + n.y + n.z;

  md = &n.x;
  v = alpha1*alpha1*alpha1;

  for (j = 0; j < 3; j++) {
    a = alpha1 - md[j];
    if (a > 0.)
      v -= a*a*a;
  }

  amax = alpha1 - amax;
  for (j = 0; j < 3; j++) {
    a = amax + md[j];
    if (a > 0.)
      v += a*a*a;
  }

  return v/(6.*n.x*n.y*n.z);
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

#if (!FTT_2D)
static gdouble plane_volume_derivative_ratio (FttVector * m, 
					      gdouble alpha, 
					      gdouble c)
{
  gdouble amax, v, vp;
  gdouble * md;
  guint j;

  md = &m->x;
  vp = alpha*alpha;
  v = alpha*vp;

  for (j = 0; j < 3; j++) {
    gdouble a = alpha - md[j];

    if (a > 0.) {
      vp -= a*a;
      v -= a*a*a;
    }
  }

  amax = alpha - m->x - m->y - m->z;
  for (j = 0; j < 3; j++) {
    gdouble a = amax + md[j];

    if (a > 0.) {
      vp += a*a;
      v += a*a*a;
    }
  }

  return (v - c)/(3.*vp);
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
  gdouble alpha, dalpha;
  FttVector n;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (c >= 0. && c <= 1., 0.);

  n.x = fabs (m->x); n.y = fabs (m->y); n.z = fabs (m->z);
  if (n.x*n.y*n.z < 1e-9)
    alpha = c;
  else {
    c *= 6.*n.x*n.y*n.z;
    alpha = (n.x + n.y + n.z)/2.;
    do {
      dalpha = plane_volume_derivative_ratio (&n, alpha, c);
      alpha -= dalpha;
    } while (fabs (dalpha) > 1e-6);
  }
  if (m->x < 0.)
    alpha += m->x;
  if (m->y < 0.)
    alpha += m->y;
  if (m->z < 0.)
    alpha += m->z;

  return alpha;
}
#endif /* 3D */

/**
 * gfs_youngs_normal:
 * @cell: a #FttCell.
 * @v: a #GfsVariable.
 * @n: a #FttVector.
 *
 * Fills @n with the Youngs-averaged gradients of @v 
 * normalised by the size of @cell.
 */
void gfs_youngs_normal (FttCell * cell, GfsVariable * v, FttVector * n)
{
  gdouble h = ftt_cell_size (cell), f[3][3];
  guint level = ftt_cell_level (cell);
  FttVector p;
  gint x, y;

  f[1][1] = GFS_VARIABLE (cell, v->i);
  ftt_cell_pos (cell, &p);
  for (x = -1; x <= 1; x++)
    for (y = -1; y <= 1; y++)
      if (x != 0 || y != 0) {
	FttVector o;
	o.x = p.x + h*x; o.y = p.y + h*y; o.z = 0.;
	FttCell * neighbor = gfs_domain_locate (v->domain, o, level);
	g_assert (neighbor);

	guint l = ftt_cell_level (neighbor);
	FttVector m;
	gdouble alpha;
	if (l == level || !gfs_vof_plane (neighbor, v, &m, &alpha))
	  f[x + 1][y + 1] = GFS_VARIABLE (neighbor, v->i);
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
	  f[x + 1][y + 1] = gfs_plane_volume (&m, 2.*alpha);
	}
      }
  n->x = (2.*f[2][1] + f[2][2] + f[2][0] - f[0][2] - 2.*f[0][1] - f[0][0])/8.;
  n->y = (2.*f[1][2] + f[2][2] + f[0][2] - f[2][0] - 2.*f[1][0] - f[0][0])/8.;
} 

void gfs_column_normal (FttCell * cell, GfsVariable * v, FttVector * n)
{
  gdouble h = ftt_cell_size (cell), f[3][3], s1, s2;
  guint level = ftt_cell_level (cell);
  FttVector p;
  gint x, y;

  f[1][1] = GFS_VARIABLE (cell, v->i);
  ftt_cell_pos (cell, &p);
  for (x = -1; x <= 1; x++)
    for (y = -1; y <= 1; y++)
      if (x != 0 || y != 0) {
	FttVector o;
	o.x = p.x + h*x; o.y = p.y + h*y; o.z = 0.;
	FttCell * neighbor = gfs_domain_locate (v->domain, o, level);
	//g_assert (neighbor);

	if (neighbor) {
	  guint l = ftt_cell_level (neighbor);
	  FttVector m;
	  gdouble alpha;
	  if (l == level || !gfs_vof_plane (neighbor, v, &m, &alpha))
	    f[x + 1][y + 1] = GFS_VARIABLE (neighbor, v->i);
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
	    f[x + 1][y + 1] = gfs_plane_volume (&m, 2.*alpha);
	  }
	}
	else
	  f[x + 1][y + 1] =  f[1][1];
      }
  s1 = f[2][0] + f[2][1] + f[2][2] - f[0][0] - f[0][1] - f[0][2];
  s2 = f[0][2] + f[1][2] + f[2][2] - f[0][0] - f[1][0] - f[2][0];
  if (fabs (s2) > fabs (s1)) {
    n->x = s1; n->y = s2 < 0. ? -2. : 2.;
  }
  else {
    n->y = s2; n->x = s1 < 0. ? - 2. : 2.;
  }
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
    FttVector q[2] = {{0., 0., 0.},{1., 1., 1.}};
    FttComponent c;
    FttVector m;
    gdouble alpha = GFS_VARIABLE (face->cell, p->alpha->i);
    
    for (c = 0; c < FTT_DIMENSION; c++)
      (&m.x)[c] = GFS_VARIABLE (face->cell, p->m[c]->i);
    if (!FTT_FACE_DIRECT (face)) {
      (&m.x)[face->d/2] = - (&m.x)[face->d/2];
      alpha += (&m.x)[face->d/2];
    }
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
    /* fixme: boundary conditions for (m, alpha) */
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
    for (c = 0; c < FTT_DIMENSION; c++)
      (&m->x)[c] /= n;
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

#define WIDTH 3

gdouble gfs_height_curvature (FttCell * cell, GfsVariable * v)
{
  gdouble f[2*WIDTH + 1][2*WIDTH + 1], f1[2*WIDTH + 1][2*WIDTH + 1], h = ftt_cell_size (cell);
  guint level = ftt_cell_level (cell);
  FttVector pos;
  gint x, y;

  ftt_cell_pos (cell, &pos);
  for (x = -WIDTH; x <= WIDTH; x++)
    for (y = -WIDTH; y <= WIDTH; y++) {
      FttVector o;
      o.x = pos.x + h*x; o.y = pos.y + h*y; o.z = 0.;
      FttCell * neighbor = gfs_domain_locate (v->domain, o, level);

      g_assert (neighbor);
      g_assert (ftt_cell_level (neighbor) == level);

      f[x + WIDTH][y + WIDTH] = GFS_VARIABLE (neighbor, v->i);
    }
  g_assert (f[WIDTH][WIDTH] > 0. && f[WIDTH][WIDTH] < 1.);

  gdouble s1 = 0.;
  for (x = WIDTH - 1; x <= WIDTH + 1; x++)
    for (y = 0; y <= WIDTH - 1; y++) {
      s1 += f[x][y];
    }
  for (x = 0; x <= 2*WIDTH; x++)
    for (y = 0; y <= 2*WIDTH; y++)
      f1[x][y] = f[x][y];

  gdouble s2 = 0.;
  for (x = WIDTH - 1; x <= WIDTH + 1; x++)
    for (y = WIDTH + 1; y <= 2*WIDTH; y++) {
      s2 += f[x][y];
    }
  if (s2 > s1) {
    s1 = s2;
    for (x = 0; x <= 2*WIDTH; x++)
      for (y = 0; y <= 2*WIDTH; y++)
	f1[x][y] = f[x][y];
  }

  gdouble s3 = 0.;
  for (y = WIDTH - 1; y <= WIDTH + 1; y++)
    for (x = 0; x <= WIDTH - 1; x++)
      s3 += f[x][y];
  if (s3 > s1) {
    s1 = s3;
    for (x = 0; x <= 2*WIDTH; x++)
      for (y = 0; y <= 2*WIDTH; y++)
	f1[x][y] = f[y][x];
  }

  gdouble s4 = 0.;
  for (y = WIDTH - 1; y <= WIDTH + 1; y++)
    for (x = WIDTH + 1; x <= 2*WIDTH; x++)
      s4 += f[x][y];
  if (s4 > s1) {
    s1 = s4;
    for (x = 0; x <= 2*WIDTH; x++)
      for (y = 0; y <= 2*WIDTH; y++)
	f1[x][y] = f[y][x];
  }

  gdouble FL = 0., FC = 0., FR = 0.;
  for (y = 0; y <= 2*WIDTH; y++) {
    FL += f1[WIDTH - 1][y];
    FC += f1[WIDTH][y];
    FR += f1[WIDTH + 1][y];
  }

  gdouble p = fabs ((FR - FL)/2.);
  p = 1 + p*p;
  return (FL + FR - 2.*FC)/(sqrt (p*p*p)*h);
}
