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

static FttCell * domain_and_boundary_locate (GfsDomain * domain, FttVector p, guint level)
{
  FttCell * cell = gfs_domain_locate (domain, p, level);
  if (cell)
    return cell;
  return gfs_domain_boundary_locate (domain, p, level);
}

static void stencil (FttCell * cell, GfsVariable * v, gdouble f[3][3])
{
  gdouble h = ftt_cell_size (cell);
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
	FttCell * neighbor = domain_and_boundary_locate (v->domain, o, level);
	if (!neighbor) /* fixme: boundary conditions */
	  f[x + 1][y + 1] = f[1][1];
	else {
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
      }  
}

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
  gdouble f[3][3];

  stencil (cell, v, f);
  n->x = (2.*f[2][1] + f[2][2] + f[2][0] - f[0][2] - 2.*f[0][1] - f[0][0])/8.;
  n->y = (2.*f[1][2] + f[2][2] + f[0][2] - f[2][0] - 2.*f[1][0] - f[0][0])/8.;
  n->z = 0.;
} 

static void column_normal (gdouble f[3][3], FttVector * n)
{
  gdouble s1, s2;

  s1 = f[2][0] + f[2][1] + f[2][2] - f[0][0] - f[0][1] - f[0][2];
  s2 = f[0][2] + f[1][2] + f[2][2] - f[0][0] - f[1][0] - f[2][0];
  if (fabs (s2) > fabs (s1)) {
    n->x = s1; n->y = s2 < 0. ? -2. : 2.;
  }
  else {
    n->y = s2; n->x = s1 < 0. ? - 2. : 2.;
  }
  n->z = 0.;
}

void gfs_column_normal (FttCell * cell, GfsVariable * v, FttVector * n)
{
  gdouble f[3][3];

  stencil (cell, v, f);
  column_normal (f, n);
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

#if FTT_3D
  g_assert_not_implemented ();
#endif

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

gdouble gfs_fit_curvature (FttCell * cell, GfsVariable * v)
{
#if !FTT_2D
  g_assert_not_implemented ();
#endif

  gdouble f[3][3];
  static gdouble c[][3] = {
    { +8.809173e-3, +8.744777e-3, +6.411956e-3 },
    { -7.302646e-2, -4.773088e-2, -1.217681e-1 },
    { +1.476918e-1, +1.911265e-1, +1.213353e-1 },
    { +8.127531e-1, +1.051815e0,  +1.091137e0  },
    { +1.662098e-1, +9.073978e-2, +3.268512e-1 },
    { -3.158622e-1, -3.656622e-1, -1.658780e-1 },
    { -7.899257e-2, -1.242422e-1, +9.856956e-2 },
    { -2.954707e-1, -3.824833e-1, -2.427756e-1 },
    { -8.268886e-1, -2.037791e-1, -1.522460e-1 },
    { +1.961225e-2, -9.380094e-1, -8.981913e-1 },
    { -1.107777e-1, -6.043060e-2, -2.178954e-1 },
    { -4.160499e-4, -1.001143e-3, -2.737967e-4 },
    { -1.507973e-2, -5.859334e-2, -2.349028e0  },
    { -1.792648e-4, -3.742121e-4, -2.355644e-5 },
    { +8.268877e-1, +2.037529e-1, +1.522527e-1 },
    { +6.325443e-1, +7.331659e-1, +3.322431e-1 },
    { -6.501391e-2, +4.143641e-1, +2.703116e-1 },
    { +1.579836e-1, +2.482401e-1, -1.982763e-1 },
    { +2.944518e-6, +3.065368e-4, +1.363379e-3 },
    { +4.704595e-6, +1.256257e-4, -3.048371e-7 }
  };

  stencil (cell, v, f);
  gdouble F = f[1][1];
  
  gdouble S = 0.;
  guint i,j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      S += f[i][j];
  
  FttVector m;
  gdouble alpha;
  g_return_val_if_fail (gfs_vof_plane (cell, v, &m, &alpha), 0.);
  alpha = (alpha + m.x + m.y)/3.;
  S -= 9.*gfs_plane_volume (&m, alpha);

  gdouble M = atan (MIN (fabs (m.x), fabs (m.y))/MAX (fabs (m.x), fabs (m.y)));

  gdouble kappa = c[0][0] + (c[1][0]*F   + c[2][0]*M   + c[3][0]*S +
			     c[4][0]*F*F + c[5][0]*M*M + c[6][0]*S*S +
			     c[7][0]*F*M + c[8][0]*F*S + c[9][0]*M*S + 
			     c[10][0]*F*F*F + c[11][0]*M*M*M + c[12][0]*S*S*S +
			     c[13][0]*F*F*M + c[14][0]*F*F*S + c[15][0]*M*M*F +
			     c[16][0]*M*M*S + c[17][0]*S*S*F + c[18][0]*S*S*M +
			     c[19][0]*F*M*S);
  if (fabs (kappa) < 0.4)
    kappa = c[0][1] + (c[1][1]*F   + c[2][1]*M   + c[3][1]*S +
		       c[4][1]*F*F + c[5][1]*M*M + c[6][1]*S*S +
		       c[7][1]*F*M + c[8][1]*F*S + c[9][1]*M*S + 
		       c[10][1]*F*F*F + c[11][1]*M*M*M + c[12][1]*S*S*S +
		       c[13][1]*F*F*M + c[14][1]*F*F*S + c[15][1]*M*M*F +
		       c[16][1]*M*M*S + c[17][1]*S*S*F + c[18][1]*S*S*M +
		       c[19][1]*F*M*S);
  if (fabs (kappa) < 0.1)
    kappa = c[0][2] + (c[1][2]*F   + c[2][2]*M   + c[3][2]*S +
		       c[4][2]*F*F + c[5][2]*M*M + c[6][2]*S*S +
		       c[7][2]*F*M + c[8][2]*F*S + c[9][2]*M*S + 
		       c[10][2]*F*F*F + c[11][2]*M*M*M + c[12][2]*S*S*S +
		       c[13][2]*F*F*M + c[14][2]*F*F*S + c[15][2]*M*M*F +
		       c[16][2]*M*M*S + c[17][2]*S*S*F + c[18][2]*S*S*M +
		       c[19][2]*F*M*S);
  return kappa/ftt_cell_size (cell);
}

static gdouble fraction (FttVector p,
			 guint level,
			 GfsVariable * v)
{
  FttCell * cell = domain_and_boundary_locate (v->domain, p, level);
  g_assert (cell);
  return gfs_vof_interpolate (cell, p, level, v);
}

static guint local_height (FttVector p,
			   FttVector o,
			   guint level,
			   GfsVariable * v,
			   FttComponent c,
			   gdouble * H)
{
  gdouble h = ftt_level_size (level);
  gdouble f1 = fraction (p, level, v);
  guint n = !GFS_IS_FULL (f1);

  *H = f1;
  FttVector i = p;
  (&i.x)[c] -= h;
  f1 = fraction (i, level, v);
  while (!GFS_IS_FULL (f1)) {
    //    fprintf (stderr, " f1: %g (%g,%g)\n", f1, i.x, i.y);
    *H += f1; n++;
    (&i.x)[c] -= h;
    f1 = fraction (i, level, v);
  }
  if (f1 > 0.5)
    *H += ((&i.x)[c] - (&o.x)[c])/h + 0.5;

  i = p;
  (&i.x)[c] += h;
  gdouble f2 = fraction (i, level, v);
  while (!GFS_IS_FULL (f2)) {
    //    fprintf (stderr, " f2: %g (%g,%g)\n", f2, i.x, i.y);
    *H += f2; n++;
    (&i.x)[c] += h;
    f2 = fraction (i, level, v);
  }
  if (f2 > 0.5) {
    if (f1 > 0.5)
      return 0;
    *H = 0.5 + *H - ((&i.x)[c] - (&o.x)[c])/h;
  }
  else if (f1 < 0.5)
    return 0;

  //  fprintf (stderr, " %d %g ", n, *H);

  return n;
}

static gboolean local_height1 (FttVector p,
			       guint level,
			       GfsVariable * v,
			       FttComponent c,
			       gdouble * H)
{
  gdouble h = ftt_level_size (level);
  FttVector i = p;
  gint x;

  for (x = -2; x <= 2; x++) {
    (&i.x)[c] = (&p.x)[c] + h*x;
    FttCell * cell = gfs_domain_locate (v->domain, i, -1);
    gdouble alpha;
    FttVector m;
    if (gfs_vof_plane (cell, v, &m, &alpha)) {
      FttComponent cp = FTT_ORTHOGONAL_COMPONENT (c);
      if ((&m.x)[c] != 0.) {
	gdouble y = (alpha - (&m.x)[cp]/2.)/(&m.x)[c];
	if (y >= 0. && y <= 1.) {
	  *H = (&i.x)[c]/h + 0.5 - y;
	  return TRUE;
	}
      }
    }
  }

  return FALSE;
}

static gboolean height (FttVector p,
			FttCell * cell,
			GfsVariable * v,
			FttDirection d,
			gdouble * H)
{
  FttVector m;
  gdouble alpha;

  if (gfs_vof_plane (cell, v, &m, &alpha)) {
    FttComponent c = d/2;
    if ((&m.x)[c] != 0.) {
      FttComponent cp = FTT_ORTHOGONAL_COMPONENT (c);
      FttVector i;
      ftt_cell_pos (cell, &i);
      gdouble h = ftt_cell_size (cell);
      gdouble x = ((&p.x)[cp] - (&i.x)[cp])/h + 0.5;
      gdouble y = (alpha - (&m.x)[cp]*x)/(&m.x)[c];
      fprintf (stderr, "y %g\n", y);
      if (y >= 0. && y <= 1.) {
	*H = (&i.x)[c] + (0.5 - y)*h;
	return TRUE;
      }
    }
  }
  return FALSE;
}

static gboolean is_interface (FttCell * cell, GfsVariable * v)
{
  gdouble f = GFS_VARIABLE (cell, v->i);
  return !GFS_IS_FULL (f);
}

static gboolean local_height2 (FttCell * cell,
			       GfsVariable * v,
			       FttDirection d,
			       gdouble * H)
{
  FttVector p;
  ftt_cell_pos (cell, &p);
  fprintf (stderr, "** %d ", is_interface (cell, v));
  if (height (p, cell, v, d, H)) {
    fprintf (stderr, "height %g\n", *H);
    return TRUE;
  }

  gboolean interface = is_interface (cell, v);
  FttCell * right = cell, * left = cell;
  guint n = 0;
  do {
    right = ftt_cell_neighbor (right, d);
    g_assert (right);
    if (height (p, right, v, d, H)) {
      fprintf (stderr, "height1 %g\n", *H);
      return TRUE;
    }
    else {
      interface = interface || is_interface (right, v);
      left = ftt_cell_neighbor (left, FTT_OPPOSITE_DIRECTION (d));
      g_assert (left);
      if (height (p, left, v, d, H)) {
	fprintf (stderr, "height2 %g\n", *H);
	return TRUE;
      }
      interface = interface || is_interface (left, v);
    }
    n++;
  } while (is_interface (right, v) || is_interface (left, v) || (!interface && n < 2));
  fprintf (stderr, "No height\n");
  return FALSE;
}

typedef struct {
  FttVector p;
  guint n;
} InterfacePoints;

guint interface_points1 (FttCell * cell, GfsVariable * v, FttDirection d, InterfacePoints ip[6])
{
  guint N = 0;

  g_return_val_if_fail (cell != NULL, 0);
  g_return_val_if_fail (v != NULL, 0);

  FttVector p;
  ftt_cell_pos (cell, &p);
  guint level = ftt_cell_level (cell);
  gdouble h = ftt_level_size (level);

  FttComponent c = d/2;
  //  for (c = 0; c < FTT_DIMENSION; c++)
  {
    FttComponent cp = FTT_ORTHOGONAL_COMPONENT (c);
    gint i;
    for (i = -1.; i <= 1; i++) {
      FttVector q = p;
      guint n;

      (&q.x)[cp] += h*i;
      if ((n = local_height1 (q, level, v, d, &((&q.x)[c])))) {
	(&q.x)[c] *= h;
	ip[N].p = q;
	ip[N++].n = n;
      }
    }
  }

  guint i;
  for (i = 0; i < N; i++) {
    ip[i].p.x = (ip[i].p.x - p.x)/h;
    ip[i].p.y = (ip[i].p.y - p.y)/h;
  }
    
  return N;
}

guint interface_points (FttCell * cell, GfsVariable * v, FttDirection d, InterfacePoints ip[6])
{
  guint N = 0;

  g_return_val_if_fail (cell != NULL, 0);
  g_return_val_if_fail (v != NULL, 0);

  FttComponent c = d/2;
  FttComponent cp = FTT_ORTHOGONAL_COMPONENT (c);

  FttVector p;
  ftt_cell_pos (cell, &p);
  gdouble h = ftt_cell_size (cell);

  if (p.x > -0.2178 && p.y > -0.1547 && p.x < -0.189 && p.y < -0.1266) {
    gdouble H;
    gboolean ret = local_height2 (cell, v, d, &H);
    fprintf (stderr, "tyty %d %d %g %d\n", d, ret, H,
	     H > (&p.x)[c] - h/2. && H < (&p.x)[c] + h/2.);
  }

  ftt_cell_pos (cell, &ip[N].p);
  if (local_height2 (cell, v, d, &((&ip[N].p.x)[c])) &&
      (&ip[N].p.x)[c] > (&p.x)[c] - h/2. && (&ip[N].p.x)[c] < (&p.x)[c] + h/2.)
    N++;
  
  FttCell * neighbor = ftt_cell_neighbor (cell, 2*cp);
  g_assert (neighbor);
  ftt_cell_pos (neighbor, &ip[N].p);
  if (local_height2 (neighbor, v, d, &((&ip[N].p.x)[c])))
    N++;

  neighbor = ftt_cell_neighbor (cell, 2*cp + 1);
  g_assert (neighbor);
  ftt_cell_pos (neighbor, &ip[N].p);
  if (local_height2 (neighbor, v, d, &((&ip[N].p.x)[c])))
    N++;
    
  guint i;
  for (i = 0; i < N; i++) {
    ip[i].p.x = (ip[i].p.x - p.x)/h;
    ip[i].p.y = (ip[i].p.y - p.y)/h;
  }

  fprintf (stderr, "N %d\n", N);

  return N;
}

static FttDirection orientation (FttCell * cell, GfsVariable * v)
{
  FttVector m;      
  gfs_youngs_normal (cell, v, &m);
    
  gdouble max = fabs (m.x);
  FttComponent i;
  FttDirection d = m.x > 0. ? 0 : 1;
  for (i = 1; i < FTT_DIMENSION; i++)
    if (fabs ((&m.x)[i]) > max) {
      max = fabs ((&m.x)[i]);
      d = (&m.x)[i] > 0. ? 2*i : 2*i + 1;
    }
  return d;
}

static gdouble curvature (FttVector * a, FttVector * b, FttVector * c)
{
  gdouble xd, yd, xe, ye;
  gdouble xad, yad, xae, yae;
  gdouble det;
  
  xd = (a->x + b->x)/2.; yd = (a->y + b->y)/2.;
  xe = (a->x + c->x)/2.; ye = (a->y + c->y)/2.;
  xad = xd - a->x; yad = yd - a->y;
  xae = xe - a->x; yae = ye - a->y;
  det = xad*yae - xae*yad;
  if (det == 0.)
    return 0.;
  gdouble x = (yae*yad*(yd - ye) + xad*yae*xd - xae*yad*xe)/det - a->x;
  gdouble y = -(xae*xad*(xd - xe) + yad*xae*yd - yae*xad*ye)/det - a->y;
  return 1./sqrt (x*x + y*y);
}

static gboolean height_curvature (FttCell * cell, GfsVariable * v, FttDirection d, gdouble * k)
{
  InterfacePoints ip[6];
  FttComponent or = d/2;
  guint n = interface_points (cell, v, d, ip);

  if (n < 3)
    return FALSE;
#if 0
  *k = curvature (&ip[0].p, &ip[1].p, &ip[2].p)/ftt_cell_size (cell);
  return TRUE;
#else
  FttVector m;
  gfs_youngs_normal (cell, v, &m);
  gts_vector_normalize (&m.x);

  guint i;
  GtsMatrix * o = gts_matrix_zero (NULL);
  GtsVector rhs = {0., 0., 0.};

  gdouble x0 = 0., y0 = 0.;
#if 1
  for (i = 0; i < n; i++) {
    x0 += ip[i].p.x;
    y0 += ip[i].p.y;
  }
  x0 /= n; y0 /= n;
#else
  (&m.x)[or] = (&m.x)[or] > 0. ? 1. : -1.;
  (&m.x)[FTT_ORTHOGONAL_COMPONENT (or)] = 0.;
#endif

  for (i = 0; i < n; i++) {
    gdouble x = - m.y*(ip[i].p.x - x0) + m.x*(ip[i].p.y - y0);
    gdouble y = - m.x*(ip[i].p.x - x0) - m.y*(ip[i].p.y - y0);
    o[0][0] += x*x*x*x; o[0][1] += x*x*x; o[0][2] += x*x;
    o[1][2] += x;
    rhs[0] += x*x*y;
    rhs[1] += x*y;
    rhs[2] += y;
  }
  o[1][0] = o[0][1];
  o[1][1] = o[2][0] = o[0][2];
  o[2][1] = o[1][2];
  o[2][2] = n;

  GtsMatrix * im = gts_matrix3_inverse (o);
  g_assert (im);
  gts_matrix_destroy (o);
  
  gdouble a = im[0][0]*rhs[0] + im[0][1]*rhs[1] + im[0][2]*rhs[2];
  gdouble b = im[1][0]*rhs[0] + im[1][1]*rhs[1] + im[1][2]*rhs[2];
  gdouble c = im[2][0]*rhs[0] + im[2][1]*rhs[1] + im[2][2]*rhs[2];
  gts_matrix_destroy (im);

  gdouble kappa = 1. + b*b;
  kappa = 2.*a/(sqrt (kappa*kappa*kappa)*ftt_cell_size (cell));

  fprintf (stderr, "unset label\n");
  fprintf (stderr, "set title 'kappa=%g'\n", kappa);
  FttVector p;
  ftt_cell_pos (cell, &p);
  gdouble h = ftt_cell_size (cell);
  //  for (i = 0; i < n; i++)
  //    fprintf (stderr, "set label \"%d\" at %g,%g\n", ip[i].n, p.x + h*ip[i].p.x, p.y + h*ip[i].p.y);
  fprintf (stderr, "plot [%g:%g][%g:%g]'/tmp/toto.gnu' w l, '-' w p\n",
	   p.x - 3.*h, p.x + 3.*h, p.y - 3.*h, p.y + 3.*h);
  for (i = 0; i < n; i++)
    fprintf (stderr, "%g %g\n", p.x + h*ip[i].p.x, p.y + h*ip[i].p.y);
  fprintf (stderr, "%g %g\n", p.x, p.y);
  fprintf (stderr, "e\npause -1\n");

  *k = kappa;
  return TRUE;
#endif
}

gdouble gfs_shahriar1_curvature (FttCell * cell, GfsVariable * v)
{
  FttDirection d = orientation (cell, v);

  FttVector p;
  ftt_cell_pos (cell, &p);
  fprintf (stderr, "%g %g %d\n", p.x, p.y, d);

  gdouble k;

  if (height_curvature (cell, v, d, &k))
    return k;
  else
    return G_MAXDOUBLE;
}

gdouble gfs_shahriar_curvature (FttCell * cell, GfsVariable * v)
{
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (v != NULL, 0.);

  FttComponent c = orientation (cell, v)/2;

  FttVector p;
  ftt_cell_pos (cell, &p);
  //  fprintf (stderr, "(%g %g) %d ", p.x, p.y, c);
  guint level = ftt_cell_level (cell);
  gdouble size = ftt_level_size (level);
  gdouble H;
  if (!local_height (p, p, level, v, c, &H))
    return G_MAXDOUBLE;
  if (H < -0.5 || H > 0.5)
    return G_MAXDOUBLE;
      
  FttComponent cp = FTT_ORTHOGONAL_COMPONENT (c);
#ifdef FTT_2D
  gdouble h[2], hxx, hx;
  FttVector q = p, m;
  gfs_youngs_normal (cell, v, &m);
  gdouble slope = rint ((&m.x)[cp]/(&m.x)[c]);

  (&q.x)[cp] += size;
  (&q.x)[c] -= slope*size;
  //  fprintf (stderr, "\n  (%g %g) ", q.x, q.y);
  if (!local_height (q, p, level, v, c, &h[0]))
    g_assert_not_reached ();

  q = p;
  (&q.x)[cp] -= size;
  (&q.x)[c] += slope*size;
  //  fprintf (stderr, "\n  (%g %g) ", q.x, q.y);
  if (!local_height (q, p, level, v, c, &h[1]))
    g_assert_not_reached ();

  hxx = h[0] - 2*H + h[1];
  hx = (h[0] - h[1])/2.;
  gdouble dnm = 1. + hx*hx;
  //  fprintf (stderr, " %g\n", hxx/(size*sqrt (dnm*dnm*dnm)));
  return hxx/(size*sqrt (dnm*dnm*dnm));
#else  /* FTT_3D */  
  g_assert_not_implemented ();
#if 0
  static FttDirection dp[3][4] =
  {{ FTT_FRONT, FTT_BACK, FTT_BOTTOM, FTT_TOP },
   { FTT_RIGHT, FTT_LEFT, FTT_BACK, FTT_FRONT },
   { FTT_TOP, FTT_BOTTOM, FTT_LEFT, FTT_RIGHT }};  
  gdouble h[8], hxx, hyy, hx, hy, hxy;
    
  for (i = 0; i < 4 ; i++) {
    n = ftt_cell_neighbor (cell, d);
    g_assert (n);
    h[i] = local_height (n, v, c1);     
    n = ftt_cell_neighbor (n, dp[c1][i]);
    g_assert (n);
    h[i + 4] = local_height (n, v, c1);
    d = (d + 1)%6;   
  }   
  
  hxx = h[0] - 2.*height + h[1];
  hyy = h[2] - 2.*height + h[3];
  hx = (h[0] - h[1])/2.;
  hy = (h[2] - h[3])/2.;
  hxy = (h[4] + h[5] - h[6] - h[7])/2.;
  dnm = 1. + hx*hx + hy*hy;
  
  return (hxx + hyy + hxx*hy*hy + hyy*hx*hx - hxy*hx*hy)
    /(ftt_cell_size (cell)*sqrt (dnm*dnm*dnm));  
#endif
#endif  
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
			     FttVector p,
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
      alpha -= (&m.x)[c]*((&q.x)[c] + H/2 - (&p.x)[c] - h/2.);
    return gfs_plane_volume (&m, alpha/h);
  }
}
