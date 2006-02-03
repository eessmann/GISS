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
 * @c1: width of the cell.
 *
 * Returns: the area of the fraction of a rectangular cell (c1,1)
 * lying under the line (@m,@alpha).
 */
gdouble gfs_line_area (FttVector * m, gdouble alpha, gdouble c1)
{
  FttVector n;
  gdouble a, v;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (m->x >= 0. && m->y >= 0., 0.);

  if (alpha <= 0.)
    return 0.;

  if (alpha >= m->x*c1 + m->y || c1 == 0.)
    return c1;

  n = *m; n.x += 1e-4; n.y += 1e-4;

  v = alpha*alpha;

  a = alpha - n.x*c1;
  if (a > 0.)
    v -= a*a;

  a = alpha - n.y;
  if (a > 0.)
    v -= a*a;

  return v/(2.*n.x*n.y);
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
 * @c1: width of the cell.
 *
 * Returns: the volume of a parallelepipedic cell (c1,1,1) lying under
 * the plane (@m,@alpha).
 */
gdouble gfs_plane_volume (FttVector * m, gdouble alpha, gdouble c1)
{
  FttVector n;
  gdouble a, amax, v;
  gdouble * md;
  guint j;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (m->x >= 0. && m->y >= 0. && m->z >= 0., 0.);
  
  if (alpha <= 0.)
    return 0.;

  if (alpha >= m->x*c1 + m->y + m->z || c1 == 0.)
    return c1;

  n = *m; n.x += 1e-4; n.y += 1e-4; n.z += 1e-4;
  amax = n.x*c1 + n.y + n.z;

  md = &n.x;
  v = alpha*alpha*alpha;

  a = alpha - n.x*c1;
  if (a > 0.)
    v -= a*a*a;
  for (j = 1; j < 3; j++) {
    a = alpha - md[j];
    if (a > 0.)
      v -= a*a*a;
  }

  amax = alpha - amax;
  a = amax + n.x*c1;
  if (a > 0.)
    v += a*a*a;
  for (j = 1; j < 3; j++) {
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

static gdouble line_area_derivative_ratio (FttVector * m, 
					   gdouble alpha, 
					   gdouble c)
{
  gdouble a, v, vp;

  vp = alpha;
  v = vp*vp;

  a = alpha - m->x;
  if (a > 0.) {
    vp -= a;
    v -= a*a;
  }

  a = alpha - m->y;
  if (a > 0.) {
    vp -= a;
    v -= a*a;
  }

  return (v - c)/(2.*vp);
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
  gdouble alpha, dalpha;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (c >= 0. && c <= 1., 0.);

  if (m->x*m->y < 1e-6)
    return c;
  c *= 2.*m->x*m->y;
  alpha = (m->x + m->y)/2.;
  
  do {
    dalpha = line_area_derivative_ratio (m, alpha, c);
    alpha -= dalpha;
  } while (fabs (dalpha) > 1e-6);

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

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (c >= 0. && c <= 1., 0.);

  if (m->x*m->y*m->z < 1e-9)
    return c;
  c *= 6.*m->x*m->y*m->z;
  alpha = (m->x + m->y + m->z)/2.;
  
  do {
    dalpha = plane_volume_derivative_ratio (m, alpha, c);
    alpha -= dalpha;
  } while (fabs (dalpha) > 1e-6);

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
  g_return_if_fail (n != NULL);

  for (i = 0; i < (FTT_DIMENSION - 1)*4; i++)
    u[i] = gfs_cell_corner_value (cell, d[i], v, -1);

#if FTT_2D
  n->x = (u[0] + u[3] - u[1] - u[2])/2.;
  n->y = (u[0] + u[1] - u[2] - u[3])/2.;
#else  /* 3D */
  n->x = (u[0] + u[3] + u[4] + u[7] - u[1] - u[2] - u[5] - u[6])/4.;
  n->y = (u[0] + u[1] + u[4] + u[5] - u[2] - u[3] - u[6] - u[7])/4.;
  n->z = (u[0] + u[1] + u[2] + u[3] - u[4] - u[5] - u[6] - u[7])/4.;
#endif /* 3D */
}

static void gfs_cell_vof_advected_face_values (FttCell * cell,
					       FttComponent c,
					       GfsAdvectionParams * par)
{
  gdouble f, dt, u_right, u_left;
  FttDirection right = 2*c, left = 2*c + 1;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (par != NULL);

  f = GFS_VARIABLE (cell, par->v->i);
  THRESHOLD (f);

  dt = par->dt/ftt_cell_size (cell);
  u_right = GFS_STATE (cell)->f[right].un*dt;
  u_left = GFS_STATE (cell)->f[left].un*dt;

  if (GFS_IS_MIXED (cell))
    GFS_VARIABLE (cell, par->fv->i) = 
      f*(GFS_STATE (cell)->solid->s[right]*u_right - 
	 GFS_STATE (cell)->solid->s[left]*u_left);
  else
    GFS_VARIABLE (cell, par->fv->i) = f*(u_right - u_left);

  if (GFS_IS_FULL (f)) {
    GFS_STATE (cell)->f[right].v = f;
    GFS_STATE (cell)->f[left].v = f;
  }
  else {
    FttVector m1, m;
    gdouble alpha, n;
    static FttComponent d[FTT_DIMENSION][FTT_DIMENSION - 1] = {
#if FTT_2D
      { FTT_Y }, { FTT_X }
#else  /* 3D */
      { FTT_Y, FTT_Z}, { FTT_Z, FTT_X }, { FTT_X, FTT_Y }
#endif /* 3D */
    };
    guint i;

    gfs_youngs_normal (cell, par->v, &m1);
    m.x = - (&m1.x)[c];
    for (i = 0; i < FTT_DIMENSION - 1; i++)
      (&m.x)[i + 1] = - (&m1.x)[d[c][i]];
    
    if (m.x < 0.) {
      n = - u_left; u_left = - u_right; u_right = n;
      m.x = - m.x;
      left = 2*c;
      right = 2*c + 1;
    }

    /* normalize */
    n = 0.;
    for (i = 0; i < FTT_DIMENSION; i++) {
      (&m.x)[i] = fabs ((&m.x)[i]);
      n += (&m.x)[i];
    }
    for (i = 0; i < FTT_DIMENSION; i++)
      (&m.x)[i] /= n;
    
    alpha = gfs_plane_alpha (&m, f);
    
    /* advect interface */
    m.x /= 1. + u_right - u_left;
    alpha += m.x*u_left;

    if (u_left < 0.)
      GFS_STATE (cell)->f[left].v =
	- gfs_plane_volume (&m, alpha - m.x*u_left, - u_left)/u_left;
    else
      GFS_STATE (cell)->f[left].v = f;

    if (u_right > 0.)
      GFS_STATE (cell)->f[right].v =
	gfs_plane_volume (&m, alpha - m.x, u_right)/u_right;
    else
      GFS_STATE (cell)->f[right].v = f;
  }
}

static void save_previous (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsVariable * vh = data[1];

  GFS_VARIABLE (cell, vh->i) = GFS_VARIABLE (cell, v->i);
}

static void average_previous (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsVariable * vh = data[1];

  GFS_VARIABLE (cell, vh->i) = (GFS_VARIABLE (cell, vh->i) +
				GFS_VARIABLE (cell, v->i))/2.;
}

static void vof_face_values (FttCell * cell, gpointer * data)
{
  GfsAdvectionParams * par = data[0];
  FttComponent * c = data[1];

  gfs_cell_vof_advected_face_values (cell, *c, par);
}

/**
 * gfs_face_vof_advection_flux:
 * @face: a #FttCellFace.
 * @par: the advection parameters.
 *
 * Adds to variable @par->fv, the value of the (conservative)
 * advection flux of the face variable through @face.
 *
 * This function assumes that the face variable has been previously
 * defined using gfs_cell_vof_advected_face_values().
 */
static void gfs_face_vof_advection_flux (const FttCellFace * face,
					 const GfsAdvectionParams * par)
{
  gdouble un, flux = 0.;

  g_return_if_fail (face != NULL);
  g_return_if_fail (par != NULL);

  un = GFS_FACE_NORMAL_VELOCITY (face);
  if (!FTT_FACE_DIRECT (face))
    un = - un;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE: case FTT_FINE_COARSE:
    flux = un >= 0. ? GFS_STATE (face->cell)->f[face->d].v :
      GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v;
    break;
  default:
    g_assert_not_reached ();
  }

  flux *= GFS_FACE_FRACTION (face)*un*par->dt/ftt_cell_size (face->cell);
  GFS_VARIABLE (face->cell, par->fv->i) -= flux;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_VARIABLE (face->neighbor, par->fv->i) += flux;
    break;
  case FTT_FINE_COARSE:
    GFS_VARIABLE (face->neighbor, par->fv->i) += flux/FTT_CELLS;
    break;
  default:
    g_assert_not_reached ();
  }
}

/**
 * gfs_tracer_vof_advection:
 * @domain: a #GfsDomain.
 * @par: the advection parameters.
 * @half: a #GfsVariable or %NULL.
 *
 * Advects the @v field of @par using the current face-centered (MAC)
 * velocity field.
 *
 * If @half is not %NULL, the half-timestep value of @par->v is
 * stored in the corresponding variable.  
 */
void gfs_tracer_vof_advection (GfsDomain * domain,
			       GfsAdvectionParams * par,
			       GfsVariable * half)
{
  gpointer data[2];
  static FttComponent cstart = 0;
  FttComponent c, c1;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);

  gfs_domain_timer_start (domain, "tracer_vof_advection");

  if (half) {
    data[0] = par->v;
    data[1] = half;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) save_previous, data);
  }
  par->fv = gfs_temporary_variable (domain);
  data[0] = par;
  data[1] = &c1;
  for (c = 0; c < FTT_DIMENSION; c++) {
    c1 = (cstart + c) % FTT_DIMENSION;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) vof_face_values, data);
    gfs_domain_face_bc (domain, c1, par->v);
    gfs_domain_face_traverse (domain, c1,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
	(FttFaceTraverseFunc) gfs_face_vof_advection_flux, par);
    gfs_domain_traverse_merged (domain,
				(GfsMergedTraverseFunc) gfs_advection_update,
				par);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, par->v);
  }
  cstart = (cstart + 1) % FTT_DIMENSION;
  gts_object_destroy (GTS_OBJECT (par->fv));
  par->fv = NULL;

  if (half) {
    data[0] = par->v;
    data[1] = half;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) average_previous, data);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, half);
  }

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
  gdouble f;
  FttCellChildren child;
  guint i;

  g_return_if_fail (parent != NULL);
  g_return_if_fail (!FTT_CELL_IS_LEAF (parent));
  g_return_if_fail (v != NULL);

  ftt_cell_children (parent, &child);

  f = GFS_VARIABLE (parent, v->i);
  THRESHOLD (f);

  if (GFS_IS_FULL (f))
    for (i = 0; i < FTT_CELLS; i++)
      GFS_VARIABLE (child.c[i], v->i) = f;
  else {
    FttVector m, m1;
    FttComponent c;
    gdouble n = 0., alpha;

    gfs_youngs_normal (parent, v, &m);
    for (c = 0; c < FTT_DIMENSION; c++) {
      (&m1.x)[c] = fabs ((&m.x)[c]);
      n += (&m1.x)[c];
    }
    for (c = 0; c < FTT_DIMENSION; c++)
      (&m1.x)[c] /= n;
    alpha = gfs_plane_alpha (&m1, f);

    for (i = 0; i < FTT_CELLS; i++) {
      gdouble alpha1 = alpha;
      FttVector p;

      ftt_cell_relative_pos (child.c[i], &p);
      for (c = 0; c < FTT_DIMENSION; c++) {
	(&p.x)[c] = (&m.x)[c] < 0. ? (&p.x)[c] + 0.25 : 0.25 - (&p.x)[c];
	alpha1 -= (&m1.x)[c]*(&p.x)[c];
      }
      GFS_VARIABLE (child.c[i], v->i) = gfs_plane_volume (&m1, 2.*alpha1, 1.);
    }
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
    FttVector m1;
    FttComponent c;
    gdouble n = 0.;

    gfs_youngs_normal (cell, v, m);
    for (c = 0; c < FTT_DIMENSION; c++) {
      (&m1.x)[c] = fabs ((&m->x)[c]);
      n += (&m1.x)[c];
    }
    for (c = 0; c < FTT_DIMENSION; c++) {
      (&m1.x)[c] /= n;
      (&m->x)[c] /= n;
    }
    *alpha = gfs_plane_alpha (&m1, f);

    for (c = 0; c < FTT_DIMENSION; c++)
      if ((&m->x)[c] < 0.)
	*alpha += (&m->x)[c];
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
