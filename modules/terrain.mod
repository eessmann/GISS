/* Gerris - The GNU Flow Solver                       (-*-C-*-)
 * Copyright (C) 2001-2008 National Institute of Water and Atmospheric Research
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
#include <glob.h>
#if GSL
# include <gsl/gsl_linalg.h>
#endif
#include "refine.h"
#include "solid.h"
#include "rsurface.h"

static gchar * default_path = ".";

/* GfsRefineTerrain: Header */

typedef struct _GfsRefineTerrain         GfsRefineTerrain;

#define NM 4

struct _GfsRefineTerrain {
  /*< private >*/
  GfsRefine parent;
  guint level;
  gboolean refined;
  GfsVariable * type;

#if !FTT_2D
  GfsVariable * min, * max;
  gdouble front;
#endif

  RSurface ** rs;
  guint nrs;

  /*< public >*/
  gchar * name, * basename;  
  GfsVariable * h[NM], * he, * hn;
  GfsFunction * criterion;
};

#define GFS_REFINE_TERRAIN(obj)            GTS_OBJECT_CAST (obj,\
					           GfsRefineTerrain,\
					           gfs_refine_terrain_class ())
#define GFS_IS_REFINE_TERRAIN(obj)         (gts_object_is_from_class (obj,\
						   gfs_refine_terrain_class ()))
     
GfsRefineClass * gfs_refine_terrain_class  (void);

/* GfsRefineTerrain: Object */

typedef struct {
  FttVector c;
  FttVector p[4];
  gdouble min[2], max[2], h;
  GfsRefineTerrain * t;
  FttCell * cell;
} Polygon;

static void map_inverse (GfsSimulation * sim, FttVector * p, gdouble min[2], gdouble max[2])
{
  gfs_simulation_map_inverse (sim, p);
  if (p->x < min[0]) min[0] = p->x;
  if (p->x > max[0]) max[0] = p->x;
  if (p->y < min[1]) min[1] = p->y;
  if (p->y > max[1]) max[1] = p->y;
}

static void polygon_init (Polygon * p, FttCell * cell, GfsRefineTerrain * t)
{
  GfsSimulation * sim = gfs_object_simulation (t);
  FttVector q;
  ftt_cell_pos (cell, &q);
  p->cell = cell;
  p->t = t;
  p->c = q;
  p->h = ftt_cell_size (cell)/2.;
  p->min[0] = p->min[1] = G_MAXDOUBLE;
  p->max[0] = p->max[1] = - G_MAXDOUBLE;
  p->p[0].x = q.x + p->h; p->p[0].y = q.y + p->h;
  map_inverse (sim, &p->p[0], p->min, p->max);
  p->p[1].x = q.x - p->h; p->p[1].y = q.y + p->h;
  map_inverse (sim, &p->p[1], p->min, p->max);
  p->p[2].x = q.x - p->h; p->p[2].y = q.y - p->h;
  map_inverse (sim, &p->p[2], p->min, p->max);
  p->p[3].x = q.x + p->h; p->p[3].y = q.y - p->h;
  map_inverse (sim, &p->p[3], p->min, p->max);
}

static gboolean right (const double a[2], const double b[2], const double c[2])
{
  return (b[0] - a[0])*(c[1] - a[1]) - (b[1] - a[1])*(c[0] - a[0]) < 0.;
}

static gboolean polygon_contains (Polygon * p, gdouble q[2])
{
  if (right (&p->p[0].x, &p->p[1].x, q))
    return FALSE;
  if (right (&p->p[1].x, &p->p[2].x, q))
    return FALSE;
  if (right (&p->p[2].x, &p->p[3].x, q))
    return FALSE;
  if (right (&p->p[3].x, &p->p[0].x, q))
    return FALSE;
  return TRUE;
}

static gboolean polygon_includes (RSurfaceRect rect, Polygon * p)
{
  gdouble q[2];
  q[0] = rect[0].l; q[1] = rect[1].l;
  if (!polygon_contains (p, q))
    return FALSE;
  q[0] = rect[0].l; q[1] = rect[1].h;
  if (!polygon_contains (p, q))
    return FALSE;
  q[0] = rect[0].h; q[1] = rect[1].l;
  if (!polygon_contains (p, q))
    return FALSE;
  q[0] = rect[0].h; q[1] = rect[1].h;
  if (!polygon_contains (p, q))
    return FALSE;
  return TRUE;
}

static gboolean polygon_intersects (RSurfaceRect rect, Polygon * p)
{
  /* fixme: this could be improved? */
  return (rect[0].l <= p->max[0] && rect[0].h >= p->min[0] &&
	  rect[1].l <= p->max[1] && rect[1].h >= p->min[1]);
}

typedef struct {
  gdouble H[NM+1], m[NM][NM];
  gdouble h[NM], he, cond, min, max;
  GfsRefineTerrain * t;
  FttCell * cell;
  Polygon * p;
  gboolean relative;
} RMS;

static void rms_init (RMS * rms, Polygon * p, gboolean relative)
{
  guint i, j;
  for (i = 0; i < NM + 1; i++)
    rms->H[i] = 0.;
  for (i = 0; i < NM; i++)
    for (j = 0; j < NM; j++)
      rms->m[i][j] = 0.;
  rms->p = p;
  rms->t = p->t;
  rms->cell = p->cell;
  rms->relative = relative;
  rms->min = G_MAXDOUBLE;
  rms->max = - G_MAXDOUBLE;
}

static gdouble rms_minimum (RMS * rms)
{
  if (rms->m[0][0] == 0.)
    return 0.;
  return sqrt (fabs (rms->h[0]*(rms->h[0]*rms->m[0][0] + 
				2.*(rms->h[1]*rms->m[0][1] + 
				    rms->h[2]*rms->m[0][2] +
				    rms->h[3]*rms->m[0][3] - rms->H[0])) +
		     rms->h[1]*(rms->h[1]*rms->m[1][1] + 
				2.*(rms->h[2]*rms->m[1][2] +
				    rms->h[3]*rms->m[1][3] - rms->H[1])) +
		     rms->h[2]*(rms->h[2]*rms->m[2][2] +
				2.*(rms->h[3]*rms->m[2][3] - rms->H[2])) +
		     rms->h[3]*(rms->h[3]*rms->m[3][3] +
				- 2.*rms->H[3]) +
		     rms->H[4])/rms->m[0][0]);
}

static void function_from_corners (gdouble h[4], gdouble H[4])
{
  h[0] = (H[0] + H[1] + H[2] + H[3])/4.;
  h[1] = (H[0] - H[1] - H[2] + H[3])/4.;
  h[2] = (H[0] + H[1] - H[2] - H[3])/4.;
  h[3] = (H[0] - H[1] + H[2] - H[3])/4.;  
}

static gdouble cell_value (FttCell * cell, GfsVariable * h[NM], FttVector p)
{
  gdouble size = ftt_cell_size (cell)/2.;
  FttVector q;
  ftt_cell_pos (cell, &q);
  p.x = (p.x - q.x)/size;
  p.y = (p.y - q.y)/size;
  return GFS_VALUE (cell, h[0]) + 
    GFS_VALUE (cell, h[1])*p.x + 
    GFS_VALUE (cell, h[2])*p.y + 
    GFS_VALUE (cell, h[3])*p.x*p.y;
}

static void cell_coefficients (FttCell * cell, GfsVariable * h[NM], gdouble hp[NM])
{
  gdouble size = ftt_cell_size (cell)/2.;
  FttVector q;
  ftt_cell_pos (cell, &q);
  hp[0] = GFS_VALUE (cell, h[0]) 
    - (GFS_VALUE (cell, h[1])*q.x + GFS_VALUE (cell, h[2])*q.y
       - GFS_VALUE (cell, h[3])*q.x*q.y/size)/size;
  hp[1] = (GFS_VALUE (cell, h[1]) - GFS_VALUE (cell, h[3])*q.y/size)/size;
  hp[2] = (GFS_VALUE (cell, h[2]) - GFS_VALUE (cell, h[3])*q.x/size)/size;
  hp[3] = GFS_VALUE (cell, h[3])/(size*size);
}

static void corners_from_parent (FttCell * cell, GfsRefineTerrain * t, gdouble H[4])
{
  gdouble size = ftt_cell_size (cell);
  FttCell * parent = ftt_cell_parent (cell);
  FttVector p;
  ftt_cell_pos (cell, &p);
  p.x += size/2.; p.y += size/2.;
  H[0] = cell_value (parent, t->h, p);
  p.x -= size;
  H[1] = cell_value (parent, t->h, p);
  p.y -= size;
  H[2] = cell_value (parent, t->h, p);
  p.x += size;
  H[3] = cell_value (parent, t->h, p);
}

static gdouble clamp (gdouble x, gdouble min, gdouble max)
{
  if (x > max) return max;
  if (x < min) return min;
  return x;
}

static void variance_check (RMS * rms)
{
  g_assert (rms->m[0][0] >= NM);
  gdouble H[4];
  if (rms->relative) {
    gdouble H0[4];
    corners_from_parent (rms->cell, rms->t, H0);
    H[0] = clamp (rms->h[0] + rms->h[1] + rms->h[2] + rms->h[3], 
		  rms->min - H0[0], rms->max - H0[0]);
    H[1] = clamp (rms->h[0] - rms->h[1] + rms->h[2] - rms->h[3], 
		  rms->min - H0[1], rms->max - H0[1]);
    H[2] = clamp (rms->h[0] - rms->h[1] - rms->h[2] + rms->h[3], 
		  rms->min - H0[2], rms->max - H0[2]);
    H[3] = clamp (rms->h[0] + rms->h[1] - rms->h[2] - rms->h[3], 
		  rms->min - H0[3], rms->max - H0[3]);
  }
  else {
    H[0] = clamp (rms->h[0] + rms->h[1] + rms->h[2] + rms->h[3], rms->min, rms->max);
    H[1] = clamp (rms->h[0] - rms->h[1] + rms->h[2] - rms->h[3], rms->min, rms->max);
    H[2] = clamp (rms->h[0] - rms->h[1] - rms->h[2] + rms->h[3], rms->min, rms->max);
    H[3] = clamp (rms->h[0] + rms->h[1] - rms->h[2] - rms->h[3], rms->min, rms->max);
  }
  function_from_corners (rms->h, H);
}

static void rms_update (RMS * rms)
{
  guint i;
  if (rms->m[0][0] == 0.) {
    for (i = 1; i < NM; i++)
      rms->h[i] = 0.;
    rms->h[0] = G_MAXDOUBLE;
    rms->he = 0.;
    rms->cond = G_MAXDOUBLE;
    return;
  }
  else if (rms->m[0][0] >= NM) {
    guint j;
    for (i = 1; i < NM; i++)
      for (j = 0; j < i; j++)
	rms->m[i][j] = rms->m[j][i];
#if GSL  
    double m[NM*NM], v[NM*NM], s[NM];
    gsl_matrix_view gv = gsl_matrix_view_array (v, NM, NM);
    gsl_vector_view gs = gsl_vector_view_array (s, NM);
    for (i = 0; i < NM; i++)
      for (j = 0; j < NM; j++)
	m[i+NM*j] = rms->m[i][j];
    gsl_matrix_view gm = gsl_matrix_view_array (m, NM, NM);
    gsl_linalg_SV_decomp_jacobi (&gm.matrix, &gv.matrix, &gs.vector);
    rms->cond = s[NM - 1] > 0. ? s[0]/s[NM - 1] : G_MAXDOUBLE;
    if (rms->cond < 10000.) {
      gsl_vector_view gH = gsl_vector_view_array (rms->H, NM);
      gsl_vector_view gh = gsl_vector_view_array (rms->h, NM);
      gsl_linalg_SV_solve (&gm.matrix, &gv.matrix, &gs.vector, &gH.vector, &gh.vector);
      variance_check (rms);
      rms->he = rms_minimum (rms);
      return;
    }
#else
    gdouble ** m = gfs_matrix_new (NM, NM, sizeof (gdouble));
    for (i = 0; i < NM; i++)
      for (j = 0; j < NM; j++)
	m[i][j] = rms->m[i][j];
    if (gfs_matrix_inverse (m, NM, 1e-5)) {
      for (i = 0; i < NM; i++) {
	rms->h[i] = 0.;
	for (j = 0; j < NM; j++)
	  rms->h[i] += m[i][j]*rms->H[j];
      }
      gfs_matrix_free (m);
      variance_check (rms);
      rms->he = rms_minimum (rms);
      return;
    }
    gfs_matrix_free (m);
#endif
  }
  rms->h[0] = rms->H[0]/rms->m[0][0];
  for (i = 1; i < NM; i++)
    rms->h[i] = 0.;
  rms->he = rms_minimum (rms);
}

#if DEBUG
static gdouble rms_value (RMS * rms, FttVector * p)
{
  return rms->h[0] + rms->h[1]*p->x + rms->h[2]*p->y + rms->h[3]*p->x*p->y;
}

static void rms_write (RMS * rms, Polygon * p)
{
  FttVector q, r;
  q.x = p->c.x + p->h; q.y = p->c.y + p->h;
  r.x = 1.; r.y = 1.;
  fprintf (stderr, "%g %g %g\n", q.x, q.y, rms_value (rms, &r)/20000.);
  q.x = p->c.x + p->h; q.y = p->c.y - p->h;
  r.x = 1.; r.y = - 1.;
  fprintf (stderr, "%g %g %g\n", q.x, q.y, rms_value (rms, &r)/20000.);
  q.x = p->c.x - p->h; q.y = p->c.y - p->h;
  r.x = - 1.; r.y = - 1.;
  fprintf (stderr, "%g %g %g\n", q.x, q.y, rms_value (rms, &r)/20000.);
  q.x = p->c.x - p->h; q.y = p->c.y + p->h;
  r.x = - 1.; r.y = 1.;
  fprintf (stderr, "%g %g %g\n", q.x, q.y, rms_value (rms, &r)/20000.);
}

static int write_points (double p[3], RMS * rms)
{
  if (polygon_contains (rms->p, p))
    fprintf (stderr, "aa %g %g %g\n", p[0], p[1], p[2]);
}
#endif

#define RAW        0. /* fitted but not continuous (C0) */
#define FAIR       1. /* fitted and C0 */
#define REFINED    2. /* non-fitted refined cell */
#define NEW_CHILD  3. /* non-fitted child extrapolated from its parent */
#define CONTAINS_SURFACE 4. /* 3D-only */
#define BOUNDARY 5.         /* 3D-only */

static int rms_add (double p[3], RMS * rms)
{
  if (polygon_contains (rms->p, p)) {
    gfs_simulation_map (gfs_object_simulation (rms->p->t), (FttVector *) p);
    if (p[2] > rms->max) rms->max = p[2];
    if (p[2] < rms->min) rms->min = p[2];
    if (rms->relative)
      p[2] -= cell_value (ftt_cell_parent (rms->p->cell), rms->p->t->h, *((FttVector *) p));
    p[0] = (p[0] - rms->p->c.x)/rms->p->h; p[1] = (p[1] - rms->p->c.y)/rms->p->h;
    rms->m[0][1] += p[0];      rms->m[0][2] += p[1];      rms->m[0][3] += p[0]*p[1];
    rms->m[1][1] += p[0]*p[0]; rms->m[1][2] += p[0]*p[1]; rms->m[1][3] += p[0]*p[0]*p[1];
                               rms->m[2][2] += p[1]*p[1]; rms->m[2][3] += p[0]*p[1]*p[1];
                                                          rms->m[3][3] += p[0]*p[0]*p[1]*p[1];
    rms->H[0] += p[2];
    rms->H[1] += p[0]*p[2];
    rms->H[2] += p[1]*p[2];
    rms->H[3] += p[0]*p[1]*p[2];
    rms->H[4] += p[2]*p[2];
    rms->m[0][0] += 1.;
  }
  return 0;
}

static void update_terrain_rms (FttCell * cell, GfsRefineTerrain * t, gboolean relative,
				RMS * rms)
{
  Polygon poly;
  polygon_init (&poly, cell, t);
  rms_init (rms, &poly, relative);
  guint i;
  for (i = 0; i < t->nrs; i++)
    r_surface_query_region (t->rs[i], poly.min, poly.max, (RSurfaceQuery) rms_add, rms);
  rms->p = NULL;
}

static void update_terrain_rms1 (FttCell * cell, GfsRefineTerrain * t, gboolean relative,
				 RMS * rms)
{
  Polygon poly;
  polygon_init (&poly, cell, t);
  rms_init (rms, &poly, relative);
  RSurfaceSum s;
  r_surface_sum_init (&s);
  guint i;
  for (i = 0; i < t->nrs; i++)
    r_surface_query_region_sum (t->rs[i],
				(RSurfaceCheck) polygon_includes,
				(RSurfaceCheck) polygon_intersects,
				&poly, &s);
  rms->m[0][0] = s.n;
  if (s.n > 0) {
    gdouble xc = poly.c.x, yc = poly.c.y;
    gdouble h = poly.h;
    /* See terrain.mac for a "maxima" derivation of the terms below */
    rms->m[0][1] = s.m01 - xc*s.n;
    rms->m[0][1] /= h;
    rms->m[0][2] = s.m02 - yc*s.n;
    rms->m[0][2] /= h;
    rms->m[0][3] = xc*yc*s.n - s.m01*yc - s.m02*xc + s.m03;
    rms->m[0][3] /= h*h;
    rms->m[1][2] = rms->m[0][3];
    rms->m[1][1] = xc*xc*s.n - 2.*s.m01*xc + s.m11;
    rms->m[1][1] /= h*h;
    rms->m[2][2] = yc*yc*s.n - 2.*s.m02*yc + s.m22;
    rms->m[2][2] /= h*h;
    rms->m[1][3] = - xc*xc*yc*s.n + 2.*s.m01*xc*yc - s.m11*yc + s.m02*xc*xc - 2.*s.m03*xc + s.m13;
    rms->m[1][3] /= h*h*h;
    rms->m[2][3] = - xc*yc*yc*s.n + s.m01*yc*yc + 2.*s.m02*xc*yc - 2.*s.m03*yc - s.m22*xc + s.m23;
    rms->m[2][3] /= h*h*h;
    rms->m[3][3] = xc*xc*yc*yc*s.n - 2.*s.m01*xc*yc*yc + s.m11*yc*yc - 2.*s.m02*xc*xc*yc 
      + 4.*s.m03*xc*yc - 2.*s.m13*yc + s.m22*xc*xc - 2.*s.m23*xc + s.m33;
    rms->m[3][3] /= h*h*h*h;  
    if (rms->relative) {
      double hp[NM];
      cell_coefficients (ftt_cell_parent (rms->cell), rms->t->h, hp);
      rms->H[0] = s.H0 - s.n*hp[0] - s.m01*hp[1] - s.m02*hp[2] - s.m03*hp[3];
      /* See terrain.mac for a "maxima" derivation of the terms below */
      rms->H[1] = (s.H1 - xc*s.H0 
		   + s.m03*(hp[3]*xc - hp[2])
		   + s.m01*(hp[1]*xc - hp[0])
		   + hp[2]*s.m02*xc
		   + hp[0]*xc*s.n
		   - hp[3]*s.m13
		   - hp[1]*s.m11);
      rms->H[2] = (s.H2 - yc*s.H0 
		   + s.m03*(hp[3]*yc - hp[1])
		   + s.m02*(hp[2]*yc - hp[0])
		   + hp[1]*s.m01*yc
		   + hp[0]*yc*s.n
		   - hp[3]*s.m23
		   - hp[2]*s.m22);
      rms->H[3] = (s.H3 - xc*s.H2 - yc*s.H1 + xc*yc*s.H0 
		   - s.m03*(hp[3]*xc*yc - hp[2]*yc - hp[1]*xc + hp[0]) 
		   + s.m13*(hp[3]*yc - hp[1]) 
		   - s.m02*xc*(hp[2]*yc - hp[0]) 
		   - s.m01*(hp[1]* xc - hp[0])*yc 
		   - hp[0]*xc*yc*s.n
		   + hp[1]*s.m11*yc 
		   + s.m23*(hp[3]*xc - hp[2]) 
		   + hp[2]*s.m22*xc 
		   - hp[3]*s.m33);
      rms->H[4] = (s.H4 - 2.*hp[3]*s.H3 - 2.*hp[2]*s.H2 - 2.*hp[1]*s.H1 - 2.*hp[0]*s.H0
		   + hp[3]*hp[3]*s.m33
		   + 2.*hp[2]*hp[3]*s.m23
		   + hp[2]*hp[2]*s.m22
		   + 2.*hp[1]*hp[3]*s.m13
		   + hp[1]*hp[1]*s.m11
		   + 2.*(hp[0]*hp[3] + hp[1]*hp[2])*s.m03
		   + 2.*hp[0]*hp[2]*s.m02
		   + 2.*hp[0]*hp[1]*s.m01
		   + hp[0]*hp[0]*s.n);
    }
    else {
      rms->H[0] = s.H0;
      rms->H[1] = s.H1 - xc*s.H0;    
      rms->H[2] = s.H2 - yc*s.H0;
      rms->H[3] = s.H3- xc*s.H2 - yc*s.H1 + xc*yc*s.H0;
      rms->H[4] = s.H4;
    }
    rms->H[1] /= h;
    rms->H[2] /= h;
    rms->H[3] /= h*h;
    rms->max = s.Hmax;
    rms->min = s.Hmin;
  }
  rms->p = NULL;
}

static void update_terrain (FttCell * cell, GfsRefineTerrain * t)
{
  RMS rms;
  guint i;
  g_assert (GFS_VALUE (cell, t->type) == REFINED);
  update_terrain_rms (cell, t, ftt_cell_parent (cell) != NULL, &rms);
#if 1
  RMS rms1;
  update_terrain_rms1 (cell, t, ftt_cell_parent (cell) != NULL, &rms1);
  {
    guint i, j;
    for (i = 0; i < NM; i++)
      for (j = 0; j <= i; j++)
	fprintf (stderr, "m%d%d: %g %g\n", j, i, rms.m[j][i], rms1.m[j][i]);
  }
  {
    guint i;
    for (i = 0; i < NM + 1; i++)
      fprintf (stderr, "H%s%d: %g %g\n", rms.relative ? "*" : "", i, rms.H[i], rms1.H[i]);
  }
  if (rms1.m[0][0] > 0.) {
    fprintf (stderr, "min: %g %g\n", rms.min, rms1.min);
    fprintf (stderr, "max: %g %g\n", rms.max, rms1.max);
  }
  rms_update (&rms1);
#endif
  rms_update (&rms);
  for (i = 0; i < NM; i++) {
    GFS_VALUE (cell, t->h[i]) = rms.h[i];
    fprintf (stderr, "h%d: %g %g\n", i, rms.h[i], rms1.h[i]);
  }
  GFS_VALUE (cell, t->he) = rms.he;
  GFS_VALUE (cell, t->hn) = rms.m[0][0];
  GFS_VALUE (cell, t->type) = RAW;
}

static void function_from_parent (FttCell * cell, GfsRefineTerrain * t, gdouble h[4])
{
  gdouble H[4];
  corners_from_parent (cell, t, H);
  function_from_corners (h, H);
}

static void cell_fine_init (FttCell * parent, GfsRefineTerrain * t)
{
  gfs_cell_fine_init (parent, GFS_DOMAIN (gfs_object_simulation (t)));
  FttCellChildren child;
  ftt_cell_children (parent, &child);
  guint i;
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
      gdouble h[NM];
      function_from_parent (child.c[i], t, h);
      guint j;
      for (j = 0; j < NM; j++)
	GFS_VALUE (child.c[i], t->h[j]) = h[j];
      GFS_VALUE (child.c[i], t->he) = GFS_VALUE (parent, t->he);
      GFS_VALUE (child.c[i], t->hn) = GFS_VALUE (parent, t->hn)/FTT_CELLS;
      GFS_VALUE (child.c[i], t->type) = NEW_CHILD;
    }
}

static gdouble corner_value (GfsRefineTerrain * t, FttVector * p, gdouble eps, guint level)
{
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (t));
  gdouble v = 0., w = 0.;
  gint i, j;
  for (i = -1; i <= 1; i += 2)
    for (j = -1; j <= 1; j += 2) {
      FttVector q;
      q.x = p->x + eps*i; q.y = p->y + eps*j; q.z = p->z;
      FttCell * cell = gfs_domain_locate (domain, q, level);
      if (cell) {
	if (ftt_cell_level (cell) < level)
	    return 0.;
	else if (GFS_VALUE (cell, t->type) == FAIR)
	  return cell_value (cell, t->h, *p);
	gdouble n = GFS_VALUE (cell, t->hn);
	if (n > 0.) {
	  g_assert (GFS_VALUE (cell, t->type) == RAW);
	  v += cell_value (cell, t->h, *p);
	  w += 1.;
	}
      }
    }
  return w > 0 ? v/w : 0.;
}

static void update_error_estimate (FttCell * cell, GfsRefineTerrain * t, gboolean relative)
{
  if (GFS_VALUE (cell, t->hn) > 0.) {
    RMS rms;
    guint i;
    update_terrain_rms (cell, t, relative, &rms);
    for (i = 0; i < NM; i++)
      rms.h[i] = GFS_VALUE (cell, t->h[i]);
#if 1
    RMS rms1;
    update_terrain_rms1 (cell, t, relative, &rms1);
    for (i = 0; i < NM; i++)
      rms1.h[i] = GFS_VALUE (cell, t->h[i]);
    fprintf (stderr, "he: %g %g\n", rms_minimum (&rms), rms_minimum (&rms1));
#endif
    GFS_VALUE (cell, t->he) = rms_minimum (&rms);
  }
  else
    GFS_VALUE (cell, t->he) = 0.;
}

static void remove_knots (FttCell * cell, GfsRefineTerrain * t)
{
  gdouble size = ftt_cell_size (cell), eps = size/1000.;
  guint level = ftt_cell_level (cell);
  FttVector p;
  ftt_cell_pos (cell, &p);
  gdouble h[4], H[4];
  p.x += size/2.; p.y += size/2.;
  H[0] = corner_value (t, &p, eps, level);
  p.x -= size;
  H[1] = corner_value (t, &p, eps, level);
  p.y -= size;
  H[2] = corner_value (t, &p, eps, level);
  p.x += size;
  H[3] = corner_value (t, &p, eps, level);
  function_from_corners (h, H);
  guint i;
  for (i = 0; i < NM; i++)
    GFS_VALUE (cell, t->h[i]) = h[i];
  GFS_VALUE (cell, t->type) = FAIR;

  update_error_estimate (cell, t, ftt_cell_parent (cell) != NULL);
}

static void update_height_and_check_for_refinement (FttCell * cell, GfsRefineTerrain * t)
{
  if (GFS_VALUE (cell, t->type) == FAIR) {
    if (ftt_cell_parent (cell)) {
      gdouble h[4];
      function_from_parent (cell, t, h);
      guint i;
      for (i = 0; i < NM; i++)
	GFS_VALUE (cell, t->h[i]) += h[i];
    }

    if (ftt_cell_level (cell) < gfs_function_value (GFS_REFINE (t)->maxlevel, cell) &&
	gfs_function_value (t->criterion, cell)) {
      g_assert (FTT_CELL_IS_LEAF (cell));
      ftt_cell_refine_single (cell, (FttCellInitFunc) cell_fine_init, t);
      FttCellChildren child;
      guint i;
      ftt_cell_children (cell, &child);
      for (i = 0; i < FTT_CELLS; i++)
	GFS_VALUE (child.c[i], t->type) = REFINED;
      t->refined = TRUE;
    }
  }
  else
    g_assert (GFS_VALUE (cell, t->type) == NEW_CHILD);
}

static void reset_terrain (FttCell * cell, GfsRefineTerrain * t)
{
  guint i;
  for (i = 0; i < NM; i++)
    GFS_VALUE (cell, t->h[i]) = 0.;
  GFS_VALUE (cell, t->type) = REFINED;
  if (FTT_CELL_IS_LEAF (cell) && ftt_cell_level (cell) < t->level)
    t->level = ftt_cell_level (cell);
}

#if FTT_2D
# define traverse_boundary(domain,order,flags,depth,func,data) \
         gfs_domain_cell_traverse(domain,order,flags,depth,func,data)
#else /* 3D */
# define traverse_boundary(domain,order,flags,depth,func,data) \
         gfs_domain_cell_traverse_boundary(domain,FTT_FRONT,order,flags,depth,func,data)

static void terrain_min_max (FttCell * cell, GfsVariable * h[NM], gdouble minmax[2])
{
  gdouble dx, dy;
  gdouble H0 = GFS_VALUE (cell, h[0]), H1 = GFS_VALUE (cell, h[1]);
  gdouble H2 = GFS_VALUE (cell, h[2]), H3 = GFS_VALUE (cell, h[3]);
  minmax[0] = G_MAXDOUBLE; minmax[1] = - G_MAXDOUBLE;
  for (dx = -1.; dx <= 1.; dx += 2.)
    for (dy = -1.; dy <= 1.; dy += 2.) {
      gdouble v = H0 + dx*H1 + dy*H2 + dx*dy*H3;
      if (v < minmax[0]) minmax[0] = v;
      if (v > minmax[1]) minmax[1] = v;
    }
}

static void min_max (FttCell * cell, GfsRefineTerrain * t)
{
  gdouble minmax[2] = { G_MAXDOUBLE, - G_MAXDOUBLE };
  if (FTT_CELL_IS_LEAF (cell)) {
    terrain_min_max (cell, t->h, minmax);
    FttVector p;
    ftt_cell_pos (cell, &p);
    if (p.z > t->front)
      t->front = p.z;
  }
  else {
    FttCellChildren child;
    guint i;
    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i]) {
	if (GFS_VALUE (child.c[i], t->max) > minmax[1])
	  minmax[1] = GFS_VALUE (child.c[i], t->max);
	if (GFS_VALUE (child.c[i], t->min) < minmax[0])
	  minmax[0] = GFS_VALUE (child.c[i], t->min);	
      }
  }
  GFS_VALUE (cell, t->min) = minmax[0];
  GFS_VALUE (cell, t->max) = minmax[1];
  GFS_VALUE (cell, t->type) = BOUNDARY;
}

static gboolean refine_terrain_from_boundary (FttCell * cell, GfsRefineTerrain * t)
{
  FttVector p;
  ftt_cell_pos (cell, &p);
  gdouble h = ftt_cell_size (cell)/2., zmin = p.z - h, zmax = p.z + h;
  p.z = t->front;
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (t));
  FttCell * boundary = gfs_domain_locate (domain, p, ftt_cell_level (cell));
  g_assert (boundary);
  if (GFS_VALUE (boundary, t->min) > zmax || GFS_VALUE (boundary, t->max) < zmin)
    return FALSE;
  GFS_VALUE (cell, t->type) = CONTAINS_SURFACE;
  return !FTT_CELL_IS_LEAF (boundary);
}

static void refine_box (GfsBox * box, GfsRefineTerrain * t)
{
  ftt_cell_refine (box->root, 
		   (FttCellRefineFunc) refine_terrain_from_boundary, t,
		   (FttCellInitFunc) gfs_cell_fine_init, gfs_box_domain (box));
}

static void init_terrain_from_boundary (FttCell * cell, GfsRefineTerrain * t)
{
  if (GFS_VALUE (cell, t->type) == CONTAINS_SURFACE) {
    FttVector p;
    ftt_cell_pos (cell, &p);
    p.z = t->front;
    GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (t));
    FttCell * boundary = gfs_domain_locate (domain, p, -1);
    g_assert (boundary);
    g_assert (ftt_cell_level (cell) == ftt_cell_level (boundary));
    guint i;
    for (i = 0; i < NM; i++)
      GFS_VALUE (cell, t->h[i]) = GFS_VALUE (boundary, t->h[i]);
    GFS_VALUE (cell, t->he) = GFS_VALUE (boundary, t->he);
    GFS_VALUE (cell, t->hn) = GFS_VALUE (boundary, t->hn);
  }
}

static gboolean coarsen_boundary (FttCell * cell, GfsRefineTerrain * t)
{
  return (GFS_VALUE (cell, t->type) != CONTAINS_SURFACE);
}

static void coarsen_box (GfsBox * box, GfsRefineTerrain * t)
{
  ftt_cell_coarsen (box->root,
		    (FttCellCoarsenFunc) coarsen_boundary, t,
		    (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
}

static void reset_empty_cell (FttCell * cell, GfsRefineTerrain * t)
{
  if (GFS_VALUE (cell, t->type) != CONTAINS_SURFACE) {
    guint i;
    for (i = 0; i < NM; i++)
      GFS_VALUE (cell, t->h[i]) = G_MAXDOUBLE;
    GFS_VALUE (cell, t->he) = G_MAXDOUBLE;
    GFS_VALUE (cell, t->hn) = G_MAXDOUBLE;
  }
}
#endif /* 3D */

#if DEBUG
static void draw_terrain (FttCell * cell, gpointer * data)
{
  GfsRefineTerrain * t = data[0];
  FILE * fp = data[1];
  gdouble h = ftt_cell_size (cell);
  FttVector p;
  ftt_cell_pos (cell, &p);
  p.x += h/2.; p.y += h/2.;
  fprintf (fp, "%g %g %g\n", p.x, p.y, cell_value (cell, t->h, p));
  p.x -= h;
  fprintf (fp, "%g %g %g\n", p.x, p.y, cell_value (cell, t->h, p));
  p.y -= h;
  fprintf (fp, "%g %g %g\n", p.x, p.y, cell_value (cell, t->h, p));
  p.x += h;
  fprintf (fp, "%g %g %g\n", p.x, p.y, cell_value (cell, t->h, p));
  p.y += h;
  fprintf (fp, "%g %g %g\n\n\n", p.x, p.y, cell_value (cell, t->h, p));
}

static void draw_level (GfsDomain * domain, GfsRefine * refine, guint level, const gchar * name)
{
  gpointer data[2];
  data[0] = refine;
  data[1] = fopen (name, "w");
  fprintf (data[1], "QUAD\n");
  traverse_boundary (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, level,
		     (FttCellTraverseFunc) draw_terrain, data);
  fclose (data[1]);
}

static void draw_all (GfsDomain * domain, GfsRefine * refine, const gchar * name)
{
  gpointer data[2];
  data[0] = refine;
  data[1] = fopen (name, "w");
  //  fprintf (data[1], "QUAD\n");
  traverse_boundary (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) draw_terrain, data);
  fclose (data[1]);
}
#endif

#define ASCII_ZERO 48 /* ASCII value for character "0" */

static void terrain_coarse_fine (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  guint len = strlen (v->name) - 1;
  gint c = v->name[len] - ASCII_ZERO;
  guint n;
  gdouble h[NM];

  g_assert (c >= 0 && c < NM);
  for (n = 0; n < NM; n++) {
    GSList * i = v->domain->variables;
    while (i && (!GFS_VARIABLE1 (i->data)->name || 
		 strncmp (v->name, GFS_VARIABLE1 (i->data)->name, len) ||
		 GFS_VARIABLE1 (i->data)->name[len] != ASCII_ZERO + n))
      i = i->next;
    g_assert (i);
    h[n] = GFS_VALUE (parent, GFS_VARIABLE1 (i->data));
  }

  ftt_cell_children (parent, &child);
  if (h[0] == G_MAXDOUBLE) {
    for (n = 0; n < FTT_CELLS; n++)
      if (child.c[n])
	GFS_VALUE (child.c[n], v) = G_MAXDOUBLE;
  }
  else {
    gdouble size = ftt_cell_size (parent)/4.;
    for (n = 0; n < FTT_CELLS; n++)
      if (child.c[n]) {
	gdouble hc[NM];
	FttVector p;
	ftt_cell_relative_pos (child.c[n], &p);
	p.x *= 2.; p.y *= 2.;
	hc[0] = h[0] + h[1]*p.x + h[2]*p.y + h[3]*p.x*p.y;
	hc[1] = (h[1] + h[3]*p.y)/2.;
	hc[2] = (h[2] + h[3]*p.x)/2.;
	hc[3] = h[3]/4.;
#if !FTT_2D
	ftt_cell_pos (child.c[n], &p);
	gdouble zmin = p.z - size, zmax = p.z + size;
	gdouble dx, dy;
	gdouble minmax[2] = { G_MAXDOUBLE, - G_MAXDOUBLE };
	for (dx = -1.; dx <= 1.; dx += 2.)
	  for (dy = -1.; dy <= 1.; dy += 2.) {
	    gdouble v = hc[0] + dx*hc[1] + dy*hc[2] + dx*dy*hc[3];
	    if (v < minmax[0]) minmax[0] = v;
	    if (v > minmax[1]) minmax[1] = v;
	  }
	if (minmax[0] > zmax || minmax[1] < zmin)
	  GFS_VALUE (child.c[n], v) = G_MAXDOUBLE;
	else
#endif
	  GFS_VALUE (child.c[n], v) = hc[c];
      }
  }
}

static void terrain_refine (GfsRefine * refine, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (refine);
  t->type = gfs_temporary_variable (domain);
  t->level = G_MAXINT/2;
  traverse_boundary (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		     (FttCellTraverseFunc) reset_terrain, refine);
  do {
    traverse_boundary (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, t->level,
		       (FttCellTraverseFunc) update_terrain, refine);
    traverse_boundary (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, t->level,
		       (FttCellTraverseFunc) remove_knots, refine);
    t->refined = FALSE;
    traverse_boundary (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, t->level,
		       (FttCellTraverseFunc) update_height_and_check_for_refinement,
		       refine);
#if DEBUG
    GfsNorm norm = gfs_domain_norm_variable (domain, t->he, NULL, FTT_TRAVERSE_LEAFS, -1);
    fprintf (stderr, "level: %d bias: %g 1: %g 2: %g inf: %g\n", 
	     t->level, norm.bias, norm.first, norm.second, norm.infty);
    fprintf (stderr, "level: %d depth: %d\n", t->level, gfs_domain_depth (domain));
    gchar name[] = "/tmp/level-x";
    name[11] = ASCII_ZERO + t->level;
    draw_level (domain, refine, t->level, name);
#endif
    t->level++;
  } while (t->refined);
#if DEBUG
  draw_all (domain, refine, "/tmp/all");
#endif
#if !FTT_2D
  /* The height field is only defined on the front boundary, we need
     to define it volumetrically */
  t->min = gfs_temporary_variable (domain);
  t->max = gfs_temporary_variable (domain);
  t->front = - G_MAXDOUBLE;
  gfs_domain_cell_traverse_boundary (domain, FTT_FRONT, FTT_POST_ORDER, FTT_TRAVERSE_ALL, -1,
				     (FttCellTraverseFunc) min_max, t);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) refine_box, t);
  gts_object_destroy (GTS_OBJECT (t->min));
  gts_object_destroy (GTS_OBJECT (t->max));
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) init_terrain_from_boundary, t);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) coarsen_box, t);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) reset_empty_cell, t);
#endif /* 3D */  
  gts_object_destroy (GTS_OBJECT (t->type));
  guint i;
  for (i = 0; i < NM; i++)
    t->h[i]->coarse_fine = terrain_coarse_fine;
}

static void refine_terrain_destroy (GtsObject * object)
{
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (object);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (object));

  gchar * dname = g_strconcat (t->name, "min", NULL);
  gfs_domain_remove_derived_variable (domain, dname);
  g_free (dname);

  dname = g_strconcat (t->name, "max", NULL);
  gfs_domain_remove_derived_variable (domain, dname);
  g_free (dname);
  
  g_free (t->name);
  g_free (t->basename);
  if (t->rs) {
    guint i;
    for (i = 0; i < t->nrs; i++)
      r_surface_close (t->rs[i]);
    g_free (t->rs);
  }
  gts_object_destroy (GTS_OBJECT (t->criterion));  
  (* GTS_OBJECT_CLASS (gfs_refine_terrain_class ())->parent_class->destroy) (object);
}

static gdouble terrain_hmin (FttCell * cell, FttCellFace * face, 
			     GfsDomain * domain, GfsRefineTerrain * t)
{
  g_return_val_if_fail (cell != NULL, 0.);

  gdouble dx, dy, min = G_MAXDOUBLE;
  gdouble H0 = GFS_VALUE (cell, t->h[0]), H1 = GFS_VALUE (cell, t->h[1]);
  gdouble H2 = GFS_VALUE (cell, t->h[2]), H3 = GFS_VALUE (cell, t->h[3]);

  for (dx = -1.; dx <= 1.; dx += 2.)
    for (dy = -1.; dy <= 1.; dy += 2.) {
      double v = H0 + dx*H1 + dy*H2 + dx*dy*H3;
      if (v < min) min = v;
    }
  return min;
}

static gdouble terrain_hmax (FttCell * cell, FttCellFace * face, 
			     GfsDomain * domain, GfsRefineTerrain * t)
{
  g_return_val_if_fail (cell != NULL, 0.);

  gdouble dx, dy, max = - G_MAXDOUBLE;
  gdouble H0 = GFS_VALUE (cell, t->h[0]), H1 = GFS_VALUE (cell, t->h[1]);
  gdouble H2 = GFS_VALUE (cell, t->h[2]), H3 = GFS_VALUE (cell, t->h[3]);

  for (dx = -1.; dx <= 1.; dx += 2.)
    for (dy = -1.; dy <= 1.; dy += 2.) {
      double v = H0 + dx*H1 + dy*H2 + dx*dy*H3;
      if (v > max) max = v;
    }
  return max;
}

static void refine_terrain_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_refine_terrain_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (name)");
    return;
  }
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (*o);
  t->name = g_strdup (fp->token->str);
  gts_file_next_token (fp);

  gchar * path = NULL;
  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_STRING, "basename", TRUE},
      {GTS_STRING, "path",      TRUE},
      {GTS_NONE}
    };
    gchar * basename = NULL;
    var[0].data = &basename;
    var[1].data = &path;
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
    if (var[0].set) { g_free (t->basename); t->basename = basename; }
    if (!var[1].set) path = g_strdup (default_path);
  }
  else
    path = g_strdup (default_path);

  if (!strcmp (t->basename, "*")) { /* file globbing */
    gchar * pattern = g_strconcat (path, "/*.Data", NULL);
    glob_t pglob;
    if (glob (pattern, GLOB_ERR, NULL, &pglob)) {
      gts_file_error (fp, "cannot find/open terrain databases in path:\n%s", pattern);
      g_free (pattern);
      g_free (path);
      return;
    }
    g_free (pattern);
    g_free (t->basename);
    t->basename = NULL;
    guint i;
    for (i = 0; i < pglob.gl_pathc; i++) {
      pglob.gl_pathv[i][strlen (pglob.gl_pathv[i]) - 5] = '\0';
      t->rs = g_realloc (t->rs, (t->nrs + 1)*sizeof (RSurface *));
      t->rs[t->nrs] = r_surface_open (pglob.gl_pathv[i], "r", -1);
      if (!t->rs[t->nrs]) {
	gts_file_error (fp, "cannot open terrain database `%s'", pglob.gl_pathv[i]);
	globfree (&pglob);
	g_free (path);
	return;
      }
      if (t->basename) {
	gchar * pathbasename = g_strconcat (t->basename, ",", pglob.gl_pathv[i], NULL);
	g_free (t->basename);
	t->basename = pathbasename;
      }
      else
	t->basename = g_strdup (pglob.gl_pathv[i]);
      t->nrs++;
    }
    globfree (&pglob);
  }
  else { /* basename is of the form: set1,set2,set3... */
    gchar ** names = g_strsplit (t->basename, ",", 0);
    if (path) {
      g_free (t->basename);
      t->basename = NULL;
    }
    gchar ** s = names;
    while (*s) {
      t->rs = g_realloc (t->rs, (t->nrs + 1)*sizeof (RSurface *));
      if (path) {
	/* search path */
	gchar ** pathes = g_strsplit (path, ":", 0);
	gchar ** spath = pathes, * fname;
	g_assert (*spath);
	do {
	  fname = (*s)[0] == '/' ? g_strdup (*s) : g_strconcat (*spath, "/", *s, NULL);
	  t->rs[t->nrs] = r_surface_open (fname, "r", -1);
	} while (t->rs[t->nrs] == NULL && *(++spath));
	g_strfreev (pathes);

	if (t->basename) {
	  gchar * pathbasename = g_strconcat (t->basename, ",", fname, NULL);
	  g_free (t->basename);
	  t->basename = pathbasename;
	  g_free (fname);
	}
	else
	  t->basename = fname;
      }
      else
	t->rs[t->nrs] = r_surface_open (*s, "r", -1);
      if (!t->rs[t->nrs]) {
	if (path)
	  gts_file_error (fp, "cannot find/open terrain database `%s' in path:\n%s", *s, path);
	else
	  gts_file_error (fp, "cannot open terrain database `%s'", *s);
	g_strfreev (names);
	g_free (path);
	return;
      }
      t->nrs++;
      s++;
    }
    g_strfreev (names);
  }
  g_free (path);

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  guint i;
  for (i = 0; i < NM; i++) {
    gchar * name = g_strdup_printf ("%s%d", t->name, i);
    t->h[i] = gfs_domain_get_or_add_variable (domain, name, "Terrain height");
    g_free (name);
  }
  gchar * name = g_strjoin (NULL, t->name, "e", NULL);
  t->he = gfs_domain_get_or_add_variable (domain, name, "Terrain RMS error");
  g_free (name);
  name = g_strjoin (NULL, t->name, "n", NULL);
  t->hn = gfs_domain_get_or_add_variable (domain, name, "Terrain samples #");
  g_free (name);

  GfsDerivedVariableInfo v;

  v.name = g_strjoin (NULL, t->name, "min", NULL);
  v.description = "Minimum terrain height";
  v.func = terrain_hmin;
  v.data = t;
  if (!gfs_domain_add_derived_variable (domain, v)) {
    gts_file_error (fp, "derived variable `%s' already defined", v.name);
    g_free (v.name);
    return;
  }
  g_free (v.name);

  v.name = g_strjoin (NULL, t->name, "max", NULL);
  v.description = "Maximum terrain height";
  v.func = terrain_hmax;
  v.data = t;
  if (!gfs_domain_add_derived_variable (domain, v)) {
    gts_file_error (fp, "derived variable `%s' already defined", v.name);
    g_free (v.name);
    return;
  }
  g_free (v.name);

  gfs_function_read (t->criterion, domain, fp);
}

static void refine_terrain_write (GtsObject * o, FILE * fp)
{
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (o);
  (* GTS_OBJECT_CLASS (gfs_refine_terrain_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s ", t->name);
  gfs_function_write (t->criterion, fp);
  fprintf (fp, " { basename = %s }", t->basename);
}

static void gfs_refine_terrain_class_init (GfsRefineClass * klass)
{
  klass->refine = terrain_refine;

  GTS_OBJECT_CLASS (klass)->destroy = refine_terrain_destroy;
  GTS_OBJECT_CLASS (klass)->read = refine_terrain_read;
  GTS_OBJECT_CLASS (klass)->write = refine_terrain_write;
}

static void gfs_refine_terrain_init (GfsRefineTerrain * t)
{
  t->criterion = gfs_function_new (gfs_function_class (), 0.);
  t->basename = g_strdup ("*");
}

GfsRefineClass * gfs_refine_terrain_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_terrain_info = {
      "GfsRefineTerrain",
      sizeof (GfsRefineTerrain),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_terrain_class_init,
      (GtsObjectInitFunc) gfs_refine_terrain_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_class ()),
				  &gfs_refine_terrain_info);
  }

  return klass;
}

/* GfsSurfaceTerrain: Header */

typedef struct _GfsSurfaceTerrain         GfsSurfaceTerrain;

struct _GfsSurfaceTerrain {
  /*< private >*/
  GfsGenericSurface parent;
  GfsVariable * h[NM];

  /*< public >*/
  gchar * name;
};

#define GFS_SURFACE_TERRAIN(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSurfaceTerrain,\
					         gfs_surface_terrain_class ())
#define GFS_IS_SURFACE_TERRAIN(obj)         (gts_object_is_from_class (obj,\
						 gfs_surface_terrain_class ()))

GfsGenericSurfaceClass * gfs_surface_terrain_class  (void);

/* GfsSurfaceTerrain: Object */

static void gfs_surface_terrain_read (GtsObject ** o, GtsFile * fp)
{
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable name");
    return;
  }
  GfsSurfaceTerrain * t = GFS_SURFACE_TERRAIN (*o);
  t->name = g_strdup (fp->token->str);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  guint i;
  for (i = 0; i < NM; i++) {
    gchar * name = g_strdup_printf ("%s%d", t->name, i);
    t->h[i] = gfs_variable_from_name (domain->variables, name);
    if (!t->h[i]) {
      gts_file_error (fp, "%s is not a valid variable name", name);
      g_free (name);
      return;
    }
    t->h[i]->coarse_fine = terrain_coarse_fine;
    g_free (name);
  }
  gts_file_next_token (fp);
}

static void gfs_surface_terrain_write (GtsObject * o, FILE * fp)
{
  fprintf (fp, " %s", GFS_SURFACE_TERRAIN (o)->name);
}

static void gfs_surface_terrain_destroy (GtsObject * object)
{
  g_free (GFS_SURFACE_TERRAIN (object)->name);
  (* GTS_OBJECT_CLASS (gfs_surface_terrain_class ())->parent_class->destroy)
    (object);
}

static GfsGenericSurface * cell_is_cut (FttCell * cell, GfsGenericSurface * s1,
					gboolean flatten, gint maxlevel)
{
  g_assert (!flatten); /* not implemented */
  if (!FTT_CELL_IS_LEAF (cell))
    return s1;
  return GFS_VALUE (cell, GFS_SURFACE_TERRAIN (s1)->h[0]) != G_MAXDOUBLE ? s1 : NULL;
}

static guint surface_segment_intersection (GfsGenericSurface * s1,
					   FttCell * cell,
					   GfsSegment * I)
{
  I->n = 0;
  I->x = 0.;
  I->inside = 0;

  FttVector pE, pD;
  pE.x = I->E->x; pE.y = I->E->y;
  pD.x = I->D->x; pD.y = I->D->y;
  gdouble vE = I->E->z - cell_value (cell, GFS_SURFACE_TERRAIN (s1)->h, pE);
  gdouble vD = I->D->z - cell_value (cell, GFS_SURFACE_TERRAIN (s1)->h, pD);
  
  if ((vE > 0. && vD <= 0.) || (vE <= 0. && vD > 0.)) {
    I->n = 1;
    I->inside = vE > 0. ? -1 : 1;
    I->x = vE/(vE - vD);
#if DEBUG
    gdouble size = ftt_cell_size (cell)/2.;
    FttVector q;
    ftt_cell_pos (cell, &q);
    pE.x = (pE.x - q.x)/size;
    pE.y = (pE.y - q.y)/size;
    pD.x = (pD.x - q.x)/size;
    pD.y = (pD.y - q.y)/size;
    fprintf (stderr, "p %g %g %g %g %g %g %g %d %g %g %g %g\n", 
	     I->D->x, I->D->y, I->D->z,
	     I->E->x, I->E->y, I->E->z,
	     I->x,
	     ftt_cell_level (cell),
	     pE.x, pE.y, pD.x, pD.y);
    fprintf (stderr, "q %g %g %g\nq %g %g %g\nq\nq\n",
	     I->D->x, I->D->y, I->D->z,
	     I->E->x, I->E->y, I->E->z);
    fprintf (stderr, "i %g %g %g\n",
	     I->E->x + I->x*(I->D->x - I->E->x),
	     I->E->y + I->x*(I->D->y - I->E->y),
	     I->E->z + I->x*(I->D->z - I->E->z));
#endif
  }
  return I->n;
}

static void surface_segment_normal (GfsGenericSurface * s1,
				    FttCell * cell,
				    GfsSegment * I,
				    GtsVector n)
{
  GfsVariable ** h = GFS_SURFACE_TERRAIN (s1)->h;
  gdouble size = ftt_cell_size (cell)/2.;
  FttVector p, q;
  ftt_cell_pos (cell, &q);
  p.x = I->E->x + I->x*(I->D->x - I->E->x);
  p.y = I->E->y + I->x*(I->D->y - I->E->y);
  p.x = (p.x - q.x)/size;
  p.y = (p.y - q.y)/size;
  n[0] = - (GFS_VALUE (cell, h[1]) + GFS_VALUE (cell, h[3])*p.y)/size;
  n[1] = - (GFS_VALUE (cell, h[2]) + GFS_VALUE (cell, h[3])*p.x)/size;
  n[2] = 1.;
}

static void gfs_surface_terrain_class_init (GfsGenericSurfaceClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_surface_terrain_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_surface_terrain_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_surface_terrain_destroy;

  klass->cell_is_cut = cell_is_cut;
  klass->segment_intersection = surface_segment_intersection;
  klass->segment_normal = surface_segment_normal;
}

GfsGenericSurfaceClass * gfs_surface_terrain_class (void)
{
  static GfsGenericSurfaceClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_surface_terrain_info = {
      "GfsSurfaceTerrain",
      sizeof (GfsSurfaceTerrain),
      sizeof (GfsGenericSurfaceClass),
      (GtsObjectClassInitFunc) gfs_surface_terrain_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_surface_class ()),
				  &gfs_surface_terrain_info);
  }

  return klass;
}

/* GfsTerrain: Header */

#define GFS_IS_TERRAIN(obj)         (gts_object_is_from_class (obj,\
						 gfs_terrain_class ()))

GfsEventClass * gfs_terrain_class  (void);

/* GfsTerrain: Object */

static void terrain_init (GfsSolid * s)
{
  gts_object_destroy (GTS_OBJECT (s->s));
  s->s = GFS_GENERIC_SURFACE (gts_object_new (GTS_OBJECT_CLASS (gfs_surface_terrain_class ())));
}

GfsEventClass * gfs_terrain_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_terrain_info = {
      "GfsTerrain",
      sizeof (GfsSolid),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) terrain_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_solid_class ()),
				  &gfs_terrain_info);
  }

  return klass;
}

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "terrain";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gchar * path = getenv ("GFS_TERRAIN_PATH");
  if (path)
    default_path = path;
  gfs_refine_terrain_class ();
  gfs_terrain_class ();
  return NULL;
}
