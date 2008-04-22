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
#include "rsurface.h"

static gchar * default_dir = "/home/popinet/Projects/GIS/rsurfaces";

/* GfsRefineTerrain: Header */

typedef struct _GfsRefineTerrain         GfsRefineTerrain;

#define NM 4

struct _GfsRefineTerrain {
  GfsRefine parent;
  guint level;
  gboolean refined;

  gchar * name, * basename;
  RSurface ** rs;
  guint nrs;

  GfsVariable * h[NM], * he, * hn, * cond;
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

typedef struct {
  gdouble H[NM+1], m[NM][NM];
  gdouble h[NM], he, cond, min, max;
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

static gdouble cell_value (FttCell * cell, GfsRefineTerrain * t, FttVector p)
{
  gdouble h = ftt_cell_size (cell)/2.;
  FttVector q;
  ftt_cell_pos (cell, &q);
  p.x = (p.x - q.x)/h;
  p.y = (p.y - q.y)/h;
  return GFS_VALUE (cell, t->h[0]) + 
    GFS_VALUE (cell, t->h[1])*p.x + 
    GFS_VALUE (cell, t->h[2])*p.y + 
    GFS_VALUE (cell, t->h[3])*p.x*p.y;
}

static void corners_from_parent (FttCell * cell, GfsRefineTerrain * t, gdouble H[4])
{
  gdouble size = ftt_cell_size (cell);
  FttCell * parent = ftt_cell_parent (cell);
  FttVector p;
  ftt_cell_pos (cell, &p);
  p.x += size/2.; p.y += size/2.;
  H[0] = cell_value (parent, t, p);
  p.x -= size;
  H[1] = cell_value (parent, t, p);
  p.y -= size;
  H[2] = cell_value (parent, t, p);
  p.x += size;
  H[3] = cell_value (parent, t, p);
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
    corners_from_parent (rms->p->cell, rms->p->t, H0);
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
    gdouble ** m = gfs_matrix_new (NM, sizeof (gdouble));
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

static int rms_add (double p[3], RMS * rms)
{
  if (polygon_contains (rms->p, p)) {
    gfs_simulation_map (gfs_object_simulation (rms->p->t), (FttVector *) p);
    if (p[2] > rms->max) rms->max = p[2];
    if (p[2] < rms->min) rms->min = p[2];
    if (rms->relative)
      p[2] -= cell_value (ftt_cell_parent (rms->p->cell), rms->p->t, *((FttVector *) p));
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

static void update_terrain (FttCell * cell, GfsRefineTerrain * t)
{
  g_assert (GFS_VALUE (cell, t->cond) == REFINED);
  Polygon p;
  polygon_init (&p, cell, t);
  RMS rms;
  rms_init (&rms, &p, TRUE);
  guint i;
  for (i = 0; i < t->nrs; i++)
    r_surface_query_region (t->rs[i], p.min, p.max, (RSurfaceQuery) rms_add, &rms);
  rms_update (&rms);
  for (i = 0; i < NM; i++)
    GFS_VALUE (cell, t->h[i]) = rms.h[i];
  GFS_VALUE (cell, t->he) = rms.he;
  GFS_VALUE (cell, t->hn) = rms.m[0][0];
  GFS_VALUE (cell, t->cond) = RAW;
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
      GFS_VALUE (child.c[i], t->cond) = NEW_CHILD;
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
      q.x = p->x + eps*i; q.y = p->y + eps*j;
      FttCell * cell = gfs_domain_locate (domain, q, level);
      if (cell) {
	if (ftt_cell_level (cell) < level)
	    return 0.;
	else if (GFS_VALUE (cell, t->cond) == FAIR)
	  return cell_value (cell, t, q);
	gdouble n = GFS_VALUE (cell, t->hn);
	if (n > 0.) {
	  g_assert (GFS_VALUE (cell, t->cond) == RAW);
	  v += cell_value (cell, t, q);
	  w += 1.;
	}
      }
    }
  return w > 0 ? v/w : 0.;
}

static void update_error_estimate (FttCell * cell, GfsRefineTerrain * t, gboolean relative)
{
  if (GFS_VALUE (cell, t->hn) > 0.) {
    Polygon poly;
    polygon_init (&poly, cell, t);
    RMS rms;
    rms_init (&rms, &poly, relative);
    guint i;
    for (i = 0; i < t->nrs; i++)
      r_surface_query_region (t->rs[i], poly.min, poly.max, (RSurfaceQuery) rms_add, &rms);
    for (i = 0; i < NM; i++)
      rms.h[i] = GFS_VALUE (cell, t->h[i]);
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
  GFS_VALUE (cell, t->cond) = FAIR;

  update_error_estimate (cell, t, TRUE);
}

static void update_height_and_check_for_refinement (FttCell * cell, GfsRefineTerrain * t)
{
  if (GFS_VALUE (cell, t->cond) == FAIR) {
    gdouble h[4];
    function_from_parent (cell, t, h);
    guint i;
    for (i = 0; i < NM; i++)
      GFS_VALUE (cell, t->h[i]) += h[i];

    if (ftt_cell_level (cell) < gfs_function_value (GFS_REFINE (t)->maxlevel, cell) &&
	gfs_function_value (t->criterion, cell)) {
      g_assert (FTT_CELL_IS_LEAF (cell));
      ftt_cell_refine_single (cell, (FttCellInitFunc) cell_fine_init, t);
      FttCellChildren child;
      guint i;
      ftt_cell_children (cell, &child);
      for (i = 0; i < FTT_CELLS; i++)
	GFS_VALUE (child.c[i], t->cond) = REFINED;
      t->refined = TRUE;
    }
  }
  else
    g_assert (GFS_VALUE (cell, t->cond) == NEW_CHILD);
}

#if DEBUG
static void draw_terrain (FttCell * cell, gpointer * data)
{
  GfsRefineTerrain * t = data[0];
  FILE * fp = data[1];
  gdouble h = ftt_cell_size (cell);
  FttVector p;
  ftt_cell_pos (cell, &p);
  p.x += h/2.; p.y += h/2.;
  fprintf (fp, "%g %g %g\n", p.x, p.y, cell_value (cell, t, p)/20000.);
  p.x -= h;
  fprintf (fp, "%g %g %g\n", p.x, p.y, cell_value (cell, t, p)/20000.);
  p.y -= h;
  fprintf (fp, "%g %g %g\n", p.x, p.y, cell_value (cell, t, p)/20000.);
  p.x += h;
  fprintf (fp, "%g %g %g\n", p.x, p.y, cell_value (cell, t, p)/20000.);
}

static void draw_level (GfsDomain * domain, GfsRefine * refine, guint level, const gchar * name)
{
  gpointer data[2];
  data[0] = refine;
  data[1] = fopen (name, "w");
  fprintf (data[1], "QUAD\n");
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, level,
			    (FttCellTraverseFunc) draw_terrain, data);
  fclose (data[1]);
}

static void draw_all (GfsDomain * domain, GfsRefine * refine, const gchar * name)
{
  gpointer data[2];
  data[0] = refine;
  data[1] = fopen (name, "w");
  fprintf (data[1], "QUAD\n");
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) draw_terrain, data);
  fclose (data[1]);
}
#endif

static void reset_terrain (FttCell * cell, GfsRefineTerrain * t)
{
  guint i;
  for (i = 0; i < NM; i++)
    GFS_VALUE (cell, t->h[i]) = 0.;
  GFS_VALUE (cell, t->cond) = REFINED;
  if (FTT_CELL_IS_LEAF (cell) && ftt_cell_level (cell) < t->level)
    t->level = ftt_cell_level (cell);
}

static void terrain_refine (GfsRefine * refine, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (refine);
  t->level = G_MAXINT/2;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) reset_terrain, refine);
  do {
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, t->level,
			      (FttCellTraverseFunc) update_terrain, refine);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, t->level,
			      (FttCellTraverseFunc) remove_knots, refine);
    t->refined = FALSE;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, t->level,
			      (FttCellTraverseFunc) update_height_and_check_for_refinement,
			      refine);
#if DEBUG
    GfsNorm norm = gfs_domain_norm_variable (domain, t->he, NULL, FTT_TRAVERSE_LEAFS, -1);
    fprintf (stderr, "level: %d bias: %g 1: %g 2: %g inf: %g\n", 
	     t->level, norm.bias, norm.first, norm.second, norm.infty);
    fprintf (stderr, "level: %d depth: %d\n", t->level, gfs_domain_depth (domain));
    gchar name[] = "/tmp/level-x";
    name[11] = 48 + t->level;
    draw_level (domain, refine, t->level, name);
#endif
    t->level++;
  } while (t->refined);
#if DEBUG
  draw_all (domain, refine, "/tmp/all");
#endif
}

static void refine_terrain_destroy (GtsObject * object)
{
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (object);
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

  gchar * dir = NULL;
  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_STRING, "basename", TRUE},
      {GTS_STRING, "dir",      TRUE},
      {GTS_NONE}
    };
    gchar * basename = NULL;
    var[0].data = &basename;
    var[1].data = &dir;
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
    if (var[0].set) { g_free (t->basename); t->basename = basename; }
    if (!var[1].set) dir = g_strdup (default_dir);
  }
  else
    dir = g_strdup (default_dir);

  if (!strcmp (t->basename, "*")) { /* file globbing */
    gchar * pattern = g_strconcat (dir, "/*.Data", NULL);
    glob_t pglob;
    if (glob (pattern, GLOB_ERR, NULL, &pglob)) {
      gts_file_error (fp, "cannot find/open RSurface files in path:\n%s", pattern);
      g_free (pattern);
      g_free (dir);
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
	gts_file_error (fp, "cannot open RSurface `%s'", pglob.gl_pathv[i]);
	globfree (&pglob);
	g_free (dir);
	return;
      }
      if (t->basename) {
	gchar * dirbasename = g_strconcat (t->basename, ",", pglob.gl_pathv[i], NULL);
	g_free (t->basename);
	t->basename = dirbasename;
      }
      else
	t->basename = g_strdup (pglob.gl_pathv[i]);
      t->nrs++;
    }
    globfree (&pglob);
  }
  else { /* basename is of the form: set1,set2,set3... */
    gchar * basename = g_strdup (t->basename);
    if (dir) {
      g_free (t->basename);
      t->basename = NULL;
    }
    gchar * s = strtok (basename, ",");
    while (s) {
      t->rs = g_realloc (t->rs, (t->nrs + 1)*sizeof (RSurface *));
      if (dir) {
	gchar * fname = s[0] == '/' ? g_strdup (s) : g_strconcat (dir, "/", s, NULL);
	t->rs[t->nrs] = r_surface_open (fname, "r", -1);
	if (t->basename) {
	  gchar * dirbasename = g_strconcat (t->basename, ",", fname, NULL);
	  g_free (t->basename);
	  t->basename = dirbasename;
	  g_free (fname);
	}
	else
	  t->basename = fname;
      }
      else
	t->rs[t->nrs] = r_surface_open (s, "r", -1);
      if (!t->rs[t->nrs]) {
	gts_file_error (fp, "cannot open RSurface `%s'", s);
	g_free (basename);
	g_free (dir);
	return;
      }
      t->nrs++;
      s = strtok (NULL, ",");
    }
    g_free (basename);
  }
  g_free (dir);

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
  name = g_strjoin (NULL, t->name, "c", NULL);
  t->cond = gfs_domain_get_or_add_variable (domain, name, "Terrain condition number");
  g_free (name);

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

/* Initialize module */

const gchar * g_module_check_init (void);
void          gfs_module_read     (GtsFile * fp);
void          gfs_module_write    (FILE * fp);

const gchar * g_module_check_init (void)
{
  gchar * dir = getenv ("GFS_RSURFACES_PATH");
  if (dir)
    default_dir = dir;
  gfs_refine_terrain_class ();
  return NULL;
}

void gfs_module_read (GtsFile * fp)
{
  g_return_if_fail (fp != NULL);
}

void gfs_module_write (FILE * fp)
{
  g_return_if_fail (fp != NULL);
}
