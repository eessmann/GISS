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
#include "config.h"
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif /* HAVE_UNISTD_H */

#include "fluid.h"
#include "graphic.h"
#include "solid.h"
#include "poisson.h"
#include "init.h"

static gdouble K = 3.;
static gdouble L = 3.;
#ifndef   FTT_2D
static gdouble M = 3.;
#endif /* not FTT_2D */

static void fix_volume (FttCell * cell)
{
  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * s = GFS_STATE (cell)->solid;
    FttVector p;
    gdouble h = ftt_cell_size (cell);
    
    ftt_cell_pos (cell, &p);
    if (s->s[0] > 0. && s->s[0] < 1. &&
	s->s[1] > 0. && s->s[1] < 1.) {
      s->ca.x = p.x;
      s->a = (s->s[0] + s->s[1])/2.;
      if (s->s[3] == 1.) {
	g_assert (s->s[2] == 0.);
	s->ca.y = p.y + h*(s->a - 0.5);
      }
      else {
	g_assert (s->s[2] == 1. && s->s[3] == 0.);
	s->ca.y = p.y - h*(s->a - 0.5);	
      }
    }
    else if (s->s[2] > 0. && s->s[2] < 1. &&
	     s->s[3] > 0. && s->s[3] < 1.) {
      s->ca.y = p.y;
      s->a = (s->s[2] + s->s[3])/2.;
      if (s->s[1] == 1.) {
	g_assert (s->s[0] == 0.);
	s->ca.x = p.x + h*(s->a - 0.5);
      }
      else {
	g_assert (s->s[0] == 1. && s->s[1] == 0.);
	s->ca.x = p.x - h*(s->a - 0.5);	
      }
    }
    else if (s->s[0] > 0. && s->s[0] < 1. &&
	     s->s[2] > 0. && s->s[2] < 1.) {
      if (s->s[1] == 1.) {
	g_assert (s->s[3] == 1.);
	s->a = 1. - (1. - s->s[0])*(1. - s->s[2])/2.;
	s->ca.x = p.x + h*s->s[2]/2.;
	s->ca.y = p.y + h*s->s[0]/2.;
      }
      else {
	g_assert (s->s[3] == 0. && s->s[1] == 0.);
	s->a = s->s[0]*s->s[2]/2.;
      	s->ca.x = p.x + h*(1. - s->s[2])/2.;
	s->ca.y = p.y + h*(1. - s->s[0])/2.;
      }
    }
    else if (s->s[0] > 0. && s->s[0] < 1. &&
	     s->s[3] > 0. && s->s[3] < 1.) {
      if (s->s[1] == 1.) {
	g_assert (s->s[2] == 1.);
	s->a = 1. - (1. - s->s[0])*(1. - s->s[3])/2.;
	s->ca.x = p.x + h*s->s[3]/2.;
	s->ca.y = p.y - h*s->s[0]/2.;
      }
      else {
	g_assert (s->s[2] == 0. && s->s[1] == 0.);
	s->a = s->s[0]*s->s[3]/2.;
      	s->ca.x = p.x + h*(1. - s->s[3])/2.;
	s->ca.y = p.y - h*(1. - s->s[0])/2.;
      }
    }

    else if (s->s[1] > 0. && s->s[1] < 1. &&
	     s->s[2] > 0. && s->s[2] < 1.) {
      if (s->s[0] == 1.) {
	g_assert (s->s[3] == 1.);
	s->a = 1. - (1. - s->s[1])*(1. - s->s[2])/2.;
	s->ca.x = p.x - h*s->s[2]/2.;
	s->ca.y = p.y + h*s->s[1]/2.;
      }
      else {
	g_assert (s->s[3] == 0. && s->s[0] == 0.);
	s->a = s->s[1]*s->s[2]/2.;
      	s->ca.x = p.x - h*(1. - s->s[2])/2.;
	s->ca.y = p.y + h*(1. - s->s[1])/2.;
      }
    }
    else if (s->s[1] > 0. && s->s[1] < 1. &&
	     s->s[3] > 0. && s->s[3] < 1.) {
      if (s->s[0] == 1.) {
	g_assert (s->s[2] == 1.);
	s->a = 1. - (1. - s->s[1])*(1. - s->s[3])/2.;
	s->ca.x = p.x - h*s->s[3]/2.;
	s->ca.y = p.y - h*s->s[1]/2.;
      }
      else {
	g_assert (s->s[2] == 0. && s->s[0] == 0.);
	s->a = s->s[1]*s->s[3]/2.;
      	s->ca.x = p.x - h*(1. - s->s[3])/2.;
	s->ca.y = p.y - h*(1. - s->s[1])/2.;
      }
    }
    else
      g_assert_not_reached ();
  }
}

static void init_sine_div (FttCell * cell)
{
  FttVector cm;
  gdouble size;

  gfs_cell_cm (cell, &cm);
  size = ftt_cell_size (cell);
#if FTT_2D
  GFS_STATE (cell)->div = - M_PI*M_PI*(K*K + L*L)*sin (M_PI*K*cm.x)*sin(M_PI*L*cm.y);
#else  /* FTT_3D */
  GFS_STATE (cell)->div = - M_PI*M_PI*(K*K + L*L + M*M)*
    sin (M_PI*K*cm.x)*sin(M_PI*L*cm.y)*sin(M_PI*M*cm[2]);
#endif /* FTT_3D */
  if (GFS_IS_MIXED (cell))
    GFS_STATE (cell)->div *= size*size*GFS_STATE (cell)->solid->a; 
  else
    GFS_STATE (cell)->div *= size*size;
}

#define EPS 1e-3

static double dx (double x, double y, double (*f) (double, double))
{
  double f0 = (*f) (x, y);
  
  return (2.*((*f) (x + EPS, y) - f0)/3. -
             ((*f) (x + 2.*EPS, y) - f0)/12. -
          2.*((*f) (x - EPS, y) - f0)/3. +
             ((*f) (x - 2.*EPS, y) - f0)/12.)/EPS;
}

static double dx2 (double x, double y, double (*f) (double, double))
{
  double f0 = (*f) (x, y);
  
  return (4.*((*f) (x + EPS, y) - f0)/3. -
             ((*f) (x + 2.*EPS, y) - f0)/12. +
          4.*((*f) (x - EPS, y) - f0)/3. -
             ((*f) (x - 2.*EPS, y) - f0)/12.)/(EPS*EPS);
}

static double dy (double x, double y, double (*f) (double, double))
{
  double f0 = (*f) (x, y);
  
  return (2.*((*f) (x, y + EPS) - f0)/3. -
             ((*f) (x, y + 2.*EPS) - f0)/12. -
          2.*((*f) (x, y - EPS) - f0)/3. +
             ((*f) (x, y - 2.*EPS) - f0)/12.)/EPS;
}

static double dy2 (double x, double y, double (*f) (double, double))
{
  double f0 = (*f) (x, y);
  
  return (4.*((*f) (x, y + EPS) - f0)/3. -
             ((*f) (x, y + 2.*EPS) - f0)/12. +
          4.*((*f) (x, y - EPS) - f0)/3. -
             ((*f) (x, y - 2.*EPS) - f0)/12.)/(EPS*EPS);
}

static double annulus_p (double x, double y)
{
  return x*cos (M_PI*(0.25 - sqrt (x*x + y*y))/(0.25 - 0.5))/sqrt (x*x + y*y);
  //  return cos (M_PI*(0.25 - sqrt (x*x + y*y))/(0.25 - 0.5));
}

static void init_annulus_div (FttCell * cell)
{
  FttVector cm;
  gdouble size;

  gfs_cell_cm (cell, &cm);
  size = ftt_cell_size (cell);
  GFS_STATE (cell)->div = (dx2 (cm.x, cm.y, annulus_p) +
			   dy2 (cm.x, cm.y, annulus_p));
  if (GFS_IS_MIXED (cell))
    GFS_STATE (cell)->div *= size*size*GFS_STATE (cell)->solid->a; 
  else
    GFS_STATE (cell)->div *= size*size;
}

extern void gfs_face_gradient_flux_centered (const FttCellFace * face,
					     GfsGradient * g,
					     guint v,
					     gint max_level);

static void cell_annulus_error (FttCell * cell, gdouble * bias)
{
#if 1
  FttVector pos;

  ftt_cell_pos (cell, &pos);
  GFS_STATE (cell)->v = GFS_STATE (cell)->dp = 
    annulus_p (pos.x, pos.y) - GFS_STATE (cell)->p -
    *bias;
#else
  FttCellFace f;
  FttVector p;
  GtsVector ca;

  f.cell = cell;
  f.d = FTT_RIGHT;
  f.neighbor = ftt_cell_neighbor (f.cell, f.d);
  ftt_face_pos (&f, &p);
  if (f.neighbor) {
    gfs_face_ca (&f, ca);
    GFS_STATE (cell)->dp = dx (ca.x, ca.y, annulus_p) - 
      GFS_STATE (cell)->f[f.d].un;
  }
  else
    GFS_STATE (cell)->dp = 0.;
#endif
}

static void init_day_pressure (FttCell * cell)
{
  FttVector p;
  gdouble r2, theta;

  ftt_cell_pos (cell, &p);
  r2 = p.x*p.x + p.y*p.y;
  theta = atan2 (p.y, p.x);
  GFS_STATE (cell)->p = r2*r2*cos (3.*theta);
}

static void init_day_div (FttCell * cell)
{
  FttVector cm;
  gdouble r2, theta;
  gdouble h = ftt_cell_size (cell);
  
  gfs_cell_cm (cell, &cm);
  r2 = cm.x*cm.x + cm.y*cm.y;
  theta = atan2 (cm.y, cm.x);
  GFS_STATE (cell)->div = 7.*r2*cos (3.*theta);
  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * s = GFS_STATE (cell)->solid;
    GtsVector n;
    gdouble a;
    gdouble x = s->ca.x;
    gdouble y = s->ca.y;

    n[0] = s->s[1] - s->s[0];
    n[1] = s->s[3] - s->s[2];
    r2 = x*x + y*y;
#if 0
    a = sqrt (n[0]*n[0] + n[1]*n[1]);
    //    fprintf (stderr, "%g %g\n", a, s->sa);
    n[0] /= a;
    n[1] /= a;
    GFS_STATE (cell)->div = GFS_STATE (cell)->div*h*h*s->a -
      ((4.*x*x*x*x - 3.*x*x*y*y - 3.*y*y*y*y)*n[0] + 
       (-x*y*(5.*x*x + 9.*y*y))*n[1])/sqrt (r2)*h*s->sa;
#else
    GFS_STATE (cell)->div = GFS_STATE (cell)->div*h*h*s->a -
      ((4.*x*x*x*x - 3.*x*x*y*y - 3.*y*y*y*y)*n[0] + 
       (-x*y*(5.*x*x + 9.*y*y))*n[1])/sqrt (r2)*h;
#endif
  }
  else {
    FttCellFace f;
    FttCellNeighbors n;

    GFS_STATE (cell)->div *= h*h;
    ftt_cell_neighbors (cell, &n);
    f.cell = cell;
    for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++)
      if (n.c[f.d] == NULL) {
	FttVector p;
	gdouble x, y, r2;
	gdouble nx[4] = { 1., -1., 0., 0.};
	gdouble ny[4] = { 0., 0., 1., -1.};

	ftt_face_pos (&f, &p);
	x = p.x; y = p.y;
	r2 = x*x + y*y;
	GFS_STATE (cell)->div -=
	  ((4.*x*x*x*x - 3.*x*x*y*y - 3.*y*y*y*y)*nx[f.d] + 
	   (-x*y*(5.*x*x + 9.*y*y))*ny[f.d])/sqrt (r2)*h;
      }
  }
}

static void cell_day_error (FttCell * cell, gdouble * bias)
{
  FttVector p;
  gdouble r2, theta;

  ftt_cell_pos (cell, &p);
  r2 = p.x*p.x + p.y*p.y;
  theta = atan2 (p.y, p.x);
  GFS_STATE (cell)->v = GFS_STATE (cell)->dp = 
    r2*r2*cos (3.*theta) - GFS_STATE (cell)->p - *bias;
}

static void cell_day_truncation (FttCell * cell, GfsNorm * norm)
{
  gdouble h = ftt_cell_size (cell);

  GFS_STATE (cell)->v = GFS_STATE (cell)->dp = GFS_STATE (cell)->res/(h*h);
  gfs_norm_add (norm, GFS_STATE (cell)->dp, h*h);
}

static void sum_volume (FttCell * cell, GtsRange * vol)
{
  gdouble size = ftt_cell_size (cell);

  if (GFS_IS_MIXED (cell))
    gts_range_add_value (vol, size*size*GFS_STATE (cell)->solid->a);
  else
    gts_range_add_value (vol, size*size);
}

static void add_ddiv (FttCell * cell, gdouble * ddiv)
{
  gdouble size = ftt_cell_size (cell);

  if (GFS_IS_MIXED (cell))
    GFS_STATE (cell)->div += size*size*GFS_STATE (cell)->solid->a*(*ddiv);
  else
    GFS_STATE (cell)->div += size*size*(*ddiv);
}

static void correct_div (FttCell * root, gboolean * verbose)
{
  GtsRange div, vol;
  gdouble ddiv;

  div = gfs_stats_variable (root, gfs_div, FTT_TRAVERSE_LEAFS, -1);
  gts_range_init (&vol);
  ftt_cell_traverse (root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) sum_volume, &vol);
  gts_range_update (&vol);
  ddiv = - div.mean/vol.mean;

  if (*verbose)
    fprintf (stderr, "Correcting divergence: %g\n", ddiv);
  ftt_cell_traverse (root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) add_ddiv, &ddiv);
}

static void cell_error (FttCell * cell, gdouble * avg)
{
  FttVector pos;

  ftt_cell_pos (cell, &pos);
#if FTT_2D
  GFS_STATE (cell)->dp = (*avg) + 
    sin (M_PI*K*pos.x)*sin (M_PI*L*pos.y) - GFS_STATE (cell)->p;
#else  /* FTT_3D */
  GFS_STATE (cell)->dp = (*avg) + 
    sin (M_PI*K*pos.x)*sin (M_PI*L*pos.y)*sin (M_PI*M*pos.z) - GFS_STATE (cell)->p;
#endif /* FTT_3D */
}

static void cell_error_ref (FttCell * cell, gpointer * data)
{
  FttCell * root = data[0], * locate;
  gdouble * bias = data[1];
  gdouble * dpmax = data[2];
  gdouble * amax = data[3], a = GFS_IS_MIXED (cell) ? 
    GFS_STATE (cell)->solid->a : 1.;
  guint level = ftt_cell_level (cell);
  FttVector pos;

  ftt_cell_pos (cell, &pos);
  locate = ftt_cell_locate (root, pos, level);
  g_assert (locate && ftt_cell_level (locate) == level);
  GFS_STATE (cell)->dp = GFS_STATE (locate)->p - *bias - GFS_STATE (cell)->p;
  if (fabs (GFS_STATE (cell)->dp) > *dpmax) {
    *dpmax = fabs (GFS_STATE (cell)->dp);
    *amax = a;
  }
}

static void cell_error_hist (FttCell * cell, gpointer * data)
{
  FttCell * root = data[0], * locate;
  gdouble * bias = data[1];
  guint level = ftt_cell_level (cell);
  FttVector pos;

  ftt_cell_pos (cell, &pos);
  locate = ftt_cell_locate (root, pos, level);
  g_assert (locate && ftt_cell_level (locate) == level);
  GFS_STATE (cell)->dp = GFS_STATE (locate)->p - *bias - GFS_STATE (cell)->p;
  printf ("%g %g\n", GFS_IS_MIXED (cell) ? 
	  GFS_STATE (cell)->solid->a : 1., GFS_STATE (cell)->dp);
}

static void weighted_residual (FttCell * cell)
{
  gdouble size = ftt_cell_size (cell);

  GFS_STATE (cell)->res /= size*size;
}

static gboolean refine (FttCell * cell, guint * minlevel)
{
  if (ftt_cell_level (cell) < *minlevel)
    return TRUE;
  return FALSE;
}

static gboolean refine_solid (FttCell * cell, guint * maxlevel)
{
  if (GFS_IS_MIXED (cell) && ftt_cell_level (cell) < *maxlevel)
    return TRUE;
  return FALSE;
}

static gboolean refine_inside_solid (FttCell * cell, gpointer * data)
{
  GNode * tree = data[0];
  guint * maxlevel = data[1];
  gboolean * is_open = data[2];

  if (ftt_cell_level (cell) < *maxlevel) {
    gdouble size = ftt_cell_size (cell)/2.;
    GtsPoint * p = gts_point_new (gts_point_class (), 0., 0., 0.);
    FttVector pos;
    static gint shift[FTT_CELLS][3] = 
#if FTT_2D
    {{-1,-1,0},{-1,1,0},{1,1,0},{1,-1,0}};
#else  /* FTT_3D */
    {{-1,-1,-1},{-1,1,-1},{1,1,-1},{1,-1,-1},
     {-1,-1, 1},{-1,1, 1},{1,1, 1},{1,-1, 1}};
#endif /* FTT_3D */
    guint i;
    
    ftt_cell_pos (cell, &pos);
    for (i = 0; i < FTT_CELLS; i++) {
      p->x = pos.x + size*shift[i][0];
      p->y = pos.y + size*shift[i][1];
      p->z = pos.z + size*shift[i][2];
      if (!gts_point_is_inside_surface (p, tree, *is_open)) {
	gts_object_destroy (GTS_OBJECT (p));
	return TRUE;
      }
    }
    gts_object_destroy (GTS_OBJECT (p));
  }
  return FALSE;
}

static gboolean refine_randomly (FttCell * cell, gpointer * data)
{
  guint * levelmax = data[0];
  gdouble * probability = data[1];

  if (ftt_cell_level (cell) < *levelmax && 
      (gdouble) rand () < *probability*(gdouble) RAND_MAX)
    return TRUE;
  return FALSE;
}

static void init_solid_fractions (FttCell * cell, gpointer * data)
{
  gboolean * is_open = data[2];
  GfsDomain * domain = data[3];

  gfs_cell_init (cell, domain);
  gfs_cell_init_solid_fractions (cell, 
				 data[0], data[1], *is_open, 
				 TRUE, (FttCellCleanupFunc) gfs_cell_cleanup,
				 NULL);
}

#define GFS_FLAG_MARKED  (1 << 5)
#define GFS_IS_MARKED(cell) (((cell)->flags & GFS_FLAG_MARKED) != 0)

typedef enum
{
  POISSON,
  RELAX,
  CUSTOM
} Method;

static GfsBox * solution (guint levels,
			  guint levelmax,
			  gboolean relative,
			  GtsSurface * solid,
			  gboolean nosolid,
			  gdouble probability,
			  gboolean verbose,
			  gboolean annulus,
                          gboolean day,
			  gboolean truncation,
			  gboolean fix_vol,
			  Method method,
			  guint vcycle,
			  guint nrelax,
			  gdouble * ttime,
			  GfsDomain ** retdomain)
{
  GfsBox * box;
  GfsDomain * domain;
  GTimer * timer;
  guint cycle;
  gdouble total_time = 0.;

  *retdomain = domain = 
   GFS_DOMAIN (gts_graph_new (GTS_GRAPH_CLASS (gfs_domain_class ()), 
			      GTS_GNODE_CLASS (gfs_box_class ()),
			      gts_gedge_class ()));
  box = gfs_box_new (gfs_box_class ());
  box->root = ftt_cell_new ((FttCellInitFunc) gfs_cell_init, domain);
  ftt_cell_refine (box->root, (FttCellRefineFunc) refine, &levels, 
		   (FttCellInitFunc) gfs_cell_init, domain);
  gts_container_add (GTS_CONTAINER (domain), GTS_CONTAINEE (box));

  if (relative)
    levelmax += levels;

  if (solid) {
    GNode * tree = gts_bb_tree_surface (solid);
    gboolean is_open = FALSE;

    if (gts_surface_volume (solid) < 0.)
      is_open = TRUE;

    if (nosolid) {
      gpointer data[3];
	
      data[0] = tree;
      data[1] = &levelmax;
      data[2] = &is_open;
      ftt_cell_refine (box->root, 
		       (FttCellRefineFunc) refine_inside_solid, data,
		       (FttCellInitFunc) gfs_cell_init, domain);
      ftt_cell_refine (box->root, 
		       (FttCellRefineFunc) ftt_refine_corner, NULL,
		       (FttCellInitFunc) gfs_cell_init, domain);
    }
    else {
      gpointer data[4];
	
      gfs_cell_init_solid_fractions (box->root, solid, tree, is_open, 
				     TRUE, 
				     (FttCellCleanupFunc) gfs_cell_cleanup,
				     NULL);

      data[0] = solid;
      data[1] = tree;
      data[2] = &is_open;
      data[3] = domain;
      if (levelmax > 0)
	ftt_cell_refine (box->root, 
			 (FttCellRefineFunc) refine_solid, &levelmax,
			 (FttCellInitFunc) init_solid_fractions, data);	
      ftt_cell_refine (box->root, 
		       (FttCellRefineFunc) ftt_refine_corner, NULL,
		       (FttCellInitFunc) init_solid_fractions, data);

      g_assert (gfs_cell_check_solid_fractions (box->root, 
						solid, tree, is_open));
      if (fix_vol)
	gfs_domain_cell_traverse (domain, 
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				  (FttCellTraverseFunc) fix_volume, NULL);
    }
    gts_bb_tree_destroy (tree, TRUE);
  }
  else {
    if (probability > 0.) {
      gpointer data[2];

      data[0] = &levelmax;
      data[1] = &probability;
      ftt_cell_refine (box->root, 
		       (FttCellRefineFunc) refine_randomly, data,
		       (FttCellInitFunc) gfs_cell_init, domain);
      ftt_cell_refine (box->root, 
		       (FttCellRefineFunc) ftt_refine_corner, NULL,
		       (FttCellInitFunc) gfs_cell_init, domain);
    }
    gfs_cell_fluid (box->root);
  }
  
  if (annulus)
    ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) init_annulus_div, NULL);
  else if (day) {
    ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) init_day_div, NULL);
    if (truncation) {
      ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			 (FttCellTraverseFunc) init_day_pressure, NULL);
      gfs_poisson_coefficients (domain, NULL, 1.);
      return box;
    }
  }
  else
    ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) init_sine_div, NULL);
  correct_div (box->root, &verbose); /* enforce solvability condition */
  if (verbose) {
    GtsRange stats = gfs_stats_variable (box->root, gfs_div, 
					 FTT_TRAVERSE_LEAFS, -1);
    fprintf (stderr, 
	     "initial div min: %g avg: %g| %g max: %g n: %7d\n",
	     stats.min, stats.mean, stats.stddev, stats.max, stats.n);
  }
  
  levels = gfs_domain_depth (domain);
  gfs_poisson_coefficients (domain, NULL, 1.);
  gfs_residual (domain, FTT_DIMENSION, FTT_TRAVERSE_LEAFS, -1, gfs_p, gfs_div, gfs_res);

  timer = g_timer_new ();
  g_timer_start (timer);

  if (method == CUSTOM) {
  }
  else {
    if (method == RELAX)
      ftt_cell_traverse (box->root, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			 (FttCellTraverseFunc) gfs_get_from_below_extensive, 
			 gfs_div);
    for (cycle = 0; cycle < vcycle; cycle++) {
      if (verbose) {
	GfsNorm norm;
	
	g_timer_stop (timer);
	total_time += g_timer_elapsed (timer, NULL);
	gfs_residual (domain, FTT_DIMENSION, FTT_TRAVERSE_LEAFS, -1, gfs_p, gfs_div, gfs_res);
	norm = gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, 1.);
	fprintf (stderr, 
		 "div n: %d time: %g first: %g second: %g infty: %g w: %g\n",
		 cycle, total_time,
		 norm.first, norm.second, norm.infty, norm.w);
	
	g_timer_start (timer);
      }
      switch (method) {
      case POISSON:
	gfs_poisson_cycle (domain, FTT_DIMENSION, 0, levels, nrelax, gfs_p, gfs_div);
	break;
      case RELAX: {
	gfs_relax (domain, FTT_DIMENSION, -1, gfs_p, gfs_div);
	break;
      }
      default :
	break;
      }
    }
  }

  g_timer_stop (timer);
  total_time += g_timer_elapsed (timer, NULL);
  *ttime = total_time;

  return box;
}

int main (int argc, char * argv[])
{
  GfsBox * box = NULL, * reference = NULL;
  GfsDomain * domain = NULL;
  guint levels, levelmax = 0, vcycle, nrelax, levels_ref = 0, l;
  GtsRange stats;
  GfsNorm norm;
  int c = 0;

  gboolean verbose = FALSE;
  gboolean output_mesh = FALSE;
  gboolean output_tree = FALSE;
  gboolean nosolid = FALSE;
  gdouble probability = -1.;
  gboolean relative = FALSE;
  gboolean histogram = FALSE;
  gboolean annulus = FALSE;
  gboolean day = FALSE;
  gboolean truncation = FALSE;
  gboolean fix_vol = FALSE;
  Method method = POISSON;
  GfsVariable * output = NULL;
  gint olevel = -1;
  gint minlevel = -1;

  GtsSurface * solid = NULL;

  gdouble total_time = 0.;

  GtsBBox bbox;

  bbox.x1 = bbox.y1 = -0.5;
  bbox.x2 = bbox.y2 = 0.5;
  bbox.z1 = bbox.z2 = 0.;

  gfs_init (&argc, &argv);

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"fixvolume", no_argument, NULL, 'f'},
      {"truncation", no_argument, NULL, 'T'},
      {"day", no_argument, NULL, 'd'},
      {"annulus", no_argument, NULL, 'A'},
      {"histogram", no_argument, NULL, 'i'},
      {"solid", required_argument, NULL, 's'},
      {"refine", required_argument, NULL, 'r'},
      {"relative", no_argument, NULL, 'R'},
      {"nosolid", no_argument, NULL, 'n'},
      {"random", required_argument, NULL, 'a'},
      {"output", required_argument, NULL, 'o'},
      {"level", required_argument, NULL, 'l'},
      {"mesh", no_argument, NULL, 'M'},
      {"method", required_argument, NULL, 'm'},
      {"tree", no_argument, NULL, 't'},
      {"error", required_argument, NULL, 'e'},
      {"order", required_argument, NULL, 'x'},
      {"kx", required_argument, NULL, 'k'},
      {"ky", required_argument, NULL, 'u'},
      {"kz", required_argument, NULL, 'w'},
      {"help", no_argument, NULL, 'h'},
      {"verbose", no_argument, NULL, 'v'},
      {"version", no_argument, NULL, 'V'}
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "hvo:m:Ml:s:Vr:nRa:te:x:k:u:w:iAdTf",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "hvo:m:Ml:s:Vr:nRa:te:x:k:u:w:iAdTf"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'f': /* fixvolume */
      fix_vol = TRUE;
      break;
    case 'T': /* truncation */
      truncation = TRUE;
      break;
    case 'd': /* day */
      day = TRUE;
      break;
    case 'A': /* annulus */
      annulus = TRUE;
      break;
    case 'i': /* histogram */
      histogram = TRUE;
      break;
    case 's': { /* solid */
      FILE * fptr = optarg[0] == '-' ? stdin : fopen (optarg, "rt");
      GtsFile * fp;
      
      if (fptr == NULL) {
	fprintf (stderr, "poisson: cannot open file `%s'\n", optarg);
	return 1;
      }
      
      solid = gts_surface_new (gts_surface_class (),
			       gts_face_class (),
			       gts_edge_class (),
			       gts_vertex_class ());
      fp = gts_file_new (fptr);
      if (gts_surface_read (solid, fp)) {
	fprintf (stderr, 
		 "poisson: file `%s' is not a valid GTS file\n"
		 "%s:%d:%d: %s\n",
		 fp->fp == stdin ? "stdin" : optarg, 
		 fp->fp == stdin ? "stdin" : optarg, 
		 fp->line, fp->pos, fp->error);
	return 1;
      }
      gts_file_destroy (fp);
      fclose (fptr);
      break;
    }
    case 'r': /* refine */
      levelmax = atoi (optarg);
      break;
    case 'R': /* relative */
      relative = TRUE;
      break;
    case 'n': /* nosolid */
      nosolid = TRUE;
      break;
    case 'a': /* random */
      probability = atof (optarg);
      break;
    case 'm': /* method */
      if (!strcmp (optarg, "POISSON"))
	method = POISSON;
      else if (!strcmp (optarg, "RELAX"))
	method = RELAX;
      else if (!strcmp (optarg, "CUSTOM"))
	method = CUSTOM;
      else {
	fprintf (stderr, "poisson: bad method option `%s'\n", optarg);
	fprintf (stderr, "Try `poisson --help' for more information.\n");
	return 1; /* failure */
      }
      break;
    case 'o': /* output */
      if (!strcmp (optarg, "GFS_P"))
	output = gfs_p;
      else if (!strcmp (optarg, "GFS_DIV"))
	output = gfs_div;
      else if (!strcmp (optarg, "GFS_RES"))
	output = gfs_res;
      else if (!strcmp (optarg, "GFS_DP"))
	output = gfs_dp;
      else {
	fprintf (stderr, "poisson: bad output option `%s'\n", optarg);
	fprintf (stderr, "Try `poisson --help' for more information.\n");
	return 1; /* failure */
      }
      break;
    case 'l': /* level */
      olevel = atoi (optarg);
      break;
    case 'M': /* mesh */
      output_mesh = TRUE;
      break;
    case 't': /* tree */
      output_tree = TRUE;
      break;
    case 'e': /* error */
      levels_ref = atoi (optarg);
      break;
    case 'x': /* order */
      minlevel = atoi (optarg);
      break;
    case 'k': /* kx */
      K = atof (optarg);
      break;
    case 'u': /* ky */
      L = atof (optarg);
      break;
#if (!FTT_2D)
    case 'w': /* kz */
      M = atof (optarg);
      break;
#endif /* not FTT_2D */
    case 'v': /* verbose */
      verbose = TRUE;
      break;
    case 'h': /* help */
      fprintf (stderr,
     "Usage: poisson [OPTION] LEVELS CYCLE NRELAX\n"
     "Solve a simple poisson equation on various grids\n"
     "\n"
     "  -f    --fixvolume   fix volume of cut cells using a piecewise linear approximation\n"
     "  -T    --truncation  compute truncation error (only for Day problem)\n"
     "  -d    --day         solves Day problem\n"
     "  -A    --annulus     solves annulus problem\n"
     "  -s F  --solid=F     read GTS file F describing solid boundaries\n"
     "                      (if F is `-' the standard input is used)\n"
     "  -r L  --refine=L    refines up to level L\n"
     "  -R    --relative    makes the refinement relative i.e. the maximum level is\n"
     "                      LEVELS + L instead of L\n"
     "  -n    --nosolid     removes the solid while keeping the generated mesh\n"
     "  -a    --random      refine the mesh randomly\n"
     "  -o V  --output=V    output a GTS file representing variable V, where V is one of:\n"
     "                      GFS_P, GFS_DIV, GFS_RES, GFS_DP\n"
     "  -l L  --level=L     output values for level L (default is leaves)\n"
     "  -M    --mesh        output a OOGL file representing the mesh\n"
     "  -t    --tree        output a GFS representation of the whole cell tree\n"
     "  -e L  --error=L     use an evaluation of the solution with a maximum level L\n"
     "                      as the reference solution for the error computation\n"
     "  -i    --histogram   outputs list of cell fractions and associated errors\n"
     "  -x L  --order=L     run several simulations of varying size (starting at L)\n"
     "                      and outputs the corresponding errors\n"
     "  -k N  --kx=N        set wave number in direction x to N (default is 3)\n"
     "  -u N  --ky=N        set wave number in direction y to N (default is 3)\n"
#if (!FTT_2D)
     "  -w N  --kz=N        set wave number in direction z to N (default is 3)\n"
#endif /* not FTT_2D */
     "  -m M  --method=M    use method M, where M is one of:\n"
     "                      POISSON, RELAX, CUSTOM\n"
     "                      (default is POISSON)\n"
     "  -v    --verbose     display info, timing etc...\n"
     "  -h    --help        display this help and exit\n"
     "  -V    --version     output version information and exit\n"
     "\n"
     "Reports bugs to %s\n",
	       FTT_MAINTAINER);
      return 0; /* success */
      break;
    case 'V': /* version */
      fprintf (stderr, 
	       "poisson: using %dD libgfs version %s\n"
	       "compiled with flags: %s\n",
	       FTT_DIMENSION,
	       GFS_VERSION,
	       GFS_COMPILATION_FLAGS);
      return 0; /* succes */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `poisson --help' for more information.\n");
      return 1; /* failure */
    }
  }

  if (optind >= argc) { /* missing LEVELS */  
    fprintf (stderr, 
	     "poisson: missing LEVELS\n"
	     "Try `poisson --help' for more information.\n");
    return 1; /* failure */
  }
  levels = atoi (argv[optind++]);

  if (optind >= argc) { /* missing CYCLE */  
    fprintf (stderr, 
	     "poisson: missing CYCLE\n"
	     "Try `poisson --help' for more information.\n");
    return 1; /* failure */
  }
  vcycle = atoi (argv[optind++]);

  if (optind >= argc) { /* missing NRELAX */  
    fprintf (stderr, 
	     "poisson: missing NRELAX\n"
	     "Try `poisson --help' for more information.\n");
    return 1; /* failure */
  }
  nrelax = atoi (argv[optind++]);
  
  if (levels_ref > levels) {
    reference = solution (levels_ref, levelmax, relative,
			  solid, nosolid, probability, 
			  FALSE, annulus, day, truncation, fix_vol,
			  method, vcycle, nrelax, 
			  &total_time, &domain);
    ftt_cell_traverse (reference->root, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
		       (FttCellTraverseFunc) gfs_get_from_below_intensive,
		       gfs_p);
  }

  if (minlevel < 0)
    minlevel = levels;
  else
    verbose = FALSE;

  for (l = minlevel; l <= levels; l++) {
    if (box)
      ftt_cell_destroy (box->root, (FttCellCleanupFunc) gfs_cell_cleanup, 
			NULL);
    box = solution (l, levelmax, relative, solid, nosolid, probability, 
		    verbose, annulus, day, truncation, fix_vol,
		    method, vcycle, nrelax, &total_time, &domain);

    gfs_residual (domain, FTT_DIMENSION, FTT_TRAVERSE_LEAFS, -1, gfs_p, gfs_div, gfs_res);
    if (verbose) {
      norm = gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, 1.);
      fprintf (stderr, 
	       "div n: %d time: %g first: %g second: %g infty: %g w: %g\n",
	       vcycle, total_time,
	       norm.first, norm.second, norm.infty, norm.w);
      
      stats = gfs_stats_variable (box->root, gfs_p, FTT_TRAVERSE_LEAFS, -1);
      fprintf (stderr, 
	       "p   min: %g avg: %g| %g max: %g n: %7d\n",
	       stats.min, stats.mean, stats.stddev, stats.max, stats.n);
      
      fprintf (stderr, "Time: %g s Speed: %.0f sites.iterations/s\n", 
	       total_time,
	       stats.n*vcycle/total_time);
    }
    
    if (verbose || minlevel < levels) {
      if (reference) {
	gpointer data[4];
	gdouble bias = 0.;
	gdouble dpmax = 0., amax;
	
	data[0] = reference->root;
	data[1] = &bias;
	data[2] = &dpmax;
	data[3] = &amax;
	ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			   (FttCellTraverseFunc) cell_error_ref, data);
	norm = gfs_norm_variable (box->root, gfs_dp, FTT_TRAVERSE_LEAFS, -1);
	data[1] = &norm.bias;
	ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			   (FttCellTraverseFunc) cell_error_ref, data);
	fprintf (stderr, "level: %d dpmax: %g amax: %g\n", 
		 l, dpmax, amax);
	norm = gfs_norm_variable (box->root, gfs_dp, FTT_TRAVERSE_LEAFS, -1);
      }
      else if (annulus) {
	norm.bias = 0.;
	ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			   (FttCellTraverseFunc) cell_annulus_error, 
			   &norm.bias);
	norm = gfs_norm_variable (box->root, gfs_dp, FTT_TRAVERSE_LEAFS, -1);
	ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			   (FttCellTraverseFunc) cell_annulus_error, 
			   &norm.bias);
	norm = gfs_norm_variable (box->root, gfs_dp, FTT_TRAVERSE_LEAFS, -1);
      }
      else if (day) {
	if (truncation) {
	  gfs_norm_init (&norm);
	  ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			     (FttCellTraverseFunc) cell_day_truncation, 
			     &norm);
	  gfs_norm_update (&norm);
	}
	else {
	  norm.bias = 0.;
	  ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			     (FttCellTraverseFunc) cell_day_error, 
			     &norm.bias);
	  norm = gfs_norm_variable (box->root, gfs_dp, FTT_TRAVERSE_LEAFS, -1);
	  ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			     (FttCellTraverseFunc) cell_day_error, 
			     &norm.bias);
	  norm = gfs_norm_variable (box->root, gfs_dp, FTT_TRAVERSE_LEAFS, -1);
	}
      }
      else {
	stats = gfs_stats_variable (box->root, gfs_p, FTT_TRAVERSE_LEAFS, -1);
	ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			   (FttCellTraverseFunc) cell_error, &stats.mean);
	norm = gfs_norm_variable (box->root, gfs_dp, FTT_TRAVERSE_LEAFS, -1);
      }
      fprintf (stderr, "total err level: %d first: %g second: %g infty: %g bias: %g w: %g\n",
	       l, norm.first, norm.second, norm.infty, norm.bias, norm.w);
      norm = gfs_norm_variable (box->root, gfs_dp, FTT_TRAVERSE_LEVEL, 
				ftt_cell_depth (box->root));
      fprintf (stderr, "refined err level: %d first: %g second: %g infty: %g bias: %g w: %g\n",
	       l, norm.first, norm.second, norm.infty, norm.bias, norm.w);
    }

    if (histogram && reference) {
      gpointer data[4];
      gdouble bias = 0.;
      gdouble dpmax = 0., amax;
      
      data[0] = reference->root;
      data[1] = &bias;
      data[2] = &dpmax;
      data[3] = &amax;
      ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			 (FttCellTraverseFunc) cell_error_ref, data);
      norm = gfs_norm_variable (box->root, gfs_dp, FTT_TRAVERSE_LEAFS, -1);
      data[1] = &norm.bias;
      ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			 (FttCellTraverseFunc) cell_error_hist, data);
    }
  }

  if (output != NULL) {
    ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		       (FttCellTraverseFunc) weighted_residual, NULL);
    if (olevel >= 0) {
      ftt_cell_traverse (box->root, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			 (FttCellTraverseFunc) gfs_get_from_below_intensive,
			 gfs_p);
      gfs_write_gts (domain, output, FTT_TRAVERSE_LEVEL, olevel, 
		     &bbox, stdout);
      //gfs_write_gnuplot (box->root, output, FTT_TRAVERSE_LEVEL, olevel, 
      //                  bbox, stdout);
    }
    else {
      GtsRange stats = gfs_domain_stats_variable (domain, output, 
						  FTT_TRAVERSE_LEAFS, -1);

      gfs_write_gts (domain, output, FTT_TRAVERSE_LEAFS, -1, 
		     &bbox, stdout);
      //gfs_write_gnuplot (box->root, output, FTT_TRAVERSE_LEAFS, -1, 
      //            &bbox, stdout);
      //      gfs_write_squares (domain, output, stats.min, stats.max, 
      //			 FTT_TRAVERSE_LEAFS, -1, NULL, stdout);
    }
  }
  else if (output_mesh) {
    if (olevel >= 0) {
      printf ("LIST {\n");
      ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, olevel,
			 (FttCellTraverseFunc) ftt_cell_draw, stdout);
      printf ("}\n");
    }
    else {
      gfs_draw_levels (box->root, stdout);
      if (solid) {
	printf ("(geometry \"solid\" { =\n");
	gts_surface_write_oogl (solid, stdout);
	printf ("})\n");
      }
    }
  }
  else if (output_tree) {
    ftt_cell_traverse (box->root, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
		       (FttCellTraverseFunc) gfs_get_from_below_intensive,
		       gfs_p);
    if (olevel >= 0)
      domain->max_depth_write = olevel;
    gts_graph_write (GTS_GRAPH (domain), stdout);
  }

  return 0;
}
