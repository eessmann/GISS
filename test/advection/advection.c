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

#include <stdlib.h>
#include <math.h>
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
#include "advection.h"
#include "init.h"

#ifndef PI
# define PI 3.14159265359
#endif

static gboolean refine (FttCell * cell, guint * maxlevel)
{
  if (ftt_cell_level (cell) < *maxlevel)
    return TRUE;
  return FALSE;
}

static gboolean refine_left (FttCell * cell, guint * maxlevel)
{
  if (ftt_cell_level (cell) < *maxlevel) {
    FttVector pos;

    ftt_cell_pos (cell, &pos);
    if (pos.x < 0.)
      return TRUE;
  }
  return FALSE;
}

static gboolean refine_right (FttCell * cell, guint * maxlevel)
{
  if (ftt_cell_level (cell) < *maxlevel) {
    FttVector pos;

    ftt_cell_pos (cell, &pos);
    if (pos.x > 0.)
      return TRUE;
  }
  return FALSE;
}

static gboolean refine_box (FttCell * cell, guint * maxlevel)
{
  if (ftt_cell_level (cell) < *maxlevel) {
    FttVector pos;

    ftt_cell_pos (cell, &pos);
    if (pos.x < 0.25 && pos.x > -0.25 &&
	pos.y < 0.25 && pos.y > -0.25)
      return TRUE;
  }
  return FALSE;
}

static void copy_fraction (FttCell * cell)
{
  if (GFS_IS_FLUID (cell))
    GFS_STATE (cell)->p = 1.;
  else
    GFS_STATE (cell)->p = GFS_STATE (cell)->solid->a;
}

static void smooth_gaussian (FttCell * cell)
{
  FttVector pos;
  gdouble r2, coeff;

  ftt_cell_pos (cell, &pos);
  r2 = pos.x*pos.x + pos.y*pos.y;
  coeff = 20. + 20000.*r2*r2*r2*r2;
  GFS_STATE (cell)->p = 
    (1. + cos(20.*pos.x)*cos(20.*pos.y))*exp (-coeff*r2)/2.;
}

static void smooth_periodic (FttCell * cell)
{
  FttVector pos;

  ftt_cell_pos (cell, &pos);
  GFS_STATE (cell)->p = (1. + cos(2.*PI*pos.x)*sin(2.*PI*pos.y))/2.;
}

static void rotate (FttCell * cell, gdouble * speed)
{
  FttVector pos;

  ftt_cell_pos (cell, &pos);

  GFS_STATE (cell)->u = - 2.*(*speed)*pos.y;
  GFS_STATE (cell)->v = 2.*(*speed)*pos.x;
}

static void translate (FttCell * cell, gdouble * speed)
{
  gdouble theta = atan (1./2.);

  GFS_STATE (cell)->u = cos (theta)*(*speed);
  GFS_STATE (cell)->v = sin (theta)*(*speed);
}

static void update_fraction (FttCell * cell)
{
  GFS_STATE (cell)->p += GFS_STATE (cell)->div;
}

static void compute_error (FttCell * cell, 
			   FttCell * ref)
{
  FttVector pos;
  FttCell * locate;
  
  ftt_cell_pos (cell, &pos);
  locate = ftt_cell_locate (ref, pos, -1);
  g_assert (locate && ftt_cell_level (locate) == ftt_cell_level (cell));
  GFS_STATE (cell)->res = GFS_STATE (cell)->p - GFS_STATE (locate)->p;
}

static void total_mass (FttCell * cell,
			gdouble * mass)
{
  gdouble size = ftt_cell_size (cell);

  *mass += GFS_STATE (cell)->p*size*size;
}

int main (int argc, char * argv[])
{
  FttCell * ref = NULL;
  GfsBox * box;
  GfsDomain * domain;
  gdouble cfl = 0.5;
  guint maxlevel = 6, t;
  guint rightrefine = 0;
  guint leftrefine = 0;
  guint boxrefine = 0;
  GtsRange stats;
  gdouble total_mass_before = 0., total_mass_after = 0.;
  guint tmax, twrite;
  GfsAdvectionParams advection_params = { 
    0.8, 1., NULL, NULL,
    gfs_center_van_leer_gradient,
    FALSE, TRUE,
    gfs_face_advection_flux
  };

  int c = 0;
  gboolean verbose = FALSE;
  gboolean display_error = FALSE;
  gboolean periodic = FALSE;
  gboolean gtsfiles = TRUE;
  void (* transform) (FttCell *, gdouble *) = rotate;
  void (* smooth) (FttCell *) = NULL;

  GtsBBox bbox;

  bbox.x1 = bbox.y1 = -0.5;
  bbox.x2 = bbox.y2 = 0.5;
  bbox.z1 = bbox.z2 = 0.;

  gfs_init (&argc, &argv);
  advection_params.v = gfs_p;
  advection_params.fv = gfs_div;

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"nogts", no_argument, NULL, 'N'},
      {"nolimiter", no_argument, NULL, 'n'},
      {"smooth", no_argument, NULL, 's'},
      {"smooth1", no_argument, NULL, 'S'},
      {"translate", no_argument, NULL, 't'},
      {"periodic", no_argument, NULL, 'p'},
      {"error", no_argument, NULL, 'e'},
      {"right", required_argument, NULL, 'R'},
      {"left", required_argument, NULL, 'r'},
      {"box", required_argument, NULL, 'b'},
      {"levels", required_argument, NULL, 'l'},
      {"cfl", required_argument, NULL, 'c'},
      {"help", no_argument, NULL, 'h'},
      {"verbose", no_argument, NULL, 'v'},
      {"version", no_argument, NULL, 'V'}
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "hvVc:l:r:eR:b:ptsSnN",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "hvVc:l:r:eR:b:ptsSnN"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'N': /* nogts */
      gtsfiles = FALSE;
      break;
    case 'n': /* nolimiter */
      advection_params.gradient = gfs_center_gradient;
      break;
    case 's': /* smooth */
      smooth = smooth_gaussian;
      break;
    case 'S': /* smooth1 */
      smooth = smooth_periodic;
      break;
    case 't': /* translate */
      transform = translate;
      break;
    case 'p': /* periodic */
      periodic = TRUE;
      break;
    case 'e': /* error */
      display_error = TRUE;
      break;
    case 'R': /* right */
      rightrefine = atoi (optarg);
      break;
    case 'r': /* left */
      leftrefine = atoi (optarg);
      break;
    case 'b': /* box */
      boxrefine = atoi (optarg);
      break;
    case 'l': /* levels */
      maxlevel = atoi (optarg);
      break;
    case 'c': /* cfl */
      cfl = atof (optarg);
      break;
    case 'v': /* verbose */
      verbose = TRUE;
      break;
    case 'h': /* help */
      fprintf (stderr,
     "Usage: advection [OPTION] < FILE\n"
     "Advects a shape defined by volume fraction. FILE is a GTS file describing the interface.\n"
     "\n"
     "  -n    --nolimiter   do not use any limiter\n"
     "  -N    --nogts       do not output GTS files\n"
     "  -s    --smooth      advects smooth (gaussian) fraction field (do not use FILE)\n"
     "  -S    --smooth1     advects smooth (periodic) fraction field (do not use FILE)\n"
     "  -t    --translate   translation (default is rotation)\n"
     "  -p    --periodic    set periodic boundaries in both directions\n"
     "  -e    --error       display error relative to initial state\n"
     "  -R L  --right=L     refine the right half of the mesh using L extra levels\n"
     "  -r L  --left=L      refine the left half of the mesh using L extra levels\n"
     "  -b L  --box=L       refine a central box using L extra levels\n"
     "  -l L  --levels=L    maximum number of levels (default is 6)\n"
     "  -c V  --cfl=V       change CFL (default is 0.5)\n"
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
	       "advection: using %dD libgfs version %s\n"
	       "compiled with flags: %s\n",
	       FTT_DIMENSION,
	       GFS_VERSION,
	       GFS_COMPILATION_FLAGS);
      return 0; /* succes */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `advection --help' for more information.\n");
      return 1; /* failure */
    }
  }

  /* create cell tree with maxlevel */
  domain = GFS_DOMAIN (gts_graph_new (GTS_GRAPH_CLASS (gfs_domain_class ()), 
				     GTS_GNODE_CLASS (gfs_box_class ()),
				     gts_gedge_class ()));
  box = gfs_box_new (gfs_box_class ());
  box->root = ftt_cell_new ((FttCellInitFunc) gfs_cell_init, domain);
  gts_container_add (GTS_CONTAINER (domain), GTS_CONTAINEE (box));
  ftt_cell_refine (box->root, (FttCellRefineFunc) refine, &maxlevel, 
		   (FttCellInitFunc) gfs_cell_init, domain);

  if (leftrefine > 0) {
    maxlevel += leftrefine;
    ftt_cell_refine (box->root, (FttCellRefineFunc) refine_left, &maxlevel, 
		     (FttCellInitFunc) gfs_cell_init, domain);
  }
  else if (rightrefine > 0) {
    maxlevel += rightrefine;
    ftt_cell_refine (box->root, (FttCellRefineFunc) refine_right, &maxlevel, 
		     (FttCellInitFunc) gfs_cell_init, domain);
  }
  else if (boxrefine > 0) {
    maxlevel += boxrefine;
    ftt_cell_refine (box->root, (FttCellRefineFunc) refine_box, &maxlevel, 
		     (FttCellInitFunc) gfs_cell_init, domain);
  }
  
  /* set periodic boundary conditions */
  if (periodic) {
    ftt_cell_set_neighbor (box->root, box->root, FTT_RIGHT, 
			   (FttCellInitFunc) gfs_cell_init, domain);
    ftt_cell_set_neighbor (box->root, box->root, FTT_LEFT, 
			   (FttCellInitFunc) gfs_cell_init, domain);
    ftt_cell_set_neighbor (box->root, box->root, FTT_TOP, 
			   (FttCellInitFunc) gfs_cell_init, domain);
    ftt_cell_set_neighbor (box->root, box->root, FTT_BOTTOM, 
			   (FttCellInitFunc) gfs_cell_init, domain);
  }

  /* refine corners (eventually) */
  ftt_cell_refine (box->root, 
		   (FttCellRefineFunc) ftt_refine_corner, NULL,
		   (FttCellInitFunc) gfs_cell_init, domain);

  if (!smooth) {
    GtsSurface * interface;
    GNode * itree;
    gboolean is_open;
    GtsFile * fp = gts_file_new (stdin);

    /* read GTS surface */
    interface = gts_surface_new (gts_surface_class (),
				 gts_face_class (),
				 gts_edge_class (),
				 gts_vertex_class ());
    if (gts_surface_read (interface, fp)) {
      fprintf (stderr, 
	       "advection: invalid GTS file\n"
	       "%d:%d: %s\n",
	       fp->line, fp->pos, fp->error);
      return 1;
    }
    
    /* compute solid fractions using surface */
    itree = gts_bb_tree_surface (interface);
    is_open = gts_surface_volume (interface) < 0. ? TRUE : FALSE;
    gfs_cell_init_solid_fractions (box->root, interface, itree, is_open,
				   FALSE, NULL, NULL);
    if (leftrefine == 0 && rightrefine == 0 && boxrefine == 0)
      g_assert (gfs_cell_check_solid_fractions (box->root, 
						interface, itree, is_open));
    gts_bb_tree_destroy (itree, TRUE);

    /* save solid fraction */
    ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		       (FttCellTraverseFunc) copy_fraction, NULL);

    /* reset all cells to fluid */
    gfs_cell_fluid (box->root);
  }
  else
    ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) smooth, NULL);

  /* initialize advection velocities */
  {
    gdouble size = ftt_level_size (maxlevel);
    gdouble speed = cfl*size;

    if (transform == rotate) {
      twrite = PI/8./speed + 1;
      speed = PI/8./twrite;
    }
    else {
      twrite = 4.472135955/8./speed + 1;
      speed = 4.472135955/8./twrite;
    }
    tmax = 8*twrite;
    ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) transform, &speed);
    /* compute MAC velocities from centered velocities */
    ftt_face_traverse (box->root,  FTT_XYZ, 
		       FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, 
		       NULL);
    ftt_face_traverse (box->root, FTT_XYZ,
		       FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
 (FttFaceTraverseFunc) gfs_face_interpolated_normal_velocity, NULL);
  }

  if (display_error)
    ref = ftt_cell_copy (box->root, (FttCellCopyFunc) gfs_cell_copy, domain);

  ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) total_mass, &total_mass_before);
  if (verbose)
    fprintf (stderr, "total mass: %g\n", total_mass_before);

  for (t = 0; t < tmax; t++) {
    if (t % twrite == 0) {
      if (verbose) {
	stats = gfs_stats_variable (box->root, gfs_p, FTT_TRAVERSE_LEAFS, -1);
	fprintf (stderr, 
		 "fraction min: %g avg: %g| %g max: %g n: %7d\n",
		 stats.min, stats.mean, stats.stddev, stats.max, stats.n);
      }

      if (gtsfiles) {
	char fname[80];
	FILE * fp;

	sprintf (fname, "advected.%05d.gts", t);
	fp = fopen (fname, "wt");
	gfs_write_gts (domain, gfs_p, FTT_TRAVERSE_LEAFS, -1, &bbox, fp);
	fclose (fp);
	if (verbose)
	  fprintf (stderr, "advected.%05d.gts written\n", t);
      }

      if (display_error) {
	GfsNorm norm;
	static guint n = 0;

	ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			   (FttCellTraverseFunc) compute_error, ref);
	norm = gfs_norm_variable (box->root, gfs_res, FTT_TRAVERSE_LEAFS, -1);
	fprintf (stderr, 
	   "Error n: %d first: %g second: %g infty: %g bias: %g w: %g\n",
		 n++, norm.first, norm.second, norm.infty, norm.bias, norm.w);
	if (gtsfiles) {
	  char fname[80];
	  FILE * fp;
	  
	  sprintf (fname, "error.%05d.gts", t);
	  fp = fopen (fname, "wt");
	  gfs_write_gts (domain, gfs_res, FTT_TRAVERSE_LEAFS, -1, &bbox, fp);
	  fclose (fp);
	  fprintf (stderr, "error.%05d.gts written\n", t);
	}
      }
    }
    ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) gfs_cell_reset, 
		       advection_params.fv);
    ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) gfs_cell_advected_face_values,
		       &advection_params);
    ftt_face_traverse (box->root, FTT_XYZ,
		       FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttFaceTraverseFunc) gfs_face_advection_flux,
		       &advection_params);
    ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) update_fraction, NULL);
  }

  if (gtsfiles) {
    char fname[80];
    FILE * fp;
    
    sprintf (fname, "advected.%05d.gts", t);
    fp = fopen (fname, "wt");
    gfs_write_gts (domain, gfs_p, FTT_TRAVERSE_LEAFS, -1, &bbox, fp);
    fclose (fp);
    if (verbose)
      fprintf (stderr, "advected.%05d.gts written\n", t);
  }

  if (display_error) {
    GfsNorm norm;
    
    ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) compute_error, ref);
    norm = gfs_norm_variable (box->root, gfs_res, FTT_TRAVERSE_LEAFS, -1);
    fprintf (stderr, 
	     "Error n: 8 first: %g second: %g infty: %g bias: %g n: %g\n",
	     norm.first, norm.second, norm.infty, norm.bias, norm.w);
  }
  
  stats = gfs_stats_variable (box->root, gfs_p, FTT_TRAVERSE_LEAFS, -1);
  ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) total_mass, &total_mass_after);
  if (verbose)
    fprintf (stderr, 
	     "fraction min: %g avg: %g| %g max: %g n: %7d\n"
	     "total mass: %g\n",
	     stats.min, stats.mean, stats.stddev, stats.max, stats.n,
	     total_mass_after);

  if (fabs (total_mass_before - total_mass_after) > 1e-10 ||
      stats.min < -1e-2 || stats.max > 1. + 1e-2)
    return 1; /* failure */
  return 0; /* success */
}
