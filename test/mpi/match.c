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
#include <string.h>
#include <math.h>
#include "config.h"
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif /* HAVE_UNISTD_H */

#include "graphic.h"
#include "domain.h"
#ifdef HAVE_MPI
#  include "mpi_boundary.h"
#else  /* not HAVE_MPI */
#  error "need MPI support for this program"
#endif /* not HAVE_MPI */
#include "init.h"

#if 0
static gboolean refine (FttCell * cell, guint * maxlevel)
{
  if (ftt_cell_level (cell) < *maxlevel)
    return TRUE;
  return FALSE;
}
#else
static gboolean refine (FttCell * cell, guint * maxlevel)
{
  gdouble probability = 0.5;

  if (ftt_cell_level (cell) < *maxlevel && 
      (gdouble) rand () < probability*(gdouble) RAND_MAX)
    return TRUE;
  return FALSE;
}
#endif

int main (int argc, char * argv[])
{
  GfsDomain * domain;
  GfsBox * box;
  GfsBoundaryMpi * bright, * bleft;
  int rank, size;
  guint maxlevel = 6;
  FILE * fp;
  char fname[80];

  gfs_init (&argc, &argv);

  atexit ((void (*)(void)) MPI_Finalize);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  if (size != 2) {
    fprintf (stderr, "match: must be run with exactly 2 processes\n");
    return 1;
  }

  domain = GFS_DOMAIN (gts_graph_new (GTS_GRAPH_CLASS (gfs_domain_class ()), 
				     GTS_GNODE_CLASS (gfs_box_class ()),
				     gts_gedge_class ()));
  box = gfs_box_new (gfs_box_class ());
  gts_container_add (GTS_CONTAINER (domain), GTS_CONTAINEE (box));

  bright = gfs_boundary_mpi_new (gfs_boundary_mpi_class (), box,
				FTT_RIGHT, rank == 0 ? 1 : 0, 0);
  bleft = gfs_boundary_mpi_new (gfs_boundary_mpi_class (), box,
			       FTT_LEFT, rank == 0 ? 1 : 0, 0);

  maxlevel += rank;
  ftt_cell_refine (box->root, (FttCellRefineFunc) refine, &maxlevel,
		   (FttCellInitFunc) gfs_cell_init, NULL);

  if (rank == 1) {
    FttCellChildren child;

    ftt_cell_children (box->root, &child);
    ftt_cell_destroy (child.c[0], (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
    ftt_cell_destroy (child.c[1], (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
  }

  gfs_domain_match (domain);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, gfs_p);

  sprintf (fname, "/tmp/match.%d", rank);
  fp = fopen (fname, "wt");
  gfs_draw_cells (box->root, FTT_TRAVERSE_LEAFS, -1, fp);
  gfs_draw_cells (GFS_BOUNDARY (bright)->root, FTT_TRAVERSE_LEAFS, -1, fp);
  gfs_draw_cells (GFS_BOUNDARY (bleft)->root, FTT_TRAVERSE_LEAFS, -1, fp);
  fclose (fp);

  return 0;
}
