/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001 National Institute for Water and Atmospheric research
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
#include "init.h"

#if 1
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

static void refine_all (GfsBox * box, guint * maxlevel)
{
  ftt_cell_refine (box->root, (FttCellRefineFunc) refine, maxlevel,
		   (FttCellInitFunc) gfs_cell_init, NULL);
}

static void draw_cells (GfsBox * box, FILE * fp)
{
  FttDirection d;

  gfs_draw_cells (box->root, FTT_TRAVERSE_LEAFS, -1, fp);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d]))
      gfs_draw_cells (GFS_BOUNDARY (box->neighbor[d])->root, 
		     FTT_TRAVERSE_LEAFS, -1, fp);
}

int main (int argc, char * argv[])
{
  GfsDomain * domain;
  FILE * fptr;
  GtsFile * fp;
  char fname[80];
  guint maxlevel = 2;
  GtsRange range;
  GfsNorm norm;
  GfsVariable * c;

  gfs_init (&argc, &argv);

  if (argc != 2) {
    fprintf (stderr, "usage: read FILE\n");
    return 1;
  }
  fptr = fopen (argv[1], "rt");
  if (fptr == NULL) {
    fprintf (stderr, "read: unable to open file `%s'\n", argv[1]);
    return 1;
  }

  fp = gts_file_new (fptr);
  if ((domain = gfs_domain_read (fp)) == NULL) {
    fprintf (stderr, "read: %s:%d:%d: %s\n",
	     argv[1], fp->line, fp->pos, fp->error);
    return 1;
  }
  gts_file_destroy (fp);
  fclose (fptr);

  maxlevel += domain->pid;
  gts_container_foreach (GTS_CONTAINER (domain),
			 (GtsFunc) refine_all, &maxlevel);
  gfs_domain_match (domain);

#if 0
  sprintf (fname, "/tmp/read.%d", domain->pid);
  fp = fopen (fname, "wt");
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) draw_cells, fp);
  //  gts_graph_write (GTS_GRAPH (domain), fp);
  fclose (fp);
#endif

  c = gfs_variable_from_name (domain->variables, "C");
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, c);

  range = gfs_domain_stats_variable (domain, c, FTT_TRAVERSE_LEAFS, -1);
  norm = gfs_domain_norm_variable (domain, c, FTT_TRAVERSE_LEAFS, -1);

  return 0;
}
