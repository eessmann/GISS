/* Gerris - The GNU Flow Solver
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

#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>

#include "config.h"
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */

#include "ftt.h"
#include "rsurface.h"

typedef struct {
  int maxdepth;
  double * ratio;
  int * size;
} Params;

static int includes (RSurfaceRect RSTrect, Params * p, int depth)
{
  if (RSTrect[0].l == RSTrect[0].h && RSTrect[1].l == RSTrect[1].h)
    p->size[depth]++;
  return 0;
}

static int intersects (RSurfaceRect RSTrect, Params * p, int depth)
{
  double w = RSTrect[0].h - RSTrect[0].l, h = RSTrect[1].h - RSTrect[1].l;
  double ratio = 1e10;
  if (w > 0. && h > 0.)
    ratio = w > h ? w/h : h/w;
  p->ratio[depth] += ratio;
  p->size[depth]++;
  return (depth < p->maxdepth);
}

static int includes_true (RSurfaceRect RSTrect, Params * p, int depth)
{
  return 1;
}

int main (int argc, char** argv)
{
  int c = 0, pagesize = 2048;
  int randomize = 0, verbose = 0;

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"pagesize", required_argument, NULL, 'p'},
      {"randomize", no_argument, NULL, 'r'},
      {"verbose", no_argument, NULL, 'v'},
      {"help", no_argument, NULL, 'h'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "p:vhr",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "p:vhr"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'p': /* pagesize */
      pagesize = atoi (optarg);
      break;
    case 'r': /* randomize */
      randomize = 1;
      break;
    case 'v': /* verbose */
      verbose = 1;
      break;
    case 'h': /* help */
      fprintf (stderr,
	       "Usage: xyz2rsurface [OPTION] BASENAME\n"
	       "\n"
	       "Converts the x, y and z coordinates on standard input to an\n"
	       "R*-tree-indexed database suitable for use with the\n"
	       "GfsRefineTerrain object of Gerris.\n"
	       "\n"
	       "  -p N  --pagesize=N  sets the pagesize in bytes (default is 2048)\n"
	       "  -r    --randomize   randomize (shuffle) the input\n"
	       "  -v    --verbose     display progress bar\n"
	       "  -h    --help        display this help and exit\n"
	       "\n"
	       "Report bugs to %s\n",
	       FTT_MAINTAINER);
      return 0; /* success */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `xyz2rsurface --help' for more information.\n");
      return 1; /* failure */
    }
  }

  if (optind >= argc) { /* missing BASENAME */
    fprintf (stderr, 
	     "xyz2rsurface: missing BASENAME\n"
	     "Try `xyz2rsurface --help' for more information.\n");
    return 1; /* failure */
  }

  RSurface * rs = r_surface_open (argv[optind], "w", pagesize);
  if (rs == NULL) {
    fprintf (stderr, "xyz2rsurface: cannot open files `%s*'\n", argv[optind]);
    return 1;
  }  

  FILE * input;
  if (!randomize)
    input = stdin;
  else {
    input = popen ("awk '{ printf(\"%5d %s\\n\", int (rand()*2**16), $0); }' "
		   "| sort -T. -n -k 1,2 | cut -c7-", "r");
    if (input == NULL) {
      fprintf (stderr, "xyz2rsurface: could not open randomization pipe\n");
      perror ("");
      return 1;
    }
  }
  
  double x[3];
  int count = 0;
  while (fscanf (input, "%lf %lf %lf", &x[0], &x[1], &x[2]) == 3) {
    if (!r_surface_insert (rs, x, 0)) {
      fprintf (stderr, "\nxyz2rsurface: error inserting point #%d (%g,%g,%g)\n",
	       count, x[0], x[1], x[2]);
      return 1;
    }
    count++;
    if (verbose && (count % 1000) == 0)
      fprintf (stderr, "\rxyz2rsurface: %9d points inserted", count);
  }
  if (verbose) {
    fprintf (stderr, "\rxyz2rsurface: %9d points inserted\n", count);
    fprintf (stderr, "xyz2rsurface: updating...\n");
  }
  r_surface_update (rs);

  if (randomize) {
    int status = pclose (input);
    if (status == -1) {
      fprintf (stderr, "xyz2rsurface: failed to close randomization pipe\n");
      perror ("");
      return 1;
    }
    status = WEXITSTATUS (status);
    if (status != 0) {
      fprintf (stderr, "xyz2rsurface: the randomization pipe returned error %d\n", status);
      return 1;
    }
  }

  if (verbose) {
    RSurfaceRect rect = {{-0.5,-0.5},{0.5,0.5}};
    RSurfaceSum s;

    r_surface_sum_init (&s);
    r_surface_query_region_sum (rs, 
				(RSurfaceCheck) includes_true, (RSurfaceCheck) includes_true, NULL,
				rect, &s);
    fprintf (stderr,
	     "Height min: %g average: %g max: %g\n", 
	     s.Hmin, s.H0/s.n, s.Hmax);

    fprintf (stderr, "Bounding box aspect ratio and average # of entries:\n");
    Params p;
    p.maxdepth = r_surface_depth (rs);
    p.size = calloc (p.maxdepth + 1, sizeof (int));
    p.ratio = calloc (p.maxdepth + 1, sizeof (double));
    r_surface_sum_init (&s);
    r_surface_query_region_sum (rs, (RSurfaceCheck) includes, (RSurfaceCheck) intersects, &p,
				rect, &s);
    int i;
    for (i = 1; i <= p.maxdepth; i++)
      if (p.size[i] > 0) {
	fprintf (stderr, "level %d: %d\n", i, p.size[i]);
	if (i < p.maxdepth && p.size[i + 1] > 0)
	  fprintf (stderr, "\taverage aspect ratio: %g average # of entries: %g\n",
		   p.ratio[i]/p.size[i],
		   p.size[i + 1]/(double) p.size[i]);
	else
	  fputc ('\n', stderr);
      }
  }

  r_surface_close (rs);

  return 0.;
}
