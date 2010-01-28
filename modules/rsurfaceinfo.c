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
  int c = 0, maxdepth = 0;

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"help", no_argument, NULL, 'h'},
      {"depth", required_argument, NULL, 'd'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "hd:",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "hd:"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'd': /* depth */
      maxdepth = atoi (optarg);
      break;
    case 'h': /* help */
      fprintf (stderr,
	       "Usage: rsurfaceinfo [OPTION] BASENAME\n"
	       "\n"
	       "Displays info about a R*-tree-indexed database.\n"
	       "\n"
	       "  -dVAL --depth=VAL   display detailed statistics for levels up to VAL\n"
	       "                      (setting VAL to -1 will display all levels)\n"
	       "  -h    --help        display this help and exit\n"
	       "\n"
	       "Report bugs to %s\n",
	       FTT_MAINTAINER);
      return 0; /* success */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `rsurfaceinfo --help' for more information.\n");
      return 1; /* failure */
    }
  }

  if (optind >= argc) { /* missing BASENAME */
    fprintf (stderr, 
	     "rsurfaceinfo: missing BASENAME\n"
	     "Try `rsurfaceinfo --help' for more information.\n");
    return 1; /* failure */
  }

  RSurface * rs = r_surface_open (argv[optind], "r", 0);
  if (rs == NULL) {
    fprintf (stderr, "rsurfaceinfo: cannot open database `%s*'\n", argv[optind]);
    return 1;
  }

  RSurfaceRect rect = {{-0.5,-0.5},{0.5,0.5}};
  RSurfaceSum s;

  r_surface_info (rs);
  r_surface_sum_init (&s);
  r_surface_query_region_sum (rs, 
			      (RSurfaceCheck) includes_true, (RSurfaceCheck) includes_true, NULL,
			      rect, &s);
  fprintf (stderr,
	   "Height min: %g average: %g max: %g\n", 
	   s.Hmin, s.H0/s.n, s.Hmax);

  if (maxdepth != 0) {
    fprintf (stderr, "Bounding box aspect ratio and average # of entries:\n");
    Params p;
    p.maxdepth = maxdepth < 0 ? r_surface_depth (rs) : maxdepth;
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
	  fprintf (stderr, "\taverage ratio: %g\taverage # of entries: %g\n",
		   p.ratio[i]/p.size[i],
		   p.size[i + 1]/(double) p.size[i]);
	else
	  fputc ('\n', stderr);
      }
  }
  r_surface_close (rs);

  return 0.;
}
