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

static int includes_true (RSurfaceRect RSTrect)
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
    RSurfaceStats * s = r_surface_stats_new (rs, maxdepth);
    fprintf (stderr, "Level\tAverage # of entries\t Aspect ratio percentiles\n");
    fprintf (stderr, "     \t                    \t  50%%   75%%   95%%   100%%\n");
    int i;
    for (i = 1; i < s->nlevel; i++) {
      fprintf (stderr, "%5d\t", i + 1);
      if (i < s->nlevel - 1)
	fprintf (stderr, "      %4.2f\t\t", s->nentries[i + 1]/(double)s->nentries[i]);
      else
	fprintf (stderr, "         -\t\t");
      fprintf (stderr, "%5.1f %5.1f %5.1f %5.1f\n",
	       1./s->aspect[i][(int) (0.50*(s->nentries[i] - 1))],
	       1./s->aspect[i][(int) (0.75*(s->nentries[i] - 1))],
	       1./s->aspect[i][(int) (0.95*(s->nentries[i] - 1))],
	       1./s->aspect[i][(int) (1.00*(s->nentries[i] - 1))]);
    }
    r_surface_stats_free (s);
  }
  r_surface_close (rs);

  return 0.;
}
