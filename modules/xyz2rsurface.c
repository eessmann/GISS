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

int main (int argc, char** argv)
{
  int c = 0, pagesize = 1024;
  int verbose = 0;

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"pagesize", required_argument, NULL, 'p'},
      {"verbose", no_argument, NULL, 'v'},
      {"help", no_argument, NULL, 'h'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "p:vh",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "p:vh"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'p': /* pagesize */
      pagesize = atoi (optarg);
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
	       "  -p N  --pagesize=N  sets the pagesize in bytes (default is 1024)\n"
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
  double x[3];
  int count = 0;
  while (scanf ("%lf %lf %lf", &x[0], &x[1], &x[2]) == 3) {
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
  r_surface_close (rs);

  return 0.;
}
