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

#include "init.h"
#include "simulation.h"
#include "refine.h"
#include "output.h"
#include "adaptive.h"

int main (int argc, char * argv[])
{
  GfsSimulation * simulation;
  FILE * fptr;
  GtsFile * fp;
  int c = 0;

  gfs_init (&argc, &argv);

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"help", no_argument, NULL, 'h'},
      {"version", no_argument, NULL, 'V'}
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "hV",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "hV"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'h': /* help */
      fprintf (stderr,
             "Usage: gerris_dir DIR1 DIR2... FILE\n"
	     "The Gerris flow solver simulation engine (multiple directories version).\n"
	     "\n"
	     "  -h    --help        display this help and exit\n"
	     "  -V    --version     output version information and exit\n"
	     "\n"
	     "Reports bugs to %s\n",
	     FTT_MAINTAINER);
      return 0; /* success */
      break;
    case 'V': /* version */
      fprintf (stderr,
     "gerris_dir: using %dD libgfs version %s\n"
     "compiled with flags: %s\n"
     "sizeof (GfsStateVector): %d sizeof (FttCell): %d sizeof (FttOct): %d\n",
	       FTT_DIMENSION,
	       GFS_VERSION,
	       GFS_COMPILATION_FLAGS,
	       sizeof (GfsStateVector),
	       sizeof (FttCell),
	       sizeof (struct _FttOct));
      return 0; /* succes */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `gerris_dir --help' for more information.\n");
      return 1; /* failure */
    }
  }

  if (optind >= argc) { /* missing DIR */
    fprintf (stderr, 
	     "gerris_dir: missing DIR\n"
	     "Try `gerris_dir --help' for more information.\n");
    return 1; /* failure */
  }

  if (optind + 1 >= argc) { /* missing FILE */
    fprintf (stderr, 
	     "gerris_dir: missing FILE\n"
	     "Try `gerris_dir --help' for more information.\n");
    return 1; /* failure */
  }
  
  simulation = gfs_simulation_new (gfs_simulation_class ());
  if (GFS_DOMAIN (simulation)->pid < 0) {
    fprintf (stderr, 
	     "gerris_dir: needs MPI support\n"
	     "Try `gerris_dir --help' for more information.\n");
    return 1; /* failure */
  }
  if (GFS_DOMAIN (simulation)->pid + optind >= argc - 1) {
    fprintf (stderr, 
	     "gerris_dir: missing DIR for PE %d\n"
	     "Try `gerris_dir --help' for more information.\n",
	     GFS_DOMAIN (simulation)->pid);
    return 1; /* failure */
  }
  if (chdir (argv[GFS_DOMAIN (simulation)->pid + optind])) {
    fprintf (stderr, 
	     "gerris_dir: cannot change to directory `%s'\n",
	     argv[GFS_DOMAIN (simulation)->pid + optind]);
    perror ("gerris_dir");
    return 1; /* failure */
  }

  fptr = fopen (argv[argc - 1], "rt");
  if (fptr == NULL) {
    fprintf (stderr, "gerris_dir: unable to open file `%s'\n", argv[argc - 1]);
    return 1;
  }

  fp = gts_file_new (fptr);
  if (!(simulation = gfs_simulation_read (fp))) {
    fprintf (stderr, 
	     "gerris_dir: file `%s' is not a valid simulation file\n"
	     "%s:%d:%d: %s\n",
	     argv[argc - 1], argv[argc - 1],
	     fp->line, fp->pos, fp->error);
    return 1;
  }
  gts_file_destroy (fp);
  fclose (fptr);

  gfs_simulation_run (simulation);

  return 0;
}







