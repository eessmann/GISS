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
#include "version.h"

static void set_box_pid (GfsBox * box, gint * pid)
{
  box->pid = *pid;
}

int main (int argc, char * argv[])
{
  GfsSimulation * simulation;
  GfsDomain * domain;
  FILE * fptr;
  GtsFile * fp;
  int c = 0;
  guint split = 0;
  guint npart = 0;
  gboolean profile = FALSE, macros = FALSE;
  gchar * m4_options = g_strdup ("-P");

  gfs_init (&argc, &argv);

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"split", required_argument, NULL, 's'},
      {"partition", required_argument, NULL, 'p'},
      {"profile", no_argument, NULL, 'P'},
      {"define", required_argument, NULL, 'D'},
      {"macros", no_argument, NULL, 'm'},
      {"help", no_argument, NULL, 'h'},
      {"version", no_argument, NULL, 'V'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "hVs:p:PD:m",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "hVs:p:PD:m"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'P': /* profile */
      profile = TRUE;
      break;
    case 'p': /* partition */
      npart = atoi (optarg);
      break;
    case 's': /* split */
      split = atoi (optarg);
      break;
    case 'D': { /* define */
      gchar * tmp = g_strjoin (" ", m4_options, "-D", optarg, NULL);
      g_free (m4_options);
      m4_options = tmp;
      /* fall through */
    }
    case 'm': /* macros */
      macros = TRUE;
      break;
    case 'h': /* help */
      fprintf (stderr,
             "Usage: gerris [OPTION] FILE\n"
	     "The Gerris flow solver simulation engine.\n"
	     "\n"
	     "  -s N   --split=N     splits the domain N times and returns\n"
             "                       the corresponding simulation\n"
             "  -p N   --partition=N partition the domain in 2^N subdomains and returns\n" 
             "                       the corresponding simulation\n"
	     "  -P     --profile     profiles calls to boundary conditions\n"
#ifdef HAVE_M4
	     "  -m     --macros      Turn macros support on\n"
	     "  -DNAME               Defines NAME as a macro expanding to VALUE\n"
	     "  -DNAME=VALUE         (macro support is implicitly turned on)\n"
	     "         --define=NAME\n"
             "         --define=NAME=VALUE\n"
#endif /* have m4 */
	     "  -h    --help        display this help and exit\n"
	     "  -V    --version     output version information and exit\n"
	     "\n"
	     "Reports bugs to %s\n",
	     FTT_MAINTAINER);
      return 0; /* success */
      break;
    case 'V': /* version */
      fprintf (stderr,
     "gerris: using %dD libgfs version %s (%s)\n"
     "compiled with flags: %s\n"
     "sizeof (GfsStateVector): %d sizeof (FttCell): %d sizeof (FttOct): %d\n",
	       FTT_DIMENSION,
	       GFS_VERSION,
	       GFS_BUILD_VERSION,
	       GFS_COMPILATION_FLAGS,
	       sizeof (GfsStateVector),
	       sizeof (FttCell),
	       sizeof (struct _FttOct));
      return 0; /* succes */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `gerris --help' for more information.\n");
      return 1; /* failure */
    }
  }

  if (optind >= argc) { /* missing FILE */
    fprintf (stderr, 
	     "gerris: missing FILE\n"
	     "Try `gerris --help' for more information.\n");
    return 1; /* failure */
  }

  if (macros) {
    const gchar awk[] = "awk -f " GFS_MODULES_DIR "/m4.awk ";
    gchar * command;
    
    if (!strcmp (argv[optind], "-"))
      command = g_strjoin (NULL, awk, "| m4 ", m4_options, NULL);
    else
      command = g_strjoin (NULL, awk, argv[optind], " | m4 ", m4_options, NULL);
    fptr = popen (command, "r");
    g_free (command);
    g_free (m4_options);
  }
  else { /* no macros */
    if (!strcmp (argv[optind], "-"))
      fptr = stdin;
    else
      fptr = fopen (argv[optind], "r");
  }

  if (fptr == NULL) {
    fprintf (stderr, "gerris: unable to open file `%s'\n", argv[optind]);
    return 1;
  }

  fp = gts_file_new (fptr);
  if (!(simulation = gfs_simulation_read (fp))) {
    fprintf (stderr, 
	     "gerris: file `%s' is not a valid simulation file\n"
	     "%s:%d:%d: %s\n",
	     argv[optind], argv[optind],
	     fp->line, fp->pos, fp->error);
    return 1;
  }
  gts_file_destroy (fp);

  if (macros)
    pclose (fptr);
  else
    if (fptr != stdin)
      fclose (fptr);

  if (npart > 0) {
    guint nmin = 1000;
    guint mmax = 10000;
    guint ntry = 10000;
    gfloat imbalance = 0.0;
    GSList * partition, * i;
    gint pid = 0;

    gts_graph_print_stats (GTS_GRAPH (simulation), stderr);
    if (gts_container_size (GTS_CONTAINER (simulation)) < pow (2., npart)) {
      fprintf (stderr,
	       "gerris: the number of boxes in the domain to partition should be >= 2^N\n"
	       "Use option '-s' to split the domain first\n"
	       "Try `gerris --help' for more information.\n");
      return 1;
    }
    partition = gts_graph_recursive_bisection (GTS_WGRAPH (simulation),
					       npart, 
					       ntry, mmax, nmin, imbalance);
    i = partition;
    while (i) {
      gts_container_foreach (GTS_CONTAINER (i->data), 
			     (GtsFunc) set_box_pid, &pid);
      pid++;
      i = i->next;
    }
    gts_graph_partition_print_stats (partition, stderr);
    gts_graph_partition_destroy (partition);
    gfs_simulation_write (simulation, -2, stdout);
    return 0;
  }

  domain = GFS_DOMAIN (simulation);
  if (split) {
    gfs_simulation_refine (simulation);
    while (split) {
      gfs_domain_split (domain, TRUE);
      split--;
    }
    gfs_simulation_write (simulation, -2, stdout);
    return 0;
  }

  domain->profile_bc = profile;

  gfs_simulation_run (simulation);

  return 0;
}
