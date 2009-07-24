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
#include "solid.h"
#include "version.h"

static void set_box_pid (GfsBox * box, gint * pid)
{
  box->pid = *pid;
}

static void setup_binary_IO (GfsDomain * domain)
{
  /* make sure that all the variables are sent */
  g_slist_free (domain->variables_io);
  domain->variables_io = NULL;
  GSList * i = domain->variables;
  while (i) {
    if (GFS_VARIABLE1 (i->data)->name)
      domain->variables_io = g_slist_append (domain->variables_io, i->data);
    i = i->next;
  }
  domain->binary = TRUE;	
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
  gboolean profile = FALSE, macros = FALSE, one_box_per_pe = TRUE, bubble = FALSE, verbose = FALSE;
  gchar * m4_options = g_strdup ("-P");
  GPtrArray * events = g_ptr_array_new ();
  gint maxlevel = -2;

  gfs_init (&argc, &argv);

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"split", required_argument, NULL, 's'},
      {"pid", no_argument, NULL, 'i'},
      {"partition", required_argument, NULL, 'p'},
      {"profile", no_argument, NULL, 'P'},
      {"define", required_argument, NULL, 'D'},
      {"macros", no_argument, NULL, 'm'},
      {"data", no_argument, NULL, 'd'},
      {"event", required_argument, NULL, 'e'},
      {"bubble", required_argument, NULL, 'b'},
      {"verbose", no_argument, NULL, 'v'},
      {"help", no_argument, NULL, 'h'},
      {"version", no_argument, NULL, 'V'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "hVs:ip:PD:mde:b:v",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "hVs:ip:PD:mde:b:v"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'P': /* profile */
      profile = TRUE;
      break;
    case 'p': /* partition */
      npart = atoi (optarg);
      break;
    case 'b': /* "bubble" partition */
      npart = atoi (optarg);
      bubble = TRUE;
      break;
    case 's': /* split */
      split = atoi (optarg);
      break;
    case 'i': /* pid */
      one_box_per_pe = FALSE;
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
    case 'd': /* data */
      maxlevel = -1;
      break;
    case 'e': /* event */
      g_ptr_array_add (events, g_strdup (optarg));
      break;
    case 'v': /* verbose */
      verbose = TRUE;
      break;
    case 'h': /* help */
      gfs_error (0,
             "Usage: gerris [OPTION] FILE\n"
	     "The Gerris flow solver simulation engine.\n"
	     "\n"
	     "  -s N   --split=N     splits the domain N times and returns\n"
             "                       the corresponding simulation\n"
	     "  -i     --pid         keep box pids when splitting\n"
             "  -p N   --partition=N partition the domain in 2^N subdomains and returns\n" 
             "                       the corresponding simulation\n"
             "  -b N   --bubble=N    partition the domain in N subdomains and returns\n" 
             "                       the corresponding simulation\n"
	     "  -d     --data        when splitting or partitioning, output all data\n"
	     "  -P     --profile     profiles calls to boundary conditions\n"
#ifdef HAVE_M4
	     "  -m     --macros      Turn macros support on\n"
	     "  -DNAME               Defines NAME as a macro expanding to VALUE\n"
	     "  -DNAME=VALUE         (macro support is implicitly turned on)\n"
	     "         --define=NAME\n"
             "         --define=NAME=VALUE\n"
#endif /* have m4 */
	     "  -eEV   --event=EV    Evaluates GfsEvent EV and returns the simulation\n"
	     "  -v     --verbose     Display more messages\n"
	     "  -h     --help        display this help and exit\n"
	     "  -V     --version     output version information and exit\n"
	     "\n"
	     "Reports bugs to %s\n",
	     FTT_MAINTAINER);
      return 0; /* success */
      break;
    case 'V': /* version */
      gfs_error (0,
	       "gerris: using %dD libgfs version %s (%s)\n"
	       "  compiled with flags: %s\n"
	       "  MPI:          %s\n"
	       "  pkg-config:   %s\n"
	       "  m4 and awk:   %s\n"
	       "Copyright (C) 2001-2009 NIWA.\n"
	       "This is free software; see the source for copying conditions.  There is NO\n"
	       "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n",
	       FTT_DIMENSION,
	       GFS_VERSION,
	       GFS_BUILD_VERSION,
	       GFS_COMPILATION_FLAGS,
#ifdef HAVE_MPI
	       "yes",
#else
	       "no",
#endif
#ifdef HAVE_PKG_CONFIG
	       "yes",
#else
	       "no",
#endif
#ifdef HAVE_M4
	       "yes"
#else
	       "no"
#endif
	       );
      return 0; /* succes */
      break;
    case '?': /* wrong options */
      gfs_error (0, "Try `gerris --help' for more information.\n");
      return 1; /* failure */
    }
  }

  if (optind >= argc) { /* missing FILE */
    gfs_error (0, 
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
  }
  else { /* no macros */
    if (!strcmp (argv[optind], "-"))
      fptr = stdin;
    else
      fptr = fopen (argv[optind], "r");
  }
  g_free (m4_options);

  if (fptr == NULL) {
    gfs_error (-1, "gerris: unable to open file `%s'\n", argv[optind]);
    return 1;
  }

  fp = gts_file_new (fptr);
  if (!(simulation = gfs_simulation_read (fp))) {
    gfs_error (-1, 
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
    guint np = bubble ? npart : pow (2., npart);
    gfloat imbalance = 0.0;
    GSList * partition, * i;
    gint pid = 0;

    if (verbose)
      gts_graph_print_stats (GTS_GRAPH (simulation), stderr);
    if (gts_container_size (GTS_CONTAINER (simulation)) < np) {
      fprintf (stderr,
	       "gerris: the number of boxes in the domain to partition should be >= %d\n"
	       "Use option '-s' to split the domain first\n"
	       "Try `gerris --help' for more information.\n",
	       np);
      return 1;
    }
    if (bubble)
      partition = gts_graph_bubble_partition (GTS_GRAPH (simulation), npart, 100, 
					      verbose ? 
					      (GtsFunc) gts_graph_partition_print_stats : NULL, 
					      stderr);
    else
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
    if (verbose)
      gts_graph_partition_print_stats (partition, stderr);
    gts_graph_partition_destroy (partition);
    gfs_simulation_write (simulation, maxlevel, stdout);
    return 0;
  }

  domain = GFS_DOMAIN (simulation);
  if (split) {
    gfs_clock_start (domain->timer);
    gfs_simulation_refine (simulation);
    gfs_clock_stop (domain->timer);
    while (split) {
      gfs_domain_split (domain, one_box_per_pe);
      split--;
    }
    gfs_simulation_write (simulation, maxlevel, stdout);
    return 0;
  }

  if (events->len > 0) {
    GSList * l = NULL;
    guint i;
    for (i = 0; i < events->len; i++) {
      GtsFile * fp = gts_file_new_from_string (g_ptr_array_index (events, i));
      if (fp->type != GTS_STRING) {
	gfs_error (-1, 
		   "gerris: invalid event: '%s'\n"
		   "expecting a GfsEvent name\n",
		   (char *) g_ptr_array_index (events, i));
	return 1;
      }
      GtsObjectClass * klass = gfs_object_class_from_name (fp->token->str);
      if (klass == NULL) {
	gfs_error (-1, "gerris: unknown event class `%s'\n", fp->token->str);
	return 1;
      }
      if (!gts_object_class_is_from_class (klass, gfs_event_class ())) {
	gfs_error (-1, "gerris: class `%s' is not a GfsEvent\n", fp->token->str);
	return 1;
      }
      GtsObject * object = gts_object_new (klass);
      gfs_object_simulation_set (object, simulation);
      g_assert (klass->read);
      (* klass->read) (&object, fp);
      if (fp->type == GTS_ERROR) {
	gfs_error (-1,
		   "gerris: invalid event: '%s'\n"
		   "%d:%d: %s\n",
		   (char *) g_ptr_array_index (events, i),
		   fp->line, fp->pos, fp->error);
	return 1;
      }
      if (GFS_IS_ADAPT (object))
	gts_container_add (GTS_CONTAINER (simulation->adapts), GTS_CONTAINEE (object));
      else if (GFS_IS_SOLID (object))
	gts_container_add (GTS_CONTAINER (simulation->solids), GTS_CONTAINEE (object));
      gts_container_add (GTS_CONTAINER (simulation->events), GTS_CONTAINEE (object));
      l = g_slist_append (l, object);
      gts_file_destroy (fp);
    }
    gfs_clock_start (domain->timer);
    g_slist_foreach (l, (GFunc) gfs_event_do, simulation);    
    gfs_clock_stop (domain->timer);
    setup_binary_IO (domain);
    gfs_simulation_write (simulation, -1, stdout);
    return 0;
  }

  domain->profile_bc = profile;

  gfs_simulation_run (simulation);

  return 0;
}
