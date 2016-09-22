/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2010 National Institute of Water and
 * Atmospheric Research
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

#include "config.h"

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <getopt.h>
#include <gfs.h>

#include "render.h"

static void read_commands (GtsFile * fp, 
			   GList ** list, 
			   GfsGlViewParams * view, 
			   GfsSimulation * sim)
{
  while (fp->type == GTS_STRING)
    if (!strcmp (fp->token->str, "View"))
      gfs_gl_view_params_read (view, fp);
    else if (!strcmp (fp->token->str, "Clear")) {
      g_list_foreach (*list, (GFunc) gts_object_destroy, NULL);
      g_list_free (*list);
      *list = NULL;
      gts_file_next_token (fp);
    }
    else if (!strcmp (fp->token->str, "Echo")) {
      gts_file_next_token (fp);
      if (fp->type != '{') {
	gts_file_error (fp, "expecting an opening brace");
	break;
      }
      guint scope = fp->scope_max;
      gint c = gts_file_getc (fp);
      while (c != EOF && fp->scope > scope) {
	putchar (c);
	c = gts_file_getc (fp);
      }
      fflush (stdout);
      if (fp->scope != scope) {
	gts_file_error (fp, "parse error");
	break;
      }
      gts_file_next_token (fp);
    }
    else if (!strcmp (fp->token->str, "Save") || !strcmp (fp->token->str, "Append")) {
      GfsGl2PSParams p;
      GfsOutputFile * out = NULL;
      FILE * fptr = NULL;

      gts_file_next_token (fp);
      if (fp->type != GTS_STRING) {
	gts_file_error (fp, "expecting a string (filename)");
	break;
      }
      if (!strcmp (fp->token->str, "Save"))
	fptr = (!strcmp (fp->token->str, "stdout") ? stdout :
		!strcmp (fp->token->str, "stderr") ? stderr :
		fopen (fp->token->str, "w"));
      else { /* Append */
	if ((out = gfs_output_file_open (fp->token->str, "w")))
	  fptr = out->fp;
      }
      if (fptr == NULL) {
	gts_file_error (fp, "cannot open file `%s'", fp->token->str);
	break;
      }
      gts_file_next_token (fp);
      gfs_gl2ps_params_read (&p, fp);
      if (fp->type != GTS_ERROR)
	gfs_gl_osmesa_render (&p, sim, view, *list, fptr, FALSE);
      if (out) {
	/* Append mode, just free memory, do not close file */
	out->refcount++;
	gfs_output_file_close (out);
      }
      else if (fptr != stdout && fptr != stderr)
	fclose (fptr);
    }
    else { /* GfsGl objects */
      GfsGl * gl;

      if ((gl = gfs_gl_new_from_file (fp))) {
	gl->p = view;
	if (sim)
	  gfs_gl_set_simulation (gl, sim);
	*list = g_list_append (*list, gl);
      }
      else if (fp->type != GTS_ERROR) {
	gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
	break;
      }
    }
}

static GfsSimulation * read_simulation (GtsFile * fp, GList * list)
{
  GfsSimulation * sim = gfs_simulation_read (fp);
  
  if (sim == NULL)
    return NULL;
  
  gfs_simulation_init (sim);
  g_list_foreach (list, (GFunc) gfs_gl_set_simulation, sim);
  return sim;
}

int main (int argc, char * argv[])
{
  GfsGlViewParams view;
  GfsSimulation * sim = NULL;
  GList * list = NULL;
  int c = 0;

  /* initialize gfs */
  gfs_init (&argc, &argv);

  /* initialize gfsgl */
  gfs_gl_init ();

  /* options */
  while (c != EOF) {
    static struct option long_options[] = {
      {"survive-broken-pipe", no_argument, NULL, 's'},
      {"help", no_argument, NULL, 'h'},
      {"version", no_argument, NULL, 'V'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "hVs", long_options, &option_index))) {
    case 's':
      /* These options are ignored, they are here for command-line compatibility with
	 GfsView interactive version */
      break;
    case 'h': /* help */
      fprintf (stderr,
             "Usage: gfsview-batch [OPTION] FILE1 FILE2 ...\n"
	     "The Gerris flow solver visualisation tool (batch mode).\n"
	     "\n"
	     "  -h    --help                display this help and exit\n"
	     "  -V    --version             output version information and exit\n"
	     "\n"
	     "Reports bugs to %s\n",
	     FTT_MAINTAINER);
      return 0; /* success */
      break;
    case 'V': /* version */
      fprintf (stderr,
	       "gfsview-batch: using %dD libgfs version %s\n",
	       FTT_DIMENSION, VERSION);
      return 0; /* succes */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `gfsview-batch --help' for more information.\n");
      return 1; /* failure */
    }
  }

  gfs_gl_view_params_init (&view);
  
  /* read files on command line */
  for (c = optind; c < argc; c++) {
    FILE * fptr = gfs_gl_popen (argv[c]);
    GtsFile * fp;

    if (fptr == NULL) {
      fprintf (stderr, "gfsview-batch: cannot open file `%s'\n", argv[c]);
      return 1;
    }
    fp = gts_file_new (fptr);
    while (fp->type != GTS_ERROR && !feof (fptr)) {
      if (fp->type == '\n')
	gts_file_next_token (fp);
      else if (fp->type == GTS_INT || !strcmp (fp->token->str, "GModule")) {
	GfsSimulation * sim1 = read_simulation (fp, list);
	if (sim1) {
	  if (sim)
	    gts_object_destroy (GTS_OBJECT (sim));
	  sim = sim1;
	}
      }
      else if (fp->type == GTS_STRING)
	read_commands (fp, &list, &view, sim);
      else
	gts_file_error (fp, "expecting an integer got %d", fp->type);
    }
    if (fp->type == GTS_ERROR) {
      fprintf (stderr, "gfsview-batch: %s:%d:%d: %s\n", argv[c], fp->line, fp->pos, fp->error);
      return 1;
    }
    gts_file_destroy (fp);
    pclose (fptr);
  }

  /* wait for parameters/simulations on standard input */
  while (1) {
    GtsFile * fp = gts_file_new (stdin);

    while (fp->type != GTS_ERROR)
      if (feof (stdin))
	return 0;
      else if (fp->type == '\n')
	gts_file_next_token (fp);
      else if (fp->type == GTS_INT) {
	GfsSimulation * sim1 = read_simulation (fp, list);
	if (sim1) {
	  if (sim)
	    gts_object_destroy (GTS_OBJECT (sim));
	  sim = sim1;	  
	}
      }
      else if (fp->type == GTS_STRING)
	read_commands (fp, &list, &view, sim);
      else
	gts_file_error (fp, "expecting an integer got %d", fp->type);

    fprintf (stderr, "gfsview-batch: <stdin>:%d:%d: %s\n", fp->line, fp->pos, fp->error);
    return 1;
  }
  
  return 0;
}
