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
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <errno.h>
#include <math.h>
#include "output.h"
#include "graphic.h"
#include "adaptive.h"
#include "solid.h"
#include "ocean.h"

/* GfsOutput: object */

typedef struct _Format Format;

typedef enum {
  ITER,
  TIME,
  PID,
  NONE
} FormatType;

struct _Format {
  gchar * s;
  FormatType t;
};

static Format * format_new (gchar * s, guint len, 
			    FormatType t)
{
  Format * f = g_malloc (sizeof (Format));
  
  f->s = g_strndup (s, len);
  f->t = t;

  return f;
}

static void format_destroy (Format * f)
{
  g_free (f->s);
  g_free (f);
}

static gchar * format_string (GSList * list, 
			      gint pid, 
			      guint niter,
			      gdouble time)
{
  gchar * s = g_strdup ("");

  while (list) {
    Format * f = list->data;
    gchar * s1, * s2 = NULL;

    switch (f->t) {
    case NONE:
      s2 = g_strconcat (s, f->s, NULL);
      break;
    case PID:
      s1 = g_strdup_printf (f->s, pid);
      s2 = g_strconcat (s, s1, NULL);
      g_free (s1);
      break;
    case ITER:
      s1 = g_strdup_printf (f->s, niter);
      s2 = g_strconcat (s, s1, NULL);
      g_free (s1);
      break;
    case TIME:
      s1 = g_strdup_printf (f->s, time);
      s2 = g_strconcat (s, s1, NULL);
      g_free (s1);
      break;
    default:
      g_assert_not_reached ();
    }
    g_free (s);
    s = s2;
    list = list->next;
  }

  return s;
}

static void output_free (GfsOutput * output)
{
  if (output->format)
    g_free (output->format);
  output->format = NULL;
  g_slist_foreach (output->formats, (GFunc) format_destroy, NULL);
  g_slist_free (output->formats);
  output->formats = NULL;
}

static void gfs_output_destroy (GtsObject * object)
{
  GfsOutput * output = GFS_OUTPUT (object);

  if (output->file)
    gfs_output_file_close (output->file);
  output_free (output);

  (* GTS_OBJECT_CLASS (gfs_output_class ())->parent_class->destroy) 
    (object);
}

static gboolean gfs_output_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* gfs_event_class ()->event) (event, sim)) {
    GfsOutput * output = GFS_OUTPUT (event);
    gchar * fname;

    if (!output->dynamic) {
      if (output->file) {
	fflush (output->file->fp);
	output->first_call = FALSE;
      }
      else {
	if (output->format[0] == '{') { /* script */
	  GString * script;
	  gint status;
	  gchar pname[L_tmpnam], * c;
	  guint i = 1, len;

	  if (!tmpnam (pname)) {
	    g_warning ("cannot create temporary name");
	    return FALSE;
	  }
	  if (mkfifo (pname, S_IWUSR|S_IRUSR)) {
	    g_warning ("cannot create named pipe: %s", strerror (errno));
	    return FALSE;
	  }
	  script = g_string_new ("sh -c \"");
	  c = output->format; c++;
	  len = strlen (output->format);
	  while (*c != '\0' && ++i < len) {
	    switch (*c) {
	    case '$': case '"':
	      g_string_append_c (script, '\\');
	    default:
	      g_string_append_c (script, *c);
	    }
	    c++;
	  }
	  g_string_append (script, "\" < ");
	  g_string_append (script, pname);
	  g_string_append (script, " &");
	  status = system (script->str);
	  g_string_free (script, TRUE);
	  if (status != -1)
	    status = WEXITSTATUS (status);
	  if (status == -1 || status != 0) {
	    g_warning ("error while executing script");
	    unlink (pname);
	    return FALSE;
	  }
	  output->file = gfs_output_file_open (pname, "w");
	  unlink (pname);
	}
	else { /* standard file */
	  fname = format_string (output->formats,
				 GFS_DOMAIN (sim)->pid,
				 sim->time.i,
				 sim->time.t);
	  output->file = gfs_output_file_open (fname, sim->time.i > 0 ? "a" : "w");
	  if (output->file == NULL)
	    g_warning ("could not open file `%s'", fname);
	  g_free (fname);
	}
      }
      return (output->file != NULL);
    }

    if (output->file)
      gfs_output_file_close (output->file);
    fname = format_string (output->formats, 
			   GFS_DOMAIN (sim)->pid,
			   sim->time.i,
			   sim->time.t);
    output->file = gfs_output_file_open (fname, "w");
    if (output->file == NULL)
      g_warning ("could not open file `%s'", fname);
    g_free (fname);
    return (output->file != NULL);
  }
  return FALSE;
}

static void gfs_output_write (GtsObject * o, FILE * fp)
{
  GfsOutput * output = GFS_OUTPUT (o);

  (* GTS_OBJECT_CLASS (gfs_output_class ())->parent_class->write) (o, fp);

  if (output->format)
    fprintf (fp, " %s", output->format);
}

static gboolean char_in_string (char c, const char * s)
{
  while (*s != '\0')
    if (*(s++) == c)
      return TRUE;
  return FALSE;
}

static void gfs_output_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutput * output;

  (* GTS_OBJECT_CLASS (gfs_output_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  output = GFS_OUTPUT (*o);
  if (output->file)
    gfs_output_file_close (output->file);
  output->file = NULL;
  if (output->format)
    g_free (output->format);
  output->format = NULL;
  output->dynamic = FALSE;
  output->first_call = TRUE;

  if (fp->type == '{') {
    GString * format;
    guint scope;
    gint c;

    format = g_string_new ("{");
    scope = fp->scope_max;
    c = gts_file_getc (fp);
    while (c != EOF && fp->scope > scope) {
      g_string_append_c (format, c);
      c = gts_file_getc (fp);
    }
    if (fp->scope != scope) {
      gts_file_error (fp, "parse error");
      g_string_free (format, TRUE);
      return;
    }
    g_string_append_c (format, '}');
    output->format = format->str;
    g_string_free (format, FALSE);
    gts_file_next_token (fp);
  }
  else if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (format)");
    return;
  }
  else {
    gchar * c, * start, * fname, * fnamebak;
    FILE * fptr;
    guint len;

    output->format = g_strdup (fp->token->str);
    gts_file_next_token (fp);
    
    if (!strcmp (output->format, "stderr")) {
      output->file = gfs_output_file_open ("stderr", "w");
      return;
    }
    
    if (!strcmp (output->format, "stdout")) {
      output->file = gfs_output_file_open ("stdout", "w");
      return;
    }
    
    start = c = output->format;
    while (*c != '\0') {
      if (*c == '%') {
	gchar * startf = c, * prev = c;
	
	len = GPOINTER_TO_UINT (startf) -  GPOINTER_TO_UINT (start);
	if (len > 0)
	  output->formats = g_slist_prepend (output->formats,
					     format_new (start, len, NONE));
	
	len = 1;
	c++;
	while (*c != '\0' && !char_in_string (*c, "diouxXeEfFgGaAcsCSpn%")) {
	  prev = c;
	  c++;
	  len++;
	}
	len++;
	if (*c == '%')
	  output->formats = g_slist_prepend (output->formats,
					     format_new ("%", 1, NONE));
	else if (char_in_string (*c, "diouxXc")) {
	  if (*prev == 'l') {
	    output->formats = g_slist_prepend (output->formats,
					       format_new (startf, len, ITER));
	    output->dynamic = TRUE;
	  }
	  else
	    output->formats = g_slist_prepend (output->formats,
					       format_new (startf, len, PID));
	}
	else if (char_in_string (*c, "eEfFgGaA")) {
	  output->formats = g_slist_prepend (output->formats,
					     format_new (startf, len, TIME));
	  output->dynamic = TRUE;
	}
	else {
	  gts_file_error (fp, 
			  "unknown conversion specifier `%c' of format `%s'",
			  *c, output->format);
	  output_free (output);
	  return;
	}
	start = c;
	start++;
      }
      c++;
    }
    len = GPOINTER_TO_UINT (c) -  GPOINTER_TO_UINT (start);
    if (len > 0)
      output->formats = g_slist_prepend (output->formats,
					 format_new (start, len, NONE));
    output->formats = g_slist_reverse (output->formats);
    
    fname = format_string (output->formats, -1, 0, 0.);
    fnamebak = g_strconcat (fname, "~", NULL);
    g_free (fname);
    fptr = fopen (fnamebak, "w");
    if (fptr == NULL) {
      gts_file_error (fp, "cannot open file specified by format `%s'",
		      output->format);
      g_free (fnamebak);
      output_free (output);
      return;
    }
    fclose (fptr);
    remove (fnamebak);
    g_free (fnamebak);
  }
}

static void gfs_output_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_event;

  GTS_OBJECT_CLASS (klass)->write = gfs_output_write;
  GTS_OBJECT_CLASS (klass)->read = gfs_output_read;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_output_destroy;
}

static void gfs_output_init (GfsOutput * object)
{
  object->file = NULL;
  object->format = NULL;
  object->formats = NULL;
  object->dynamic = FALSE;
  object->first_call = TRUE;
}

GfsOutputClass * gfs_output_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_info = {
      "GfsOutput",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_class_init,
      (GtsObjectInitFunc) gfs_output_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_output_info);
  }

  return klass;
}

/**
 * gfs_output_mute:
 * @output: a #GfsOutput.
 *
 * "Mutes" the output defined by @output, the event associated with
 * @output still takes place but the output itself is redirected to
 * /dev/null.
 */
void gfs_output_mute (GfsOutput * output)
{
  g_return_if_fail (output != NULL);

  output->dynamic = FALSE;
  if (output->file)
    gfs_output_file_close (output->file);
  output->file = gfs_output_file_open ("/dev/null", "w");
}

static GHashTable * gfs_output_files = NULL;

/**
 * gfs_output_file_open:
 * @name: the name of the file to open.
 * @mode: the fopen mode.
 *
 * Checks whether @name has already been opened. If it has, its
 * reference count is incremented and the corresponding #GfsOutputFile
 * is returned. If it has not, it is created and opened for writing.
 *
 * Returns: the #GfsOutputFile of file @name.  
 */
GfsOutputFile * gfs_output_file_open (const gchar * name, const gchar * mode)
{
  GfsOutputFile * file;
  FILE * fp;

  g_return_val_if_fail (name != NULL, NULL);

  if (gfs_output_files == NULL) {
    gfs_output_files = g_hash_table_new (g_str_hash, g_str_equal);
    file = g_malloc (sizeof (GfsOutputFile));
    file->refcount = 2;
    file->name = g_strdup ("stderr");
    file->fp = stderr;
    g_hash_table_insert (gfs_output_files, file->name, file);
    file = g_malloc (sizeof (GfsOutputFile));
    file->refcount = 2;
    file->name = g_strdup ("stdout");
    file->fp = stdout;
    g_hash_table_insert (gfs_output_files, file->name, file);
  }

  if ((file = g_hash_table_lookup (gfs_output_files, name))) {
    file->refcount++;
    return file;
  }

  fp = fopen (name, mode);
  if (fp == NULL)
    return NULL;

  file = g_malloc (sizeof (GfsOutputFile));
  file->refcount = 1;
  file->name = g_strdup (name);
  file->fp = fp;
  g_hash_table_insert (gfs_output_files, file->name, file);

  return file;  
}

/**
 * gfs_output_file_close:
 * @file: a #GfsOutputFile.
 * 
 * Decreases the reference count of @file. If it reaches zero the file
 * corresponding to @file is closed and @file is freed.
 */
void gfs_output_file_close (GfsOutputFile * file)
{
  g_return_if_fail (file);

  file->refcount--;
  if (file->refcount == 0) {
    g_hash_table_remove (gfs_output_files, file->name);
    fclose (file->fp);
    g_free (file->name);
    g_free (file);
  }
}

/* GfsOutputTime: Object */

static gboolean time_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "step: %7u t: %15.8f dt: %13.6e\n",
	     sim->time.i, sim->time.t, 
	     sim->advection_params.dt);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_time_class_init (GfsEventClass * klass)
{
  klass->event = time_event;
}

GfsOutputClass * gfs_output_time_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_time_info = {
      "GfsOutputTime",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_time_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_time_info);
  }

  return klass;
}

/* GfsOutputProgress: Object */

static gboolean progress_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    gdouble idone = sim->time.i/(gdouble) sim->time.iend;
    gdouble tdone = sim->time.t/sim->time.end;

    if (idone > tdone) tdone = idone;
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "\r%3.0f%% done",
	     100.*tdone);
    if (tdone > 0.) {
      gdouble remaining = GFS_DOMAIN (sim)->timestep.sum*(1. - tdone)/tdone;
      gdouble hours = floor (remaining/3600.);
      gdouble mins = floor ((remaining - 3600.*hours)/60.);
      gdouble secs = floor (remaining - 3600.*hours - 60.*mins);
      fprintf (GFS_OUTPUT (event)->file->fp,
	       ", %02.0f:%02.0f:%02.0f remaining ",
	       hours, mins, secs);
    }
    if (tdone == 1.)
      fputc ('\n', GFS_OUTPUT (event)->file->fp);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_progress_class_init (GfsEventClass * klass)
{
  klass->event = progress_event;
}

GfsOutputClass * gfs_output_progress_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_progress_info = {
      "GfsOutputProgress",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_progress_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_progress_info);
  }

  return klass;
}

/* GfsOutputProjectionStats: Object */

static gdouble rate (gdouble a, gdouble b, guint n)
{
  if (a > 0. && b > 0. && n > 0)
    return exp (log (b/a)/n);
  return 0.;
}

static void multilevel_stats_write (GfsMultilevelParams * p,
				    FILE * fp)
{
  fprintf (fp,
	   "    niter: %4d\n"
	   "    residual.bias:   % 10.3e % 10.3e\n"
	   "    residual.first:  % 10.3e % 10.3e %6.2g\n"
	   "    residual.second: % 10.3e % 10.3e %6.2g\n"
	   "    residual.infty:  % 10.3e % 10.3e %6.2g\n",
	   p->niter,
	   p->residual_before.bias,
	   p->residual.bias,
	   p->residual_before.first,
	   p->residual.first,
	   rate (p->residual.first,
		 p->residual_before.first,
		 p->niter),
	   p->residual_before.second,
	   p->residual.second,
	   rate (p->residual.second,
		 p->residual_before.second,
		 p->niter),
	   p->residual_before.infty,
	   p->residual.infty,
	   rate (p->residual.infty,
		 p->residual_before.infty,
		 p->niter));
}

static gboolean projection_stats_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    FILE * fp = GFS_OUTPUT (event)->file->fp;

    if (sim->projection_params.niter > 0) {
      fprintf (fp, "MAC projection        before     after       rate\n");
      multilevel_stats_write (&sim->projection_params, fp);
    }
    fprintf (fp, "Approximate projection\n");
    multilevel_stats_write (&sim->approx_projection_params, fp);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_projection_stats_class_init (GfsEventClass * klass)
{
  klass->event = projection_stats_event;
}

GfsOutputClass * gfs_output_projection_stats_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_projection_stats_info = {
      "GfsOutputProjectionStats",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_projection_stats_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_projection_stats_info);
  }

  return klass;
}

/* GfsOutputDiffusionStats: Object */

static gboolean diffusion_stats_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    FILE * fp = GFS_OUTPUT (event)->file->fp;

    fprintf (fp, "Velocity diffusion    before     after       rate\n");
    multilevel_stats_write (&sim->diffusion_params, fp);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_diffusion_stats_class_init (GfsEventClass * klass)
{
  klass->event = diffusion_stats_event;
}

GfsOutputClass * gfs_output_diffusion_stats_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_diffusion_stats_info = {
      "GfsOutputDiffusionStats",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_diffusion_stats_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_diffusion_stats_info);
  }

  return klass;
}

/* GfsOutputSolidStats: Object */

static gboolean gfs_output_solid_stats_event (GfsEvent * event, 
					     GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_solid_stats_class ())->parent_class)->event) (event, sim)) {
    GtsRange stats = gfs_domain_stats_solid (GFS_DOMAIN (sim));
    GtsRange ma, mn;

    gfs_domain_stats_merged (GFS_DOMAIN (sim), &ma, &mn);
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "Solid volume fraction\n"
	     "    min: %10.3e avg: %10.3e | %10.3e max: %10.3e n: %10d\n"
	     "Total merged solid volume fraction\n"
	     "    min: %10.3e avg: %10.3e | %10.3e max: %10.3e n: %10d\n"
	     "Number of cells merged per merged cell\n"
	     "    min: %10.0f avg: %10.3f | %10.3f max: %10.0f n: %10d\n",
	     stats.min, stats.mean, stats.stddev, stats.max, stats.n,
	     ma.min, ma.mean, ma.stddev, ma.max, ma.n,
	     mn.min, mn.mean, mn.stddev, mn.max, mn.n);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_solid_stats_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_solid_stats_event;
}

GfsOutputClass * gfs_output_solid_stats_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_solid_stats_info = {
      "GfsOutputSolidStats",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_solid_stats_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_solid_stats_info);
  }

  return klass;
}

/* GfsOutputAdaptStats: Object */

static gboolean gfs_output_adapt_stats_event (GfsEvent * event, 
					      GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_adapt_stats_class ())->parent_class)->event) (event, sim)) {
    gfs_adapt_stats_update (&sim->adapts_stats);
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "Adaptive mesh refinement statistics\n"
	     "  Cells removed\n"
	     "    min: %10.0f avg: %10.3f | %10.3f max: %10.0f n: %10d\n"
	     "  Cells created\n"
	     "    min: %10.0f avg: %10.3f | %10.3f max: %10.0f n: %10d\n"
	     "  Maximum cost\n"
	     "    min: %10.3e avg: %10.3e | %10.3e max: %10.3e n: %10d\n"
	     "  Number of cells\n"
	     "    min: %10.0f avg: %10.3f | %10.3f max: %10.0f n: %10d\n",
	     sim->adapts_stats.removed.min,
	     sim->adapts_stats.removed.mean,
	     sim->adapts_stats.removed.stddev,
	     sim->adapts_stats.removed.max,
	     sim->adapts_stats.removed.n,
	     sim->adapts_stats.created.min,
	     sim->adapts_stats.created.mean,
	     sim->adapts_stats.created.stddev,
	     sim->adapts_stats.created.max,
	     sim->adapts_stats.created.n,
	     sim->adapts_stats.cmax.min,
	     sim->adapts_stats.cmax.mean,
	     sim->adapts_stats.cmax.stddev,
	     sim->adapts_stats.cmax.max,
	     sim->adapts_stats.cmax.n,
	     sim->adapts_stats.ncells.min,
	     sim->adapts_stats.ncells.mean,
	     sim->adapts_stats.ncells.stddev,
	     sim->adapts_stats.ncells.max,
	     sim->adapts_stats.ncells.n);
    gfs_adapt_stats_init (&sim->adapts_stats);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_adapt_stats_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_adapt_stats_event;
}

GfsOutputClass * gfs_output_adapt_stats_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_adapt_stats_info = {
      "GfsOutputAdaptStats",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_adapt_stats_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_adapt_stats_info);
  }

  return klass;
}

/* GfsOutputTiming: Object */

static void timing_print (GtsRange * r, gdouble total, FILE * fp)
{
  fprintf (fp, 
	   "      min: %9.3f avg: %9.3f (%4.1f%%) | %7.3f max: %9.3f\n",
	   r->min,
	   r->mean, total > 0. ? 100.*r->sum/total : 0.,
	   r->stddev, 
	   r->max);	   
}

static void timing_bc_print (GtsRange * r, gdouble total, guint n, FILE * fp)
{
  fprintf (fp,
	   "      calls: %7d avg: %9.3f (%4.1f%%)\n",
	   r->n,
	   n > 0 ? r->sum/n : 0.,
	   total > 0. ? 100.*r->sum/total : 0.);
}

static void timer_print (gchar * name, GfsTimer * t, gpointer * data)
{
  FILE * fp = data[0];
  GfsDomain * domain = data[1];

  fprintf (fp, "  %s:\n", name);
  timing_print (&t->r, domain->timestep.sum, fp);
}

static gboolean timing_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    FILE * fp = GFS_OUTPUT (event)->file->fp;

    if (domain->timestep.mean > 0.) {
      gpointer data[2];

      fprintf (fp,
	       "Timing summary: %u timesteps %.0f node.timestep/s\n"
	       "  timestep:\n"
	       "      min: %9.3f avg: %9.3f         | %7.3f max: %9.3f\n"
               "  domain size:\n"
	       "      min: %9.0f avg: %9.0f         | %7.0f max: %9.0f\n",
	       domain->timestep.n,
	       domain->size.mean/domain->timestep.mean,
	       domain->timestep.min,
	       domain->timestep.mean,
	       domain->timestep.stddev, 
	       domain->timestep.max,
	       domain->size.min,
	       domain->size.mean,
	       domain->size.stddev, 
	       domain->size.max);
      data[0] = fp;
      data[1] = domain;
      g_hash_table_foreach (domain->timers, (GHFunc) timer_print, data);
      if (domain->mpi_messages.n > 0)
	fprintf (fp,
		 "Message passing summary\n"
		 "  n: %10d size: %10.0f bytes\n",
		 domain->mpi_messages.n,
		 domain->mpi_messages.sum);
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_timing_class_init (GfsEventClass * klass)
{
  klass->event = timing_event;
}

GfsOutputClass * gfs_output_timing_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_timing_info = {
      "GfsOutputTiming",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_timing_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_timing_info);
  }

  return klass;
}

/* GfsOutputBalance: Object */

static gboolean gfs_output_balance_event (GfsEvent * event, 
					  GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_balance_class ())->parent_class)->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    FILE * fp = GFS_OUTPUT (event)->file->fp;
    GtsRange size, boundary, mpiwait;
    
    gfs_domain_stats_balance (domain, &size, &boundary, &mpiwait);
    fprintf (fp, 
	     "Balance summary: %u PE\n"
	     "  domain size:\n"
	     "      min: %9.0f avg: %9.0f         | %7.0f max: %9.0f\n",
	     size.n,
	     size.min, size.mean, size.stddev, size.max);
    if (boundary.max > 0.)
      fprintf (fp, 
	       "  parallel boundary size:\n"
	       "      min: %9.0f avg: %9.0f         | %7.0f max: %9.0f\n",
	       boundary.min, boundary.mean, boundary.stddev, boundary.max);
    if (mpiwait.max > 0.)
      fprintf (fp,
	       "  average timestep MPI wait time:\n"
	       "      min: %9.3f avg: %9.3f         | %7.3f max: %9.3f\n",
	       mpiwait.min, mpiwait.mean, mpiwait.stddev, mpiwait.max);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_balance_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_balance_event;
}

GfsOutputClass * gfs_output_balance_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_balance_info = {
      "GfsOutputBalance",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_balance_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_balance_info);
  }

  return klass;
}

/* GfsOutputSolidForce: Object */

static gboolean gfs_output_solid_force_event (GfsEvent * event, 
					      GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_solid_force_class ())->parent_class)->event)
      (event, sim) &&
      sim->advection_params.dt > 0.) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    FILE * fp = GFS_OUTPUT (event)->file->fp;
    FttVector pf, vf;

    if (GFS_OUTPUT (event)->first_call)
      fputs ("# 1: T (2,3,4): Pressure force (5,6,7): Viscous force\n", fp);
    
    gfs_domain_solid_force (domain, &pf, &vf);
    fprintf (fp, "%g %g %g %g %g %g %g\n",
	     sim->time.t,
	     pf.x, pf.y, pf.z,
	     vf.x, vf.y, vf.z);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_solid_force_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_solid_force_event;
}

GfsOutputClass * gfs_output_solid_force_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_solid_force_info = {
      "GfsOutputSolidForce",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_solid_force_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_solid_force_info);
  }

  return klass;
}

/* GfsOutputLocation: Object */

static void gfs_output_location_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputLocation * l = GFS_OUTPUT_LOCATION (*o);

  if (GTS_OBJECT_CLASS (gfs_output_location_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_output_location_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.x)");
    return;
  }
  l->p.x = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.y)");
    return;
  }
  l->p.y = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.z)");
    return;
  }
  l->p.z = atof (fp->token->str);
  gts_file_next_token (fp);  
}

static void gfs_output_location_write (GtsObject * o, FILE * fp)
{
  GfsOutputLocation * l = GFS_OUTPUT_LOCATION (o);

  if (GTS_OBJECT_CLASS (gfs_output_location_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_output_location_class ())->parent_class->write) 
      (o, fp);
  fprintf (fp, " %g %g %g", l->p.x, l->p.y, l->p.z);
}

static gboolean gfs_output_location_event (GfsEvent * event, 
					   GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_location_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsOutputLocation * location = GFS_OUTPUT_LOCATION (event);
    FttCell * cell = gfs_domain_locate (domain, location->p, -1);

    if (GFS_OUTPUT (event)->first_call) {
      FILE * fp = GFS_OUTPUT (event)->file->fp;
      GfsVariable * v = domain->variables;
      guint nv = 5;

      fputs ("# 1:T 2:X 3:Y 4:Z", fp);
      while (v) {
	if (v->name)
	  fprintf (fp, " %d:%s", nv++, v->name);
	v = v->next;
      }
      fputc ('\n', fp);
    }
    if (cell != NULL) {
      FILE * fp = GFS_OUTPUT (event)->file->fp;
      GfsVariable * v = domain->variables;
   
      fprintf (fp, "%g %g %g %g", 
	       sim->time.t,
	       location->p.x, location->p.y, location->p.z);
      while (v) {
	if (v->name)
	  fprintf (fp, " %g", gfs_interpolate (cell, location->p, v->i));
	v = v->next;
      }
      fputc ('\n', fp);
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_location_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_location_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_output_location_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_output_location_write;
}

static void gfs_output_location_init (GfsOutputLocation * object)
{
  object->p.x = object->p.y = object->p.z = 0.;
}

GfsOutputClass * gfs_output_location_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_location_info = {
      "GfsOutputLocation",
      sizeof (GfsOutputLocation),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_location_class_init,
      (GtsObjectInitFunc) gfs_output_location_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_location_info);
  }

  return klass;
}

/* GfsOutputSimulation: Object */

static void output_simulation_destroy (GtsObject * object)
{
  GfsOutputSimulation * output = GFS_OUTPUT_SIMULATION (object);

  gfs_variable_list_destroy (output->var);

  (* GTS_OBJECT_CLASS (gfs_output_simulation_class ())->parent_class->destroy) 
    (object);
}

static gboolean output_simulation_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsVariable * var = domain->variables_io;
    gboolean binary = domain->binary;

    domain->variables_io = GFS_OUTPUT_SIMULATION (event)->var;
    domain->binary = GFS_OUTPUT_SIMULATION (event)->binary;
    gfs_simulation_write (sim,
			  GFS_OUTPUT_SIMULATION (event)->max_depth,
			  GFS_OUTPUT (event)->file->fp);
    domain->variables_io = var;
    domain->binary = binary;
    fflush (GFS_OUTPUT (event)->file->fp);
    return TRUE;
  }
  return FALSE;
}

static void output_simulation_write (GtsObject * o, FILE * fp)
{
  GfsOutputSimulation * output = GFS_OUTPUT_SIMULATION (o);
  GfsVariable * v = output->var;

  (* GTS_OBJECT_CLASS (gfs_output_simulation_class ())->parent_class->write) (o, fp);

  fputs (" {", fp);
  if (output->max_depth != -1)
    fprintf (fp, " depth = %d", output->max_depth);
  if (v != NULL) {
    fprintf (fp, " variables = %s", v->name);
    v = v->next;
    while (v) {
      if (v->name)
	fprintf (fp, ",%s", v->name);
      v = v->next;
    }
  }
  if (output->binary)
    fputs (" binary = 1", fp);
  fputs (" }", fp);
}

static void output_simulation_read (GtsObject ** o, GtsFile * fp)
{
  GtsFileVariable var[] = {
    {GTS_INT,    "depth",     TRUE},
    {GTS_STRING, "variables", TRUE},
    {GTS_INT,    "binary",    TRUE},
    {GTS_NONE}
  };
  gchar * variables = NULL;
  GfsOutputSimulation * output = GFS_OUTPUT_SIMULATION (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (output));

  (* GTS_OBJECT_CLASS (gfs_output_simulation_class ())->parent_class->read) 
    (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  var[0].data = &output->max_depth;
  var[1].data = &variables;
  var[2].data = &output->binary;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR) {
    g_free (variables);
    return;
  }

  if (variables != NULL) {
    gchar * error = NULL;
    GfsVariable * vars = gfs_variables_from_list (domain->variables, 
						  variables, &error);

    if (vars == NULL) {
      gts_file_variable_error (fp, var, "variables",
			       "unknown variable `%s'", error);
      g_free (variables);
      return;
    }
    if (output->var)
      gfs_variable_list_destroy (output->var);
    output->var = vars;
    g_free (variables);
  }
  else if (output->var == NULL)
    output->var = gfs_variable_list_copy (domain->variables, 
					  GTS_OBJECT (domain));
}

static void gfs_output_simulation_class_init (GfsEventClass * klass)
{
  klass->event = output_simulation_event;
  GTS_OBJECT_CLASS (klass)->destroy = output_simulation_destroy;
  GTS_OBJECT_CLASS (klass)->read = output_simulation_read;
  GTS_OBJECT_CLASS (klass)->write = output_simulation_write;
}

static void gfs_output_simulation_init (GfsOutputSimulation * object)
{
  object->max_depth = -1;
  object->var = NULL;
}

GfsOutputClass * gfs_output_simulation_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_simulation_info = {
      "GfsOutputSimulation",
      sizeof (GfsOutputSimulation),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_simulation_class_init,
      (GtsObjectInitFunc) gfs_output_simulation_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_simulation_info);
  }

  return klass;
}

/* GfsOutputBoundaries: Object */

static gboolean output_boundaries_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    FILE * fp = GFS_OUTPUT (event)->file->fp;
    
    gfs_draw_refined_boundaries (domain, fp);
    gfs_draw_solid_boundaries (domain, fp);
    gfs_draw_boundary_conditions (domain, fp);
    fflush (fp);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_boundaries_class_init (GfsEventClass * klass)
{
  klass->event = output_boundaries_event;
}

GfsOutputClass * gfs_output_boundaries_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_boundaries_info = {
      "GfsOutputBoundaries",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_boundaries_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_boundaries_info);
  }

  return klass;
}

/* GfsOutputScalar: Object */

static void gfs_output_scalar_destroy (GtsObject * o)
{
  if (GFS_OUTPUT_SCALAR (o)->box)
    gts_object_destroy (GTS_OBJECT (GFS_OUTPUT_SCALAR (o)->box));
  (* GTS_OBJECT_CLASS (gfs_output_scalar_class ())->parent_class->destroy) (o);
}

static void gfs_output_scalar_read (GtsObject ** o, GtsFile * fp)
{
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "min",      TRUE},
    {GTS_DOUBLE, "max",      TRUE},
    {GTS_STRING, "v",        TRUE},
    {GTS_INT,    "maxlevel", TRUE},
    {GTS_STRING, "box",      TRUE},
    {GTS_NONE}
  };
  GfsOutputScalar * output;
  GfsDomain * domain;
  gchar * vname = NULL, * box = NULL;

  if (GTS_OBJECT_CLASS (gfs_output_scalar_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_output_scalar_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  output = GFS_OUTPUT_SCALAR (*o);
  domain = GFS_DOMAIN (gfs_object_simulation (output));
  output->autoscale = TRUE;

  var[0].data = &output->min;
  var[1].data = &output->max;
  var[2].data = &vname;
  var[3].data = &output->maxlevel;
  var[4].data = &box;

  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (vname != NULL) {
    gfs_derived_last->next = domain->variables;
    output->v = gfs_variable_from_name (gfs_derived_first, vname);
    if (output->v == NULL) {
      gts_file_variable_error (fp, var, "v", "unknown scalar `%s'", vname);
      g_free (vname);
      g_free (box);
      return;
    }
    g_free (vname);
  }
  
  if (box != NULL) {
    gchar * s = strtok (box, ",");

    output->box = GTS_BBOX (gts_object_new (GTS_OBJECT_CLASS (gts_bbox_class ())));
    if (s == NULL) {
      gts_file_variable_error (fp, var, "box", "expecting a number (x1)");
      g_free (box);
      return;
    }
    output->box->x1 = atof (s);
    s = strtok (NULL, ",");
    if (s == NULL) {
      gts_file_variable_error (fp, var, "box", "expecting a number (y1)");
      g_free (box);
      return;
    }
    output->box->y1 = atof (s);
    s = strtok (NULL, ",");
#if (!FTT_2D)
    if (s == NULL) {
      gts_file_variable_error (fp, var, "box", "expecting a number (z1)");
      g_free (box);
      return;
    }
    output->box->z1 = atof (s);
    s = strtok (NULL, ",");
#endif /* 3D */
    if (s == NULL) {
      gts_file_variable_error (fp, var, "box", "expecting a number (x2)");
      g_free (box);
      return;
    }
    output->box->x2 = atof (s);
    if (output->box->x2 < output->box->x1) {
      gts_file_variable_error (fp, var, "box", "x2 must be larger than x1");
      g_free (box);
      return;
    }
    s = strtok (NULL, ",");
    if (s == NULL) {
      gts_file_variable_error (fp, var, "box", "expecting a number (y2)");
      g_free (box);
      return;
    }
    output->box->y2 = atof (s);
    if (output->box->y2 < output->box->y1) {
      gts_file_variable_error (fp, var, "box", "y2 must be larger than y1");
      g_free (box);
      return;
    }
#if (!FTT_2D)
    s = strtok (NULL, ",");
    if (s == NULL) {
      gts_file_variable_error (fp, var, "box", "expecting a number (z2)");
      g_free (box);
      return;
    }
    output->box->z2 = atof (s);
    if (output->box->z2 < output->box->z1) {
      gts_file_variable_error (fp, var, "box", "z2 must be larger than z1");
      g_free (box);
      return;
    }
#endif /* 3D */
    g_free (box);
  }

  if (var[0].set || var[1].set)
    output->autoscale = FALSE;

  if (var[0].set && output->min > output->max) {
    gts_file_variable_error (fp, var, "min", 
	     "min `%g' must be smaller than or equal to max `%g'", 
			     output->min, output->max);
    return;
  }

  if (var[1].set && output->max < output->min) {
    gts_file_variable_error (fp, var, "max", 
	     "max `%g' must be larger than or equal to min `%g'", 
			     output->max, output->min);
    return;
  }
}

static void gfs_output_scalar_write (GtsObject * o, FILE * fp)
{
  GfsOutputScalar * output = GFS_OUTPUT_SCALAR (o);

  if (GTS_OBJECT_CLASS (gfs_output_scalar_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_output_scalar_class ())->parent_class->write) 
      (o, fp);

  fprintf (fp, " { v = %s", output->v->name);
  if (output->maxlevel >= 0)
    fprintf (fp, " maxlevel = %d", output->maxlevel);
  if (output->box != NULL)
#if FTT_2D
    fprintf (fp, " box = %g,%g,%g,%g", 
	     output->box->x1, output->box->y1, output->box->x2, output->box->y2);
#else  /* 3D */
    fprintf (fp, " box = %g,%g,%g,%g,%g,%g",
	     output->box->x1, output->box->y1, output->box->z1,
	     output->box->x2, output->box->y2, output->box->z2);
#endif /* 3D */
  if (!output->autoscale)
    fprintf (fp, " min = %g max = %g }", output->min, output->max);
  else
    fputs (" }", fp);
}

static gboolean gfs_output_scalar_event (GfsEvent * event,
					 GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_scalar_class ())->parent_class)->event)
      (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    GfsDomain * domain = GFS_DOMAIN (sim);
    
    if (output->v->derived) {
      gfs_variable_set_parent (output->v, domain);
      gfs_domain_cell_traverse (domain,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) output->v->derived, 
				output->v);
    }
    if (output->maxlevel >= 0)
        gfs_domain_cell_traverse (domain,
				  FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
				  (FttCellTraverseFunc) output->v->fine_coarse,
				  output->v);
    if (output->autoscale) {
      GtsRange stats = gfs_domain_stats_variable (domain, output->v, 
	     FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, output->maxlevel);

      output->min = stats.min;
      output->max = stats.max;
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_scalar_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_scalar_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_output_scalar_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_output_scalar_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_output_scalar_destroy;
}

static void gfs_output_scalar_init (GfsOutputScalar * object)
{
  object->v = gfs_p;
  object->min = object->max = 0.;
  object->autoscale = TRUE;
  object->maxlevel = -1;
  object->box = NULL;
}

GfsOutputClass * gfs_output_scalar_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_scalar_info = {
      "GfsOutputScalar",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_scalar_class_init,
      (GtsObjectInitFunc) gfs_output_scalar_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_scalar_info);
  }

  return klass;
}

/* GfsOutputScalarNorm: Object */

static gboolean gfs_output_scalar_norm_event (GfsEvent * event, 
					      GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_scalar_norm_class ())->parent_class)->event) (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    GfsNorm norm = gfs_domain_norm_variable (GFS_DOMAIN (sim), 
					     output->v,
					     FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
					     output->maxlevel);
    
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "%s time: %g first: % 10.3e second: % 10.3e infty: % 10.3e\n",
	     output->v->name,
	     sim->time.t,
	     norm.first, norm.second, norm.infty);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_scalar_norm_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_scalar_norm_event;
}

GfsOutputClass * gfs_output_scalar_norm_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_scalar_norm_info = {
      "GfsOutputScalarNorm",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_scalar_norm_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_scalar_norm_info);
  }

  return klass;
}

/* GfsOutputScalarStats: Object */

static gboolean gfs_output_scalar_stats_event (GfsEvent * event, 
					     GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_scalar_stats_class ())->parent_class)->event) (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    GtsRange stats = gfs_domain_stats_variable (GFS_DOMAIN (sim), 
						output->v,
						FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
						output->maxlevel);
    
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "%s time: %g min: %10.3e avg: %10.3e | %10.3e max: %10.3e\n",
	     output->v->name,
	     sim->time.t,
	     stats.min, stats.mean, stats.stddev, stats.max);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_scalar_stats_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_scalar_stats_event;
}

GfsOutputClass * gfs_output_scalar_stats_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_scalar_stats_info = {
      "GfsOutputScalarStats",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_scalar_stats_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_scalar_stats_info);
  }

  return klass;
}

/* GfsOutputScalarSum: Object */

static void add (FttCell * cell, gpointer * data)
{
  GfsSolidVector * solid = GFS_STATE (cell)->solid;
  gdouble vol = (solid ? solid->a : 1.)*ftt_cell_volume (cell);
  GfsVariable * v = data[0];
  gdouble * sum = data[1];

  *sum += vol*GFS_VARIABLE (cell, v->i);
}

static gboolean gfs_output_scalar_sum_event (GfsEvent * event, 
					     GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_scalar_sum_class ())->parent_class)->event) (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    gpointer data[2];
    gdouble sum = 0.;

    data[0] = output->v;
    data[1] = &sum;
    gfs_domain_cell_traverse (GFS_DOMAIN (sim),
			      FTT_PRE_ORDER, 
			      FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,
			      output->maxlevel,
			      (FttCellTraverseFunc) add, data);
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "%s time: %g sum: % 15.6e\n",
	     output->v->name, sim->time.t, sum);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_scalar_sum_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_scalar_sum_event;
}

GfsOutputClass * gfs_output_scalar_sum_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_scalar_sum_info = {
      "GfsOutputScalarSum",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_scalar_sum_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_scalar_sum_info);
  }

  return klass;
}

/* GfsOutputEnergy: Object */

static void add_energy (FttCell * cell, gpointer * data)
{
  GfsSolidVector * solid = GFS_STATE (cell)->solid;
  GfsStateVector * s = GFS_STATE (cell);
  GtsRange * ps = data[2];
  gdouble vol = (solid ? solid->a : 1.)*ftt_cell_volume (cell);
  gdouble * ke = data[0];
  gdouble * pe = data[1];

  *ke += vol*(s->u*s->u + s->v*s->v);;
  *pe += vol*(s->p - ps->mean)*(s->p - ps->mean);
}

static gboolean gfs_output_energy_event (GfsEvent * event,
					 GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_energy_class ())->parent_class)->event) 
      (event, sim)) {
    GfsOutput * output = GFS_OUTPUT (event);
    gpointer data[3];
    gdouble ke = 0., pe = 0.;
    GtsRange stats = gfs_domain_stats_variable (GFS_DOMAIN (sim), 
						gfs_p,
						FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
						-1);

    data[0] = &ke;
    data[1] = &pe;
    data[2] = &stats;
    gfs_domain_cell_traverse (GFS_DOMAIN (sim),
			      FTT_PRE_ORDER, 
			      FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,
			      -1,
			      (FttCellTraverseFunc) add_energy, data);
    fprintf (output->file->fp,
	     "Energy time: %g kinetic: %10.3e potential: %10.3e\n",
	     sim->time.t, ke, pe/sim->physical_params.g);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_energy_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_energy_event;
}

GfsOutputClass * gfs_output_energy_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_energy_info = {
      "GfsOutputEnergy",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_energy_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_energy_info);
  }

  return klass;
}

/* GfsOutputErrorNorm: Object */

static void output_error_norm_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_OUTPUT_ERROR_NORM (o)->s));

  (* GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class->destroy) 
    (o);
}

static void output_error_norm_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputErrorNorm * n;

  if (GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  n = GFS_OUTPUT_ERROR_NORM (*o);
  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return;
  }
  fp->scope_max++;
  gts_file_next_token (fp);
  while (fp->type != GTS_ERROR && fp->type != '}') {
    if (fp->type == '\n') {
      gts_file_next_token (fp);
      continue;
    }
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a parameter");
      return;
    }
    else if (!strcmp (fp->token->str, "unbiased")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);
      if (fp->type != GTS_INT) {
	gts_file_error (fp, "expecting an integer");
	return;
      }
      n->unbiased = atoi (fp->token->str);
      gts_file_next_token (fp);
    }
    else if (!strcmp (fp->token->str, "s")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);
      gfs_function_read (n->s, fp);
      if (fp->type == GTS_ERROR)
	return;
    }
    else if (!strcmp (fp->token->str, "v")) {
      GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);
      if (fp->type != GTS_STRING) {
	gts_file_error (fp, "expecting a variable name");
	return;
      }
      n->v = gfs_variable_from_name (domain->variables, fp->token->str);
      if (!n->v)
	n->v = gfs_domain_add_variable (domain, fp->token->str);
      g_assert (n->v);
      gts_file_next_token (fp);
    }
    else {
      gts_file_error (fp, "unknown identifier `%s'", fp->token->str);
      return;
    }
  }
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);
}

static void output_error_norm_write (GtsObject * o, FILE * fp)
{
  GfsOutputErrorNorm * n = GFS_OUTPUT_ERROR_NORM (o);

  if (GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class->write) 
      (o, fp);
  fputs (" { s = ", fp);
  gfs_function_write (n->s, fp);
  fprintf (fp, " unbiased = %d", n->unbiased);
  if (n->v != gfs_div)
    fprintf (fp, " v = %s }", n->v->name);
  else
    fputs (" }", fp);
}

static void compute_error (FttCell * cell, GfsOutputScalar * o)
{
  FttVector p;
  GfsSimulation * sim = gfs_object_simulation (o);

  if (o->v->centered)
    ftt_cell_pos (cell, &p);
  else
    gfs_cell_cm (cell, &p);
  GFS_VARIABLE (cell, GFS_OUTPUT_ERROR_NORM (o)->v->i) = GFS_VARIABLE (cell, o->v->i) -
    gfs_function_value (GFS_OUTPUT_ERROR_NORM (o)->s, &p, sim->time.t); 
}

static void remove_bias (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsNorm * norm = data[1];
  GFS_VARIABLE (cell, v->i) -= norm->bias;
}

static gboolean gfs_output_error_norm_event (GfsEvent * event, 
					     GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class)->event) (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    GfsVariable * v = GFS_OUTPUT_ERROR_NORM (event)->v;
    GfsNorm norm;

    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, 
			      FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,  
			      output->maxlevel,
			      (FttCellTraverseFunc) compute_error, output);
    norm = gfs_domain_norm_variable (GFS_DOMAIN (sim), v,
				     FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
				     output->maxlevel);
    if (GFS_OUTPUT_ERROR_NORM (event)->unbiased) {
      gpointer data[2];

      data[0] = v;
      data[1] = &norm;
      gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, 
				FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,  
				output->maxlevel,
				(FttCellTraverseFunc) remove_bias, data);
      norm = gfs_domain_norm_variable (GFS_DOMAIN (sim), v,
				       FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
				       output->maxlevel);
    }
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "%s time: %g first: % 10.3e second: % 10.3e infty: % 10.3e bias: %10.3e\n",
	     output->v->name, sim->time.t,
	     norm.first, norm.second, norm.infty, norm.bias);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_error_norm_class_init (GfsOutputClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = output_error_norm_destroy;
  GTS_OBJECT_CLASS (klass)->read = output_error_norm_read;
  GTS_OBJECT_CLASS (klass)->write = output_error_norm_write;
  GFS_EVENT_CLASS (klass)->event = gfs_output_error_norm_event;
}

static void output_error_norm_init (GfsOutputErrorNorm * e)
{
  e->s = gfs_function_new (gfs_function_class (), 0.);
  e->v = gfs_div;
}

GfsOutputClass * gfs_output_error_norm_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_error_norm_info = {
      "GfsOutputErrorNorm",
      sizeof (GfsOutputErrorNorm),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_error_norm_class_init,
      (GtsObjectInitFunc) output_error_norm_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_error_norm_info);
  }

  return klass;
}

/* GfsOutputCorrelation: Object */

static void compute_correlation (FttCell * cell, gpointer * data)
{
  GfsOutputScalar * o = data[0];
  gdouble * bias = data[1];
  gdouble * sum = data[2];
  gdouble * sumref = data[3];
  gdouble v, ref, w;
  FttVector p;
  GfsSimulation * sim = gfs_object_simulation (o);

  if (o->v->centered)
    ftt_cell_pos (cell, &p);
  else
    gfs_cell_cm (cell, &p);
  ref = gfs_function_value (GFS_OUTPUT_ERROR_NORM (o)->s, &p, sim->time.t);
  v = GFS_VARIABLE (cell, o->v->i) - *bias;
  w = ftt_cell_volume (cell)*(GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.);
  *sumref += ref*ref*w;
  *sum += v*ref*w;
}

static gboolean gfs_output_correlation_event (GfsEvent * event, 
					      GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_error_norm_class ())->parent_class)->event)
      (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    GfsVariable * v = GFS_OUTPUT_ERROR_NORM (event)->v;
    gdouble bias = 0., sum = 0., sumref = 0.;
    gpointer data[4];

    if (GFS_DOMAIN (sim)->pid != -1)
      g_assert_not_implemented ();

    if (GFS_OUTPUT_ERROR_NORM (event)->unbiased) {
      GfsNorm enorm;

      gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER,
				FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,
				output->maxlevel,
				(FttCellTraverseFunc) compute_error, output);
      enorm = gfs_domain_norm_variable (GFS_DOMAIN (sim), v,
					FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
					output->maxlevel);
      bias = enorm.bias;
    }
    data[0] = output;
    data[1] = &bias;
    data[2] = &sum;
    data[3] = &sumref;
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER,
			      FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,
			      output->maxlevel,
			      (FttCellTraverseFunc) compute_correlation, data);
    fprintf (GFS_OUTPUT (event)->file->fp,
	     "%s time: %g %10.3e\n",
	     output->v->name, sim->time.t, sumref > 0. ? sum/sumref : 0.);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_correlation_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_correlation_event;
}

GfsOutputClass * gfs_output_correlation_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_correlation_info = {
      "GfsOutputCorrelation",
      sizeof (GfsOutputErrorNorm),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_correlation_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_error_norm_class ()),
				  &gfs_output_correlation_info);
  }

  return klass;
}

/* GfsOutputSquares: Object */

static gboolean gfs_output_squares_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_squares_class ())->parent_class)->event) (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    
    gfs_write_squares (GFS_DOMAIN (sim), 
		       output->v, output->min, output->max,
		       FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL,
		       output->maxlevel, NULL, 
		       GFS_OUTPUT (event)->file->fp);
    fflush (GFS_OUTPUT (event)->file->fp);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_squares_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_squares_event;
}

GfsOutputClass * gfs_output_squares_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_squares_info = {
      "GfsOutputSquares",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_squares_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_squares_info);
  }

  return klass;
}

/* GfsOutputStreamline: Object */

static void gfs_output_streamline_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputStreamline * l = GFS_OUTPUT_STREAMLINE (*o);

  if (GTS_OBJECT_CLASS (gfs_output_streamline_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_output_streamline_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.x)");
    return;
  }
  l->p.x = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.y)");
    return;
  }
  l->p.y = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.z)");
    return;
  }
  l->p.z = atof (fp->token->str);
  gts_file_next_token (fp);
}

static void gfs_output_streamline_write (GtsObject * o, FILE * fp)
{
  GfsOutputStreamline * l = GFS_OUTPUT_STREAMLINE (o);

  if (GTS_OBJECT_CLASS (gfs_output_streamline_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_output_streamline_class ())->parent_class->write) 
      (o, fp);
  fprintf (fp, " %g %g %g", l->p.x, l->p.y, l->p.z);
}

static gboolean gfs_output_streamline_event (GfsEvent * event, 
					    GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_streamline_class ())->parent_class)->event) (event,sim)) {
    GSList * stream = gfs_streamline_new (GFS_DOMAIN (sim),
					  GFS_OUTPUT_STREAMLINE (event)->p,
					  GFS_OUTPUT_SCALAR (event)->v,
					  0., 0.,
					  TRUE);

    gfs_streamline_write (stream, GFS_OUTPUT (event)->file->fp);
    fflush (GFS_OUTPUT (event)->file->fp);
    gfs_streamline_destroy (stream);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_streamline_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_streamline_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_output_streamline_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_output_streamline_write;
}

GfsOutputClass * gfs_output_streamline_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_streamline_info = {
      "GfsOutputStreamline",
      sizeof (GfsOutputStreamline),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_streamline_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_streamline_info);
  }

  return klass;
}

/* GfsOutputStreakline: Object */

static void gfs_output_streakline_destroy (GtsObject * o)
{
  GfsOutputStreakline * l = GFS_OUTPUT_STREAKLINE (o);

  gts_object_destroy (GTS_OBJECT (l->p));
  g_slist_foreach (l->streak, (GFunc) gts_object_destroy, NULL);
  g_slist_free (l->streak);
  
  (* GTS_OBJECT_CLASS (gfs_output_streakline_class ())->parent_class->destroy) 
    (o);
}

static void gfs_output_streakline_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputStreakline * l = GFS_OUTPUT_STREAKLINE (*o);

  if (GTS_OBJECT_CLASS (gfs_output_streakline_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_output_streakline_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.x)");
    return;
  }
  l->p->x = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.y)");
    return;
  }
  l->p->y = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (p.z)");
    return;
  }
  l->p->z = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (ds)");
    return;
  }
  l->ds = atof (fp->token->str);
  gts_file_next_token (fp);

  l->streak = g_slist_prepend (l->streak, 
			       gts_point_new (gts_point_class (), 
					      l->p->x, l->p->y, l->p->z));
}

static void gfs_output_streakline_write (GtsObject * o, FILE * fp)
{
  GfsOutputStreakline * l = GFS_OUTPUT_STREAKLINE (o);

  if (GTS_OBJECT_CLASS (gfs_output_streakline_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_output_streakline_class ())->parent_class->write) 
      (o, fp);
  fprintf (fp, " %g %g %g %g", l->p->x, l->p->y, l->p->z, l->ds);
}

static void advect_point (GtsPoint * p, GfsSimulation * sim)
{
  gfs_domain_advect_point (GFS_DOMAIN (sim), p, sim->advection_params.dt);
}

static gboolean gfs_output_streakline_event (GfsEvent * event, 
					    GfsSimulation * sim)
{
  GfsOutputStreakline * l = GFS_OUTPUT_STREAKLINE (event);
  gboolean ret = FALSE;

  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_streakline_class ())->parent_class)->event) (event,sim)) {
    FILE * fp = GFS_OUTPUT (event)->file->fp;
    GSList * i = l->streak;

    fprintf (fp, "%s %u\n", 
	     GTS_OBJECT (event)->klass->info.name, 
	     g_slist_length (i));
    while (i) {
      fprintf (fp, "%g %g %g\n", 
	       GTS_POINT (i->data)->x, 
	       GTS_POINT (i->data)->y, 
	       GTS_POINT (i->data)->z);
      i = i->next;
    }
    l->started = ret = TRUE;
  }
  if (l->started) {
    if (gts_point_distance2 (l->p, l->streak->data) >= l->ds*l->ds)
      l->streak = g_slist_prepend (l->streak, 
				   gts_point_new (gts_point_class (), 
						  l->p->x, l->p->y, l->p->z));
    g_slist_foreach (l->streak, (GFunc) advect_point, sim);
  }
  return ret;
}

static void gfs_output_streakline_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_streakline_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_output_streakline_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_output_streakline_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_output_streakline_destroy;
}

static void gfs_output_streakline_init (GfsOutputStreakline * l)
{
  l->started = FALSE;
  l->ds = G_MAXDOUBLE;
  l->p = gts_point_new (gts_point_class (), 0., 0., 0.);
  l->streak = NULL;
}

GfsOutputClass * gfs_output_streakline_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_streakline_info = {
      "GfsOutputStreakline",
      sizeof (GfsOutputStreakline),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_streakline_class_init,
      (GtsObjectInitFunc) gfs_output_streakline_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS 
				  (gfs_output_scalar_class ()),
				  &gfs_output_streakline_info);
  }

  return klass;
}

/* GfsOutputPPM: Object */

static void gfs_output_ppm_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_output_ppm_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_output_ppm_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
#if (!FTT_2D)
  if (!GFS_IS_OCEAN (gfs_object_simulation (*o))) {
    gts_file_error (fp, 
		    "In more than two dimensions PPM output is possible\n"
		    "only for GfsOcean simulations");
    return;
  }
#endif /* 2D3 or 3D */
}

static gboolean gfs_output_ppm_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_ppm_class ())->parent_class)->event) 
      (event, sim)) {
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    GfsDomain * domain = GFS_IS_OCEAN (sim) ? GFS_OCEAN (sim)->toplayer : GFS_DOMAIN (sim);

    gfs_write_ppm (domain,
		   output->box,
		   output->v, output->min, output->max,
		   FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, output->maxlevel,
		   GFS_OUTPUT (event)->file->fp);
    fflush (GFS_OUTPUT (event)->file->fp);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_ppm_class_init (GfsOutputClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_output_ppm_read;
  GFS_EVENT_CLASS (klass)->event = gfs_output_ppm_event;
}

GfsOutputClass * gfs_output_ppm_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_ppm_info = {
      "GfsOutputPPM",
      sizeof (GfsOutputScalar),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_ppm_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_ppm_info);
  }

  return klass;
}
