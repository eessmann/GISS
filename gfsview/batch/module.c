/* Gerris - The GNU Flow Solver                       (-*-C-*-)
 * Copyright (C) 2010 National Institute of Water and Atmospheric Research
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

#include <errno.h>
#include <string.h>
#include <gerris/output.h>

#include "render.h"

/* GfsOutputView: Header */

typedef struct _GfsOutputView         GfsOutputView;

struct _GfsOutputView {
  /*< private >*/
  GfsOutput parent;

  /*< public >*/
  GfsGlViewParams view;
  GfsGl2PSParams p;
  GList * list;
  gchar * fname;
};

#define GFS_OUTPUT_VIEW(obj)            GTS_OBJECT_CAST (obj,\
					         GfsOutputView,\
					         gfs_output_view_class ())
#define GFS_IS_OUTPUT_VIEW(obj)         (gts_object_is_from_class (obj,\
						 gfs_output_view_class ()))

GfsEventClass * gfs_output_view_class  (void);

/* GfsOutputView: Object */

static void gfs_output_view_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_view_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type == '{') {
    gfs_gl2ps_params_read (&GFS_OUTPUT_VIEW (*o)->p, fp);
    if (fp->type == GTS_ERROR)
      return;
  }

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsView parameter file)");
    return;
  }
  FILE * f = fopen (fp->token->str, "r");
  if (f == NULL) {
    gts_file_error (fp, "could not open file '%s'\n%s", fp->token->str, strerror (errno));
    return;
  }
  GfsOutputView * output = GFS_OUTPUT_VIEW (*o);
  GtsFile * fp1 = gts_file_new (f);
  gfs_gl_view_params_read (&output->view, fp1);
  while (fp1->type == '\n')
    gts_file_next_token (fp1);
  while (fp1->type == GTS_STRING) {
    GfsGl * gl;
    if ((gl = gfs_gl_new_from_file (fp1))) {
      gl->p = &output->view;
      output->list = g_list_append (output->list, gl);
    }
    else if (fp1->type != GTS_ERROR) {
      gts_file_error (fp1, "unknown keyword `%s'", fp1->token->str);
      break;
    }
    while (fp1->type == '\n')
      gts_file_next_token (fp1);
  }
  if (fp1->type == GTS_ERROR) {
    gts_file_error (fp, "not a valid GfsView file\n%s:%d:%d: %s", 
		    fp->token->str, fp1->line, fp1->pos, fp1->error);
    gts_file_destroy (fp1);
    fclose (f);
  }
  else {
    gts_file_destroy (fp1);
    fclose (f);
    g_free (output->fname);
    output->fname = g_strdup (fp->token->str);
    gts_file_next_token (fp);
  }
}

static void gfs_output_view_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_view_class ())->parent_class->write) (o, fp);

  gfs_gl2ps_params_write (&GFS_OUTPUT_VIEW (o)->p, fp);
  fprintf (fp, " %s", GFS_OUTPUT_VIEW (o)->fname);
}

static void gfs_output_view_destroy (GtsObject * object)
{
  GfsOutputView * output = GFS_OUTPUT_VIEW (object);
  g_list_foreach (output->list, (GFunc) gts_object_destroy, NULL);
  g_list_free (output->list);
  g_free (output->fname);

  (* GTS_OBJECT_CLASS (gfs_output_view_class ())->parent_class->destroy) (object);
}

static gboolean gfs_output_view_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)) {
    GfsOutputView * output = GFS_OUTPUT_VIEW (event);
    g_list_foreach (output->list, (GFunc) gfs_gl_set_simulation, sim);
    gfs_gl_osmesa_render (&output->p, sim, &output->view, output->list, 
			  GFS_OUTPUT (event)->file->fp, 
			  !GFS_OUTPUT (event)->parallel);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_view_class_init (GtsObjectClass * klass)
{
  klass->read = gfs_output_view_read;
  klass->write = gfs_output_view_write;
  klass->destroy = gfs_output_view_destroy;
  GFS_EVENT_CLASS (klass)->event = gfs_output_view_event;
}

static void gfs_output_view_init (GfsOutputView * object)
{
  gfs_gl2ps_params_init (&object->p);
  gfs_gl_view_params_init (&object->view);
}

GfsEventClass * gfs_output_view_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_view_info = {
      "GfsOutputView",
      sizeof (GfsOutputView),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_output_view_class_init,
      (GtsObjectInitFunc) gfs_output_view_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_view_info);
  }

  return klass;
}

/* Initialize module */
const gchar   gfs_module_name[] = "gfsview";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_gl_init ();
  gfs_output_view_class ();
  return NULL;
}
