/*
 * Copyright (C) 1998 Janne L.bÅˆf <jlof@mail.student.oulu.fi>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
/*
 * Modified on 10th June 2002, for porting this program
 * to the 'gtkglext-0.1.0` extension of gtk-2.0.
 *
 * Alif Wahid, <awah005@users.sourceforge.net>
 */
/*
 * Improved mouse operation.
 *
 * Naofumi Yasufuku  <naofumi@users.sourceforge.net>
 */
/*
 * Adapted for gfsview. May 2004.
 *
 * Stephane Popinet <popinet@users.sf.net>
 */

#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>

#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <gdk/gdkx.h>

#include <gtk/gtkgl.h>

#define SN_API_NOT_YET_FROZEN
#include <libsn/sn.h>

#if defined(__APPLE__)
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

#include "gfkgl.h"

#include "gl/trackball.h"
#include "glade/interface.h"
#include "glade/support.h"
#include "glade/callbacks.h"

static void pick (GtkWidget * widget, int winx, int winy, gboolean motion)
{
  GLdouble model[16], proj[16];
  GLint viewport[4];
  GLdouble x, y, z;
  GfsGlRay r;

  glGetDoublev (GL_MODELVIEW_MATRIX, model);
  glGetDoublev (GL_PROJECTION_MATRIX, proj);
  glGetIntegerv (GL_VIEWPORT, viewport);
  winy = widget->allocation.height - winy;
  g_return_if_fail (gluUnProject (winx, winy, 0., model, proj, viewport, &x, &y, &z));
  r.a.x = x; r.a.y = y; r.a.z = z;
  g_return_if_fail (gluUnProject (winx, winy, 1., model, proj, viewport, &x, &y, &z));
  r.b.x = x; r.b.y = y; r.b.z = z;
  
  gfk_gl_view_pick (lookup_widget (widget, "view"), &r, motion);
}

static GList * get_symmetries (GtkTreeModel * list)
{
  GtkTreeIter iter;  
  gboolean valid = gtk_tree_model_get_iter_first (list, &iter);
  GList * symmetry = NULL;

  while (valid) {
    gboolean visible;
    GfkGl * gl;

    gtk_tree_model_get (list, &iter, GL_COLUMN, &gl, VISIBLE_COLUMN, &visible, -1);
    if (visible && GFK_IS_GL_SYMMETRY (gl))
      symmetry = g_list_append (symmetry, gl->gl);
    valid = gtk_tree_model_iter_next (list, &iter);
  }
  return symmetry;
}

static gboolean
expose(GtkWidget      * widget,
       GdkEventExpose * event)
{
  GdkGLContext * glcontext = gtk_widget_get_gl_context (widget);
  GdkGLDrawable * gldrawable = gtk_widget_get_gl_drawable (widget);
  GfsGlViewParams * info = g_object_get_data (G_OBJECT (widget), "GfsGlViewParams");

  /* draw only last expose */
  if (event->count > 0)
    return TRUE;

  if (!gdk_gl_drawable_gl_begin(gldrawable, glcontext))
    return TRUE;

  /* basic initialization */
  if (info->do_init == TRUE) {
    gfs_gl_init_gl ();
    info->do_init = FALSE;
  }

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  GfsDomain * domain = g_object_get_data (G_OBJECT (widget), "sim");
  GtkWidget * view = lookup_widget (widget, "view");
  GtkTreeModel * list = gtk_tree_view_get_model (GTK_TREE_VIEW (lookup_widget (view, "gl_list")));
  GList * symmetries = get_symmetries (list);
  gdouble max = gfs_gl_domain_extent (domain, symmetries);
  g_list_free (symmetries);
  gluPerspective (info->fov, widget->allocation.width/(float) widget->allocation.height,
		  1., 1. + 2.*max);
  glMatrixMode (GL_MODELVIEW);

  /* draw object */
  glLoadIdentity ();
  glTranslatef (info->tx, info->ty, - (1. + max));
  gfs_gl_add_quats (info->dquat, info->quat, info->quat);
  GLfloat m[4][4];
  gfs_gl_build_rotmatrix (m, info->quat);
  glMultMatrixf (&m[0][0]);
  glScalef (info->sx, info->sy, info->sz);

  /* draw object */
  glClearColor (info->bg.r, info->bg.g, info->bg.b, 1);
  glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  gfk_gl_view_draw (view, GFSGL_SCREEN);

  /* swap backbuffer to front */
  if (gdk_gl_drawable_is_double_buffered (gldrawable))
    gdk_gl_drawable_swap_buffers (gldrawable);
  else
    glFlush ();

  gdk_gl_drawable_gl_end (gldrawable);

  return FALSE;
}

static gboolean
configure(GtkWidget         *widget,
          GdkEventConfigure *event)
{
  GdkGLContext *glcontext;
  GdkGLDrawable *gldrawable;

  g_return_val_if_fail(widget && event, FALSE);

  glcontext = gtk_widget_get_gl_context(widget);
  gldrawable = gtk_widget_get_gl_drawable(widget);

  /*** OpenGL BEGIN ***/
  if (!gdk_gl_drawable_gl_begin(gldrawable, glcontext))
    return TRUE;

  glViewport (0, 0, widget->allocation.width, widget->allocation.height);

  gdk_gl_drawable_gl_end(gldrawable);
  /*** OpenGL END ***/
  return TRUE;
}

static void
destroy(GtkWidget *widget)
{
  /* delete mesh info */
  GfsGlViewParams *info = (GfsGlViewParams*)g_object_get_data(G_OBJECT(widget), "GfsGlViewParams");

  g_free(info);
}

static gboolean
button_press(GtkWidget      *widget,
             GdkEventButton *event)
{
  if (event->button == 1 && event->state & GDK_CONTROL_MASK) {
    pick (widget, event->x, event->y, FALSE);
    return FALSE;
  }
  else {
    GfsGlViewParams * info = g_object_get_data (G_OBJECT (widget), "GfsGlViewParams");
    
    info->dquat[0] = 0.0;
    info->dquat[1] = 0.0;
    info->dquat[2] = 0.0;
    info->dquat[3] = 1.0;
    
    /* beginning of drag, reset mouse position */
    info->beginx = event->x;
    info->beginy = event->y;
    info->motion = TRUE;
    
    return FALSE;
  }
}

static gboolean
button_release(GtkWidget      *widget,
               GdkEventButton *event)
{
  GfsGlViewParams *info = (GfsGlViewParams*)g_object_get_data(G_OBJECT(widget), "GfsGlViewParams");

  info->dquat[0] = 0.0;
  info->dquat[1] = 0.0;
  info->dquat[2] = 0.0;
  info->dquat[3] = 1.0;

  info->dx = 0.0;
  info->dy = 0.0;
  info->motion = FALSE;
  if (info->res != info->base_res)
    gdk_window_invalidate_rect (widget->window, &widget->allocation, FALSE);

  return FALSE;
}

static gboolean
motion_notify(GtkWidget      *widget,
              GdkEventMotion *event)
{
  int x = 0;
  int y = 0;
  GdkModifierType state = 0;
  float width, height;
  gboolean redraw = FALSE;
  GfsGlViewParams *info = (GfsGlViewParams*)g_object_get_data(G_OBJECT(widget), "GfsGlViewParams");

  if (event->is_hint)
    gdk_window_get_pointer(event->window, &x, &y, &state);
  else {
    x = event->x;
    y = event->y;
    state = event->state;
  }

  width = widget->allocation.width;
  height = widget->allocation.height;

  if (state & GDK_CONTROL_MASK && state & GDK_BUTTON1_MASK) {
    if (x >= 0 && y >= 0 && x < width && y < height)
      pick (widget, x, y, TRUE);
    return TRUE;
  }

  if (state & GDK_BUTTON1_MASK) {
    /* drag in progress, simulate trackball */
    gfs_gl_trackball( info->dquat,
		      (2.0*info->beginx -            width) / width,
		      (          height - 2.0*info->beginy) / height,
		      (           2.0*x -            width) / width,
		      (          height -            2.0*y) / height );
#if FTT_2D
    if (!(event->state & GDK_SHIFT_MASK)) {
      /* constrain the rotations in the plane */
      info->dquat[0] = info->dquat[1] = 0.;
      gdouble n = sqrt (info->dquat[2]*info->dquat[2] + info->dquat[3]*info->dquat[3]);
      info->dquat[2] /= n; info->dquat[3] /= n;
    }
#endif
    info->dx = x - info->beginx;
    info->dy = y - info->beginy;
    
    /* orientation has changed, redraw mesh */
    redraw = TRUE;
  }

  if (state & GDK_BUTTON2_MASK) {
    /* zooming drag */
    info->fov += (0.1 + 3.*info->fov)*(y - info->beginy)/height;
    if (info->fov > 100.)
      info->fov = 100.;
    else if (info->fov < 0.01)
      info->fov = 0.01;
    /* zoom has changed, redraw mesh */
    redraw = TRUE;
  }

  if (state & GDK_BUTTON3_MASK) {
    /* translate drag */
    info->tx += (x - info->beginx)/width*0.02*(0.01 + 3.*info->fov);
    info->ty -= (y - info->beginy)/height*0.02*(0.01 + 3.*info->fov);
    
    /* translate has changed, redraw mesh */
    redraw = TRUE;
  }

  info->beginx = x;
  info->beginy = y;

  if (redraw)
    gdk_window_invalidate_rect (widget->window, &widget->allocation, FALSE);

  return TRUE;
}

static gboolean
scroll_notify(GtkWidget      *widget,
	      GdkEventScroll *event)
{
  GfsGlViewParams *info = (GfsGlViewParams*)g_object_get_data(G_OBJECT(widget), "GfsGlViewParams");

  switch (event->direction) {
  case GDK_SCROLL_UP:
    info->fov -= 0.02*(0.1 + 3.*info->fov);
    if (info->fov > 100.)
      info->fov = 100.;
    else if (info->fov < 0.01)
      info->fov = 0.01;
    /* zoom has changed, redraw mesh */
    gdk_window_invalidate_rect (widget->window, &widget->allocation, FALSE);
    break;
  case GDK_SCROLL_DOWN:
    info->fov += 0.02*(0.1 + 3.*info->fov);
    if (info->fov > 100.)
      info->fov = 100.;
    else if (info->fov < 0.01)
      info->fov = 0.01;
    /* zoom has changed, redraw mesh */
    gdk_window_invalidate_rect (widget->window, &widget->allocation, FALSE);
    break;
  default:
    ;
  }
  return FALSE;
}

static gboolean
key_press_event(GtkWidget   *widget,
                GdkEventKey *event)
{
  GfsGlViewParams *info = (GfsGlViewParams*)g_object_get_data(G_OBJECT(widget), "GfsGlViewParams");

  switch (event->keyval) {
  case GDK_plus:
    info->fov -= 0.02*(0.1 + 3.*info->fov);
    if (info->fov > 100.)
      info->fov = 100.;
    else if (info->fov < 0.01)
      info->fov = 0.01;
    /* zoom has changed, redraw mesh */
    gdk_window_invalidate_rect (widget->window, &widget->allocation, FALSE);
    break;
    
  case GDK_minus:
    info->fov += 0.02*(0.1 + 3.*info->fov);
    if (info->fov > 100.)
      info->fov = 100.;
    else if (info->fov < 0.01)
      info->fov = 0.01;
    /* zoom has changed, redraw mesh */
    gdk_window_invalidate_rect (widget->window, &widget->allocation, FALSE);
    break;
    
  default:
    return FALSE;
  }
  
  return TRUE;
}

static GtkWidget * gl_area_new (GdkGLConfig * glconfig)
{
  GtkWidget * glarea;
  
  /* create new OpenGL widget */
  glarea = gtk_drawing_area_new ();
  if (glarea == NULL) {
    fprintf (stderr, "gfsview: cannot create GtkDrawingArea widget\n");
    return NULL;
  }

  /* Set OpenGL-capability to the widget. */
  gtk_widget_set_gl_capability (GTK_WIDGET (glarea),
				glconfig,
				NULL,
				TRUE,
				GDK_GL_RGBA_TYPE);

  /* set up events and signals for OpenGL widget */
  gtk_widget_set_events (glarea,
			 GDK_EXPOSURE_MASK|
			 GDK_BUTTON_PRESS_MASK|
			 GDK_BUTTON_RELEASE_MASK|
			 GDK_POINTER_MOTION_MASK|
			 GDK_POINTER_MOTION_HINT_MASK);

  g_signal_connect (G_OBJECT (glarea), "expose_event",
		    G_CALLBACK (expose), NULL);
  g_signal_connect (G_OBJECT (glarea), "motion_notify_event",
		    G_CALLBACK (motion_notify), NULL);
  g_signal_connect (G_OBJECT (glarea), "button_press_event",
		    G_CALLBACK (button_press), NULL);
  g_signal_connect (G_OBJECT (glarea), "button_release_event",
		    G_CALLBACK (button_release), NULL);
  g_signal_connect (G_OBJECT (glarea), "scroll_event",
		    G_CALLBACK (scroll_notify), NULL);
  g_signal_connect (G_OBJECT (glarea), "configure_event",
		    G_CALLBACK (configure), NULL);
  g_signal_connect (G_OBJECT (glarea), "destroy",
		    G_CALLBACK (destroy), NULL);
  g_signal_connect (G_OBJECT (glarea), "key_press_event",
		    G_CALLBACK (key_press_event), glarea);

  return glarea;
}

typedef struct {
  GtkWidget * view;
  gboolean survive_broken_pipe;
} ScriptingArgs;

G_LOCK_DEFINE (main_loop_started);

static void send_scripting_message (GfkScriptingEvent event, GtkWidget * view, gpointer data)
{
  G_LOCK (scripting_pending);
  GfkScriptingMessage * msg = g_malloc (sizeof (GfkScriptingMessage));
  msg->event = event;
  msg->view = view;
  msg->data = data;
  if (!g_idle_add (gfk_receive_scripting_message, msg)) {
    g_warning ("could not send scripting message");
    g_free (msg->data);
    g_free (msg);
    G_UNLOCK (scripting_pending);
  }
}

static gpointer scripting (ScriptingArgs * s)
{
  gboolean scripting_off;
  fd_set rfds;

  G_LOCK (main_loop_started); /* make sure main loop has been started */
  G_UNLOCK (main_loop_started);

  FD_ZERO (&rfds);
  FD_SET (0, &rfds);
  while (select (1, &rfds, NULL, NULL, NULL) > 0) {
    GtsFile * fp = gts_file_new (stdin);

    while (fp->type != GTS_ERROR) {
      if (feof (stdin)) {
	if (!s->survive_broken_pipe) {
	  G_LOCK (scripting_pending); /* wait for pending scripting events */
	  gdk_threads_enter ();
	  gtk_main_quit ();
	  gdk_threads_leave ();
	}
	return NULL;
      }
      else if (fp->type == '\n')
	gts_file_next_token (fp);
      else if (fp->type == GTS_INT || !strcmp (fp->token->str, "GModule")) {
	GfsSimulation * sim = gfs_simulation_read (fp);

	if (sim == NULL)
	  break;

	G_LOCK (gfk_gl_scripting);
	scripting_off = !gfk_gl_scripting;
	G_UNLOCK (gfk_gl_scripting);
	if (scripting_off)
	  gts_object_destroy (GTS_OBJECT (sim));
	else {
	  gfs_simulation_init (sim);

	  gdk_threads_enter ();
	  gfk_gl_view_set_simulation (s->view, sim, "<stdin>");
	  gdk_threads_leave ();
	}
      }
      else if (fp->type == GTS_STRING) {
	if (!strcmp (fp->token->str, "Save") || !strcmp (fp->token->str, "Append")) {
	  gboolean append = !strcmp (fp->token->str, "Append");
	  GfsOutputFile * out = NULL;
	  GfsGl2PSParams * p;
	  gchar * fname;

	  gts_file_next_token (fp);
	  if (fp->type != GTS_STRING) {
	    gts_file_error (fp, "expecting a string (filename)");
	    break;
	  }
	  fname = g_strdup (fp->token->str);
	  gts_file_next_token (fp);
	  p = g_malloc (sizeof (GfsGl2PSParams));
	  gfs_gl2ps_params_read (p, fp);
	  if (fp->type == GTS_ERROR) {
	    g_free (fname);
	    g_free (p);
	    break;
	  }

	  G_LOCK (gfk_gl_scripting);
	  scripting_off = !gfk_gl_scripting;
	  G_UNLOCK (gfk_gl_scripting);

	  if (!scripting_off) {
	    if (append) {
	      if ((out = gfs_output_file_open (fp->token->str, "w")))
		p->fp = out->fp;
	    }
	    else /* Save */
	      p->fp = (!strcmp (fname, "stdout") ? stdout :
		       !strcmp (fname, "stderr") ? stderr :
		       fopen (fname, "w"));
	    if (p->fp == NULL)
	      fprintf (stderr, "gfsview: <stdin>: cannot open file `%s'\n", fname);
	    else {
	      send_scripting_message (append ? GFS_APPEND_EVENT : GFS_SAVE_EVENT, s->view, p);
	      /* p is freed by the receiver of the message */
	      p = NULL;
	      if (out) {
		/* Append mode, just free memory, do not close file */
		out->refcount++;
		gfs_output_file_close (out);
	      }
	    }
	  }
	  g_free (fname);
	  g_free (p);
	}
	else if (!strcmp (fp->token->str, "View")) {
	  G_LOCK (gfk_gl_scripting);
	  scripting_off = !gfk_gl_scripting;
	  G_UNLOCK (gfk_gl_scripting);

	  gdk_threads_enter ();
	  if (!gfk_gl_view_read_parameters (s->view, fp, scripting_off)) {
	    gdk_threads_leave ();
	    break;
	  }
	  gdk_threads_leave ();
	}
	else if (!strcmp (fp->token->str, "Clear")) {
	  gdk_threads_enter ();
	  gfk_gl_view_clear (s->view);
	  gdk_threads_leave ();
	  gts_file_next_token (fp);
	}
	else if (!strcmp (fp->token->str, "Echo")) {
	  gts_file_next_token (fp);
	  if (fp->type != '{') {
	    gts_file_error (fp, "expecting an opening brace");
	    break;
	  }
	  GString * echo = g_string_new ("");
	  guint scope = fp->scope_max;
	  gint c = gts_file_getc (fp);
	  while (c != EOF && fp->scope > scope) {
	    g_string_append_c (echo, c);
	    c = gts_file_getc (fp);
	  }
	  if (fp->scope != scope) {
	    g_string_free (echo, TRUE);
	    gts_file_error (fp, "parse error");
	    break;
	  }
	  gts_file_next_token (fp);

	  send_scripting_message (GFS_ECHO_EVENT, s->view, echo->str);
	  /* echo->str is freed by the receiver of the message */
	  g_string_free (echo, FALSE);
	}
	else {
	  gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
	  break;
	}
      }
      else {
	gts_file_error (fp, "expecting an integer got %d", fp->type);
	break;
      }
    }
    gdk_threads_enter ();
    gfk_gl_view_set_scripting (s->view, FALSE);
    GtkWidget * msg = gtk_message_dialog_new (GTK_WINDOW (s->view), 0,
					      GTK_MESSAGE_WARNING,
					      GTK_BUTTONS_CLOSE,
					      "GfsView: error in scripting thread\n"
					      "%d:%d: %s\n\n"
					      "Scripting is disabled",
					      fp->line, fp->pos, fp->error);
    gtk_dialog_run (GTK_DIALOG (msg));
    gtk_widget_destroy (msg);
    gdk_threads_leave ();
    gts_file_destroy (fp);
    while (fgetc (stdin) != EOF)
      ;
  }
  return NULL;
}

static void error_trap_push (SnDisplay * display,
			     Display   * xdisplay)
{
}

static void error_trap_pop (SnDisplay * display,
			    Display   * xdisplay)
{
  XSync (xdisplay, False); /* get all errors out of the queue */
}

static gboolean unlock_main_loop (gpointer data)
{
  G_UNLOCK (main_loop_started);  
  return FALSE;
}

int main (int argc, char * argv[])
{
  GdkGLConfig * glconfig;
  ScriptingArgs s;
  GtkWidget * glarea;
  int c = 0;

  /* initialize multithreading */
  g_thread_init (NULL);
  gdk_threads_init ();

  /* initialize gtk */
  gtk_set_locale ();
  gtk_init (&argc, &argv);

  add_pixmap_directory (PACKAGE_DATA_DIR "/pixmaps");

  /* initialize gtkglext */
  gtk_gl_init(&argc, &argv);

  /* initialize gfs */
  gfs_init (&argc, &argv);

  /* OpenGL drivers seem to often generate floating-point
     exceptions... turn them off so that people don't blame
     gfsview! */
  gfs_disable_floating_point_exceptions ();

  /* options */
  s.survive_broken_pipe = FALSE;
  while (c != EOF) {
    static struct option long_options[] = {
      {"survive-broken-pipe", no_argument, NULL, 's'},
      {"help", no_argument, NULL, 'h'},
      {"version", no_argument, NULL, 'V'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "hVs", long_options, &option_index))) {
    case 's': /* survive-broken-pipe */
      s.survive_broken_pipe = TRUE;
      break;
    case 'h': /* help */
      fprintf (stderr,
             "Usage: gfsview [OPTION] FILE1 FILE2 ...\n"
	     "The Gerris flow solver visualisation tool.\n"
	     "\n"
             "  -s    --survive-broken-pipe GfsView will not terminate\n"
	     "                              if the standard input pipe is broken\n"
	     "  -h    --help                display this help and exit\n"
	     "  -V    --version             output version information and exit\n"
	     "\n"
	     "Reports bugs to %s\n",
	     FTT_MAINTAINER);
      return 0; /* success */
      break;
    case 'V': /* version */
      fprintf (stderr,
	       "gfsview: %dD version %s\n",
	       FTT_DIMENSION, VERSION);
      return 0; /* succes */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `gfsview --help' for more information.\n");
      return 1; /* failure */
    }
  }

  /* Configure OpenGL-capable visual. */

  /* Try double-buffered visual */
  glconfig = gdk_gl_config_new_by_mode (GDK_GL_MODE_RGB |
					GDK_GL_MODE_DEPTH |
					GDK_GL_MODE_DOUBLE);
  if (glconfig == NULL) {
    g_print ("*** Cannot find the double-buffered visual.\n");
    g_print ("*** Trying single-buffered visual.\n");
    
    /* Try single-buffered visual */
    glconfig = gdk_gl_config_new_by_mode (GDK_GL_MODE_RGB |
					  GDK_GL_MODE_DEPTH);
    if (glconfig == NULL) {
      g_print ("*** No appropriate OpenGL-capable visual found.\n");
      exit(1);
    }
  }

  /* startup notification */
  Display * xdisplay = GDK_DISPLAY ();
  SnDisplay * display = sn_display_new (xdisplay, error_trap_push, error_trap_pop);
  SnLauncheeContext * 
    launched = sn_launchee_context_new_from_environment (display, DefaultScreen (xdisplay));

  /* Create view */
  glarea = gl_area_new (glconfig);
  gtk_widget_show (glarea);
  s.view = gfk_gl_view (glarea);
  
  /* Register scripting thread */
  G_LOCK (main_loop_started);
  if (!isatty (STDIN_FILENO) && g_thread_supported ()) {
    GError * error;
    if (g_thread_create ((GThreadFunc) scripting, &s, FALSE, &error))
      gfk_gl_view_set_scripting (s.view, TRUE);
    else {
      GtkWidget * msg = gtk_message_dialog_new (NULL, 0,
						GTK_MESSAGE_WARNING,
						GTK_BUTTONS_CLOSE,
						"GfsView could not start scripting thread:\n\n"
						"%s\n\n"
						"Scripting is disabled\n",
						error->message);
      gtk_dialog_run (GTK_DIALOG (msg));
      gtk_widget_destroy (msg);
      g_clear_error (&error);
    }
  }

  /* Read files on command line */
  gtk_widget_show (s.view);
  for (c = optind; c < argc; c++)
    gfk_gl_simulation_read (argv[c], s.view, TRUE);

  /* startup finished */
  if (launched)
    sn_launchee_context_complete (launched);

  g_timeout_add (0, unlock_main_loop, NULL);
  gtk_main ();

  return 0;
}
