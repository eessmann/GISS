#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <math.h>
#include <glob.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gtk/gtk.h>

#include "callbacks.h"
#include "interface.h"
#include "support.h"
#include "../gfkgl.h"
#include "gl/trackball.h"

gpointer lookup_gl (gpointer widget);

gpointer lookup_gl (gpointer widget)
{
  GtkWidget * parent;
  GfkGl * gl;

  for (;;)
    {
      if (GTK_IS_MENU (widget))
        parent = gtk_menu_get_attach_widget (GTK_MENU (widget));
      else
        parent = GTK_WIDGET (widget)->parent;
      if (parent == NULL)
        break;
      widget = parent;
    }

  gl = g_object_get_data (G_OBJECT (widget), "GfkGl");
  if (!gl)
    g_warning ("GfkGl not found");
  return gl;
}

gpointer lookup_widget_params (gpointer widget, const gchar * widget_name);

gpointer lookup_widget_params (gpointer widget, const gchar * widget_name)
{
  GtkWidget *parent, *found_widget;

  for (;;)
    {
      if (GTK_IS_MENU (widget))
        parent = gtk_menu_get_attach_widget (GTK_MENU (widget));
      else
        parent = GTK_WIDGET (widget)->parent;
      if (parent == NULL || GTK_IS_NOTEBOOK (parent))
        break;
      widget = parent;
    }

  found_widget = (GtkWidget*) g_object_get_data (G_OBJECT (widget),
                                                 widget_name);
  return found_widget;
}

void
on_amin_toggled                        (GtkToggleButton *togglebutton,
					gpointer user_data)
{
  GfsGlScalar * gl = GFS_GL_SCALAR (GFK_GL (lookup_gl (togglebutton))->gl);
  GtkWidget * smin = lookup_widget_params (togglebutton, "spinbuttonmin");

  gl->amin = !gl->amin;
  gtk_widget_set_sensitive (smin, !gl->amin);
  if (gl->amin) {
    gtk_spin_button_set_value (GTK_SPIN_BUTTON (smin), gl->aminv);
    gl->min = gl->aminv;
  }
}

void
on_amax_toggled                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GfsGlScalar * gl = GFS_GL_SCALAR (GFK_GL (lookup_gl (togglebutton))->gl);
  GtkWidget * smax = lookup_widget_params (togglebutton, "spinbuttonmax");
  
  gl->amax = !gl->amax;
  gtk_widget_set_sensitive (smax, !gl->amax);
  if (gl->amax) {
    gtk_spin_button_set_value (GTK_SPIN_BUTTON (smax), gl->amaxv);
    gl->max = gl->amaxv;
  }
}

void
on_spinbuttonmin_changed               (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkSpinButton * min = GTK_SPIN_BUTTON (editable);
  GfkGl * gl = lookup_gl (editable);
  GfsGlScalar * gls = GFS_GL_SCALAR (gl->gl);
  gdouble v = gtk_spin_button_get_value_as_float (min);

  if (gls->min != v) {
    GtkSpinButton * max = lookup_widget_params (editable, "spinbuttonmax");
    GtkAdjustment * amax = gtk_spin_button_get_adjustment (max);
    gls->min = amax->lower = v;
    gtk_spin_button_set_adjustment (max, amax);
    gfk_gl_expose (gl);
  }
}

void
on_spinbuttonmax_changed               (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkSpinButton * max = GTK_SPIN_BUTTON (editable);
  GfkGl * gl = lookup_gl (editable);
  GfsGlScalar * gls = GFS_GL_SCALAR (gl->gl);
  gdouble v = gtk_spin_button_get_value_as_float (max);

  if (gls->max != v) {
    GtkSpinButton * min = lookup_widget_params (editable, "spinbuttonmin");
    GtkAdjustment * amin = gtk_spin_button_get_adjustment (min);
    gls->max = amin->upper = v;
    gtk_spin_button_set_adjustment (min, amin);
    gfk_gl_expose (gl);
  }
}

void
expose_gl                              (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  gfk_gl_expose (lookup_gl (togglebutton));
}

void
on_scale_changed                       (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (editable);
  GfsGlVectors * glv = GFS_GL_VECTORS (gl->gl);
  gdouble scale = gtk_spin_button_get_value_as_float (GTK_SPIN_BUTTON (editable));

  if (scale != glv->scale) {
    glv->scale = scale;
    gfk_gl_expose (gl);
  }
}

void
on_open1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gtk_widget_show (lookup_widget (GTK_WIDGET (menuitem), "filew"));
}


void
on_save1_activate                    (GtkMenuItem     *menuitem,
				      gpointer         user_data)
{
  gtk_widget_show (lookup_widget (GTK_WIDGET (menuitem), "gl2ps"));
}


void
on_delete_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget * list = lookup_widget (GTK_WIDGET (menuitem), "gl_list");
  GtkTreeSelection * select = gtk_tree_view_get_selection (GTK_TREE_VIEW (list));
  GtkTreeModel * model;
  GtkTreeIter iter;

  if (gtk_tree_selection_get_selected (select, &model, &iter)) {
    GfkGl * gl;

    gtk_tree_model_get (model, &iter, GL_COLUMN, &gl, -1);
    gtk_list_store_remove (GTK_LIST_STORE (model), &iter);
    if (g_object_get_data (G_OBJECT (list), "former") == gl)
      g_object_set_data (G_OBJECT (list), "former", NULL);
    gtk_widget_set_sensitive (lookup_widget (list, "properties1"), FALSE);
    gtk_widget_set_sensitive (lookup_widget (list, "delete3"), FALSE);
    gfk_gl_expose (gl);
    gts_object_destroy (GTS_OBJECT (gl));
  }
}


void
on_color_clicked                       (GtkButton       *button,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (button);
  gtk_widget_show (gl->color_selector);
}


void
on_maxlevel_changed                    (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (editable);
  guint maxlevel = gtk_spin_button_get_value_as_int (GTK_SPIN_BUTTON (editable));

  if (maxlevel != gl->gl->maxlevel) {
    gl->gl->maxlevel = maxlevel;
    gfk_gl_expose (gl);
  }
}


void
on_linewidth_value_changed             (GtkSpinButton    *spinbutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (spinbutton);
  gfloat lw = gtk_spin_button_get_value (spinbutton);

  if (lw != gl->gl->line_width) {
    gl->gl->line_width = lw;
    gfk_gl_expose (gl);
  }
}


void
on_finest_toggled                      (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (togglebutton);
  GtkWidget * widget = lookup_widget_params (togglebutton, "maxlevel");
    
  if (gl->gl->maxlevel == -1) {
    gl->gl->maxlevel = gtk_spin_button_get_value_as_int (GTK_SPIN_BUTTON (widget));
    gtk_widget_set_sensitive (widget, TRUE);
  }
  else {
    GtkAdjustment * amax = gtk_spin_button_get_adjustment (GTK_SPIN_BUTTON (widget));
    gtk_spin_button_set_value (GTK_SPIN_BUTTON (widget), amax->upper);
    gl->gl->maxlevel = -1;
    gtk_widget_set_sensitive (widget, FALSE);
  }
  gfk_gl_expose (gl);
}


void
on_constant1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_flat1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_smooth1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_csmooth1_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_color_selector_clicked              (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget * color_selection = lookup_widget (GTK_WIDGET (button), "colorselection1");
  GfkGl * gl = lookup_gl (button);
  GtkWidget * color_button = lookup_widget_params (gl->properties, "default_color");
  GdkColor c;

  gtk_widget_hide (lookup_widget (GTK_WIDGET (button), "color_selector"));
  gtk_color_selection_get_current_color (GTK_COLOR_SELECTION (color_selection), &c);
  gtk_widget_modify_bg (color_button, GTK_STATE_NORMAL, &c);
  gtk_widget_modify_bg (color_button, GTK_STATE_PRELIGHT, &c);
  gtk_widget_modify_bg (color_button, GTK_STATE_ACTIVE, &c);
  gl->gl->lc.r = c.red/(gdouble)65535;
  gl->gl->lc.g = c.green/(gdouble)65535;
  gl->gl->lc.b = c.blue/(gdouble)65535;
  gfk_gl_expose (gl);
}

void
on_clear1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gfk_gl_view_clear (GTK_WIDGET (menuitem));
}


void
on_properties1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget * list = lookup_widget (GTK_WIDGET (menuitem), "gl_list");
  GtkTreeSelection * select = gtk_tree_view_get_selection (GTK_TREE_VIEW (list));
  GtkTreeModel * model;
  GtkTreeIter iter;

  if (gtk_tree_selection_get_selected (select, &model, &iter)) {
    GfkGl * gl, * former = g_object_get_data (G_OBJECT (list), "former");

    gtk_tree_model_get (model, &iter, GL_COLUMN, &gl, -1);
    if (former != gl) {
      if (former) {
	gint x, y;
	gtk_window_get_position (GTK_WINDOW (former->params), &x, &y);
	gtk_window_move (GTK_WINDOW (gl->params), x, y);
	gtk_widget_show (gl->params);
	gtk_widget_hide (former->params);
      }
      else
	gtk_widget_show (gl->params);
      g_object_set_data (G_OBJECT (list), "former", gl);
    }
  }
}


void
on_preferences1_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gtk_widget_show (lookup_widget (GTK_WIDGET (menuitem), "preferences"));
}


void
on_resolution_changed                  (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkWidget * view = lookup_widget (GTK_WIDGET (editable), "view");
  GtkWidget * glarea = lookup_widget (GTK_WIDGET (view), "glarea");
  GfsGlViewParams * p = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");
  gdouble res = gtk_spin_button_get_value_as_float (GTK_SPIN_BUTTON (editable));

  if (res != p->base_res) {
    p->res = p->base_res = res;
    gdk_window_invalidate_rect (glarea->window, &glarea->allocation, FALSE);
  }
}


void
on_lines_closer_changed                (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkWidget * view = lookup_widget (GTK_WIDGET (editable), "view");
  GtkWidget * glarea = lookup_widget (GTK_WIDGET (view), "glarea");
  GfsGlViewParams * p = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");
  gdouble lc = gtk_spin_button_get_value_as_float (GTK_SPIN_BUTTON (editable));

  if (lc != p->lc) {
    p->lc = lc;
    gdk_window_invalidate_rect (glarea->window, &glarea->allocation, FALSE);
  }
}

void
on_reactivity_changed                  (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkWidget * view = lookup_widget (GTK_WIDGET (editable), "view");
  GtkWidget * glarea = lookup_widget (GTK_WIDGET (view), "glarea");
  GfsGlViewParams * p = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");
  
  p->reactivity = gtk_spin_button_get_value_as_float (GTK_SPIN_BUTTON (editable));
}

void
on_bgcolor_clicked                     (GtkButton       *button,
                                        gpointer         user_data)
{
  gtk_widget_show (g_object_get_data (G_OBJECT (button), "colorsel"));
}


void
on_bg_color_selector_clicked           (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget * color_selection = lookup_widget (GTK_WIDGET (button), "colorselection2");
  GtkWidget * color_button = lookup_widget (GTK_WIDGET (button), "bgcolor");
  GtkWidget * view = lookup_widget (GTK_WIDGET (color_button), "view");
  GtkWidget * glarea = lookup_widget (GTK_WIDGET (view), "glarea");
  GfsGlViewParams * p = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");
  GdkColor c;

  gtk_widget_hide (lookup_widget (GTK_WIDGET (button), "bg_color_selector"));
  gtk_color_selection_get_current_color (GTK_COLOR_SELECTION (color_selection), &c);
  gtk_widget_modify_bg (color_button, GTK_STATE_NORMAL, &c);
  gtk_widget_modify_bg (color_button, GTK_STATE_PRELIGHT, &c);
  gtk_widget_modify_bg (color_button, GTK_STATE_ACTIVE, &c);
  p->bg.r = c.red/(gdouble)65535;
  p->bg.g = c.green/(gdouble)65535;
  p->bg.b = c.blue/(gdouble)65535;
  gdk_window_invalidate_rect (glarea->window, &glarea->allocation, FALSE);
}

void gl2ps_ppm_set_sensitive (GtkWidget * w, gboolean s, gboolean s1);

void gl2ps_ppm_set_sensitive (GtkWidget * w, gboolean s, gboolean s1)
{
  gtk_widget_set_sensitive (lookup_widget (w, "orientationlabel"), s);
  gtk_widget_set_sensitive (lookup_widget (w, "orientationmenu"), s);
  gtk_widget_set_sensitive (lookup_widget (w, "label25"), s);
  gtk_widget_set_sensitive (lookup_widget (w, "spinbutton1"), s);

  if (s)
    gtk_widget_show (lookup_widget (w, "gl2ps_options"));
  else
    gtk_widget_hide (lookup_widget (w, "gl2ps_options"));
  if (s1)
    gtk_widget_show (lookup_widget (w, "ppm_options"));
  else
    gtk_widget_hide (lookup_widget (w, "ppm_options"));
}

#define gl2psparams(x) ((GfsGl2PSParams *) g_object_get_data (G_OBJECT (lookup_widget (GTK_WIDGET (x), "gl2ps")), "GfsGl2PSParams"))

static void set_extension (GtkMenuItem * item, const gchar * ext)
{
  GtkEntry * entry = GTK_ENTRY (lookup_widget (GTK_WIDGET (item), "entry2"));
  const gchar * fname = gtk_entry_get_text (entry);

  if (fname != NULL) {
    gchar * fname1 = g_strdup (fname), * c = fname1, * e = NULL;

    while (*c != '\0') {
      if (*c == '.')
	e = c;
      c++;
    }
    if (e != NULL) {
      gchar * fname2;

      *e = '\0';
      fname2 = g_strconcat (fname1, ".", ext, NULL);
      gtk_entry_set_text (entry, fname2);
      g_free (fname2);
    }
    g_free (fname1);
  }
}

void
on_portable_pixmap1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2ps_ppm_set_sensitive (GTK_WIDGET (menuitem), FALSE, TRUE);
  set_extension (menuitem, "ppm");
  gl2psparams (menuitem)->format = GFSGL_PPM;
}


void
on_postscript1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2ps_ppm_set_sensitive (GTK_WIDGET (menuitem), TRUE, FALSE);
  set_extension (menuitem, "ps");
  gl2psparams (menuitem)->format = GL2PS_PS;
}


void
on_encapsulated_postscript1_activate   (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2ps_ppm_set_sensitive (GTK_WIDGET (menuitem), TRUE, FALSE);
  set_extension (menuitem, "eps");
  gl2psparams (menuitem)->format = GL2PS_EPS;
}


void
on_portable_document_format1_activate  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2ps_ppm_set_sensitive (GTK_WIDGET (menuitem), TRUE, FALSE);
  set_extension (menuitem, "pdf");
  gl2psparams (menuitem)->format = GL2PS_PDF;
}


void
on_scalable_vector_graphics_activate (GtkMenuItem     *menuitem,
				      gpointer         user_data)
{
  gl2ps_ppm_set_sensitive (GTK_WIDGET (menuitem), TRUE, FALSE);
  set_extension (menuitem, "svg");
  gl2psparams (menuitem)->format = GL2PS_SVG;
}


void
on_gnuplot1_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2ps_ppm_set_sensitive (GTK_WIDGET (menuitem), TRUE, FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "orientationlabel"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "orientationmenu"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "label25"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "spinbutton1"), FALSE);  
  set_extension (menuitem, "gnu");
  gl2psparams (menuitem)->format = GFSGL_GNUPLOT;
}


void
on_wavefront_obj1_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2ps_ppm_set_sensitive (GTK_WIDGET (menuitem), TRUE, FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "orientationlabel"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "orientationmenu"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "label25"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "spinbutton1"), FALSE);  
  set_extension (menuitem, "obj");
  gl2psparams (menuitem)->format = GFSGL_OBJ;
}


void
on_kml_obj_activate             (GtkMenuItem     *menuitem,
				 gpointer         user_data)
{
  gl2ps_ppm_set_sensitive (GTK_WIDGET (menuitem), TRUE, FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "orientationlabel"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "orientationmenu"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "label25"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "spinbutton1"), FALSE);  
  set_extension (menuitem, "kml");
  gl2psparams (menuitem)->format = GFSGL_KML;
}


void
on_latex1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2ps_ppm_set_sensitive (GTK_WIDGET (menuitem), FALSE, FALSE);
  set_extension (menuitem, "tex");
  gl2psparams (menuitem)->format = GL2PS_TEX;
}


static void gl2pstoggle (gpointer w, GLint option)
{
  GfsGl2PSParams * q = gl2psparams (w);

  if (q->options & option)
    q->options &= ~option;
  else
    q->options |= option;
}


void
on_none1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2psparams (menuitem)->sort = GL2PS_NO_SORT;
}


void
on_barycenter1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2psparams (menuitem)->sort = GL2PS_SIMPLE_SORT;
}


void
on_bsp_tree1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2psparams (menuitem)->sort = GL2PS_BSP_SORT;
}


void
on_portrait1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2psparams (menuitem)->options &= ~GL2PS_LANDSCAPE;
}


void
on_landscape1_activate                 (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2psparams (menuitem)->options |= GL2PS_LANDSCAPE;
}


void
on_line_width_changed                  (GtkSpinButton    *editable,
                                        gpointer         user_data)
{
  gl2psparams (editable)->lw = gtk_spin_button_get_value_as_float (editable);
}


void
on_line_offset_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  gl2pstoggle (togglebutton, GL2PS_SIMPLE_LINE_OFFSET);
}


void
on_best_root_toggled                   (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  gl2pstoggle (togglebutton, GL2PS_BEST_ROOT);
}


void
on_no_text_toggled                     (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  gl2pstoggle (togglebutton, GL2PS_NO_TEXT);
}


void
on_no_pixmap_toggled                   (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  gl2pstoggle (togglebutton, GL2PS_NO_PIXMAP);
}


void
on_no_ps3_toggled                      (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  gl2pstoggle (togglebutton, GL2PS_NO_PS3_SHADING);
}

static gint browse_ok (GtkWidget * widget, GtkWidget * filew)
{
  GtkWidget * gl2ps = g_object_get_data (G_OBJECT (filew), "gl2ps");

  gtk_entry_set_text (GTK_ENTRY (lookup_widget (gl2ps, "entry2")),
		      gtk_file_selection_get_filename (GTK_FILE_SELECTION (filew)));
  gtk_widget_destroy (filew);
  return TRUE;
}

void
on_browse_clicked                      (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget * filew = gtk_file_selection_new ("Save as");

  gtk_file_selection_set_filename (GTK_FILE_SELECTION (filew),
				   gtk_entry_get_text (GTK_ENTRY (lookup_widget 
								  (GTK_WIDGET (button), 
								   "entry2"))));
  g_object_set_data (G_OBJECT (filew), "gl2ps", lookup_widget (GTK_WIDGET (button), "gl2ps"));
  g_signal_connect (G_OBJECT (GTK_FILE_SELECTION (filew)->ok_button), "clicked",
		    G_CALLBACK (browse_ok), filew);

  g_signal_connect_swapped (G_OBJECT (GTK_FILE_SELECTION (filew)->cancel_button), "clicked", 
			    G_CALLBACK (gtk_widget_destroy), filew);
  gtk_widget_show (filew);
}


static void refine_cell_corner (FttCell * cell, GfsDomain * domain)
{
  if (FTT_CELL_IS_LEAF (cell) && ftt_refine_corner (cell))
    ftt_cell_refine_single (cell, (FttCellInitFunc) gfs_cell_fine_init, domain);
}

void
on_save_as_button_clicked              (GtkButton       *button,
                                        gpointer         user_data)
{
  GfsGl2PSParams * p = gl2psparams (button);
  GtkWidget * view = lookup_widget (GTK_WIDGET (button), "view");
  const gchar * fname = 
    gtk_entry_get_text (GTK_ENTRY (lookup_widget (GTK_WIDGET (button), "entry2")));
  FILE * fp;

  if (fname == NULL || fname[0] == '\0') {
    GtkWidget * msg = gtk_message_dialog_new (GTK_WINDOW (view),
					      GTK_DIALOG_DESTROY_WITH_PARENT,
					      GTK_MESSAGE_WARNING,
					      GTK_BUTTONS_CLOSE,
					      "You must first select a filename");
    gtk_dialog_run (GTK_DIALOG (msg));
    gtk_widget_destroy (msg);
    return;
  }
  fp = fopen (fname, "w");
  if (!fp) {
    GtkWidget * msg = gtk_message_dialog_new (GTK_WINDOW (view),
					      GTK_DIALOG_DESTROY_WITH_PARENT,
					      GTK_MESSAGE_WARNING,
					      GTK_BUTTONS_CLOSE,
					      "Cannot open file `%s`",
					      fname);
    gtk_dialog_run (GTK_DIALOG (msg));
    gtk_widget_destroy (msg);
    return;
  }
  switch (p->format) {
  case GFSGL_GERRIS: {
    GtkWidget * glarea = lookup_widget (view, "glarea");
    GfsSimulation * sim = g_object_get_data (G_OBJECT (glarea), "sim");
    GfsDomain * domain = GFS_DOMAIN (sim);
    gint l;
    for (l = gfs_domain_depth (domain) - 2; l >= 0; l--)
      gfs_domain_cell_traverse (domain, 
				FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, l,
				(FttCellTraverseFunc) refine_cell_corner, domain);
    gfs_simulation_write (sim, -1, fp);
    /* fall through */
  }
  case GFSGL_GFSVIEW: {
    GtkWidget * glarea = lookup_widget (view, "glarea");
    GtkTreeModel * list = 
      gtk_tree_view_get_model (GTK_TREE_VIEW (lookup_widget (view, "gl_list")));
    GtkTreeIter iter;
    gboolean valid = gtk_tree_model_get_iter_first (list, &iter);

    fprintf (fp, "# GfsView %dD\n", FTT_DIMENSION);
    gfs_gl_view_params_write (g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams"), fp);
    fputc ('\n', fp);
    while (valid) {
      gboolean visible;
      GtsObject * gl;
      
      gtk_tree_model_get (list, &iter, GL_COLUMN, &gl, VISIBLE_COLUMN, &visible, -1);
      if (visible) {
	(* gl->klass->write) (gl, fp);
	fputc ('\n', fp);
      }
      valid = gtk_tree_model_iter_next (list, &iter);
    }
    break;
  }
  default:
    gfs_gl2ps (p, fp, fname, view);
  }
  fclose (fp);
  gtk_widget_hide (lookup_widget (view, "gl2ps"));
}

void
on_about1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gtk_widget_show (lookup_widget (GTK_WIDGET (menuitem), "about"));
}

static 
void 
rotate_view                               (GtkMenuItem     *menuitem, 
					   float * n, 
					   float alpha)
{
  GtkWidget * glarea = lookup_widget (GTK_WIDGET (menuitem), "glarea");
  GfsGlViewParams * p = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");

  p->dquat[0] = p->dquat[1] = p->dquat[2] = 0; p->dquat[3] = 1.;
  gfs_gl_axis_to_quat (n, alpha, p->quat);
  gdk_window_invalidate_rect (glarea->window, &glarea->allocation, FALSE);
}

void
on_front1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  float n[3] = { 0., 0., 1. };
  rotate_view (menuitem, n, 0.);
}


void
on_back1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  float n[3] = { 0., 1., 0. };
  rotate_view (menuitem, n, M_PI);
}


void
on_top1_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  float n[3] = { 1., 0., 0. };
  rotate_view (menuitem, n, - M_PI/2.);
}


void
on_bottom1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  float n[3] = { 1., 0., 0. };
  rotate_view (menuitem, n, M_PI/2.);
}


void
on_left1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  float n[3] = { 0., 1., 0. };
  rotate_view (menuitem, n, - M_PI/2.);
}


void
on_right1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  float n[3] = { 0., 1., 0. };
  rotate_view (menuitem, n, M_PI/2.);
}

void
on_entry2_changed                      (GtkEditable     *editable,
                                        gpointer         user_data)
{
  const gchar * fname = gtk_entry_get_text (GTK_ENTRY (editable));

  if (fname != NULL) {
    const gchar * c = fname, * e = NULL;

    while (*c != '\0') {
      if (*c == '.')
	e = c;
      c++;
    }
    if (e != NULL) {
      GtkOptionMenu * menu = GTK_OPTION_MENU (lookup_widget (GTK_WIDGET (editable), 
							     "format"));

      if (!strcmp (e, ".ppm")) {
	gl2ps_ppm_set_sensitive (GTK_WIDGET (menu), FALSE, TRUE);
	gtk_option_menu_set_history (menu, 0);
	gl2psparams (menu)->format = GFSGL_PPM;
      }
      else if (!strcmp (e, ".ps")) {
	gl2ps_ppm_set_sensitive (GTK_WIDGET (menu), TRUE, FALSE);
	gtk_option_menu_set_history (menu, 1);
	gl2psparams (menu)->format = GL2PS_PS;
      }
      else if (!strcmp (e, ".eps")) {
	gl2ps_ppm_set_sensitive (GTK_WIDGET (menu), TRUE, FALSE);
	gtk_option_menu_set_history (menu, 2);
	gl2psparams (menu)->format = GL2PS_EPS;
      }
      else if (!strcmp (e, ".pdf")) {
	gl2ps_ppm_set_sensitive (GTK_WIDGET (menu), TRUE, FALSE);
	gtk_option_menu_set_history (menu, 3);
	gl2psparams (menu)->format = GL2PS_PDF;
      }
      else if (!strcmp (e, ".svg")) {
	gl2ps_ppm_set_sensitive (GTK_WIDGET (menu), TRUE, FALSE);
	gtk_option_menu_set_history (menu, 4);
	gl2psparams (menu)->format = GL2PS_SVG;
      }
      else if (!strcmp (e, ".gnu")) {
	gl2ps_ppm_set_sensitive (GTK_WIDGET (menu), TRUE, FALSE);
	gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menu), "orientationlabel"), FALSE);
	gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menu), "orientationmenu"), FALSE);
	gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menu), "label25"), FALSE);
	gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menu), "spinbutton1"), FALSE);  
	gtk_option_menu_set_history (menu, 5);
	gl2psparams (menu)->format = GFSGL_GNUPLOT;
      }
      else if (!strcmp (e, ".obj")) {
	gl2ps_ppm_set_sensitive (GTK_WIDGET (menu), TRUE, FALSE);
	gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menu), "orientationlabel"), FALSE);
	gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menu), "orientationmenu"), FALSE);
	gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menu), "label25"), FALSE);
	gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menu), "spinbutton1"), FALSE);  
	gtk_option_menu_set_history (menu, 6);
	gl2psparams (menu)->format = GFSGL_OBJ;
      }
      else if (!strcmp (e, ".kml")) {
	gl2ps_ppm_set_sensitive (GTK_WIDGET (menu), TRUE, FALSE);
	gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menu), "orientationlabel"), FALSE);
	gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menu), "orientationmenu"), FALSE);
	gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menu), "label25"), FALSE);
	gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menu), "spinbutton1"), FALSE);  
	gtk_option_menu_set_history (menu, 7);
	gl2psparams (menu)->format = GFSGL_KML;
      }
      else if (!strcmp (e, ".tex")) {
	gl2ps_ppm_set_sensitive (GTK_WIDGET (menu), FALSE, FALSE);
	gtk_option_menu_set_history (menu, 8);
	gl2psparams (menu)->format = GL2PS_TEX;
      }
      else if (!strcmp (e, ".gfv")) {
	gl2ps_ppm_set_sensitive (GTK_WIDGET (menu), FALSE, FALSE);
	gtk_option_menu_set_history (menu, 9);
	gl2psparams (menu)->format = GFSGL_GFSVIEW;
      }
      else if (!strcmp (e, ".gfs")) {
	gl2ps_ppm_set_sensitive (GTK_WIDGET (menu), FALSE, FALSE);
	gtk_option_menu_set_history (menu, 10);
	gl2psparams (menu)->format = GFSGL_GERRIS;
      }
    }
  }
}

void
on_spinbuttonx_changed                 (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (editable);
  GfkGl2D * gl2 = GFK_GL2D (gl);
  gdouble nx = gtk_spin_button_get_value (GTK_SPIN_BUTTON (editable));

  if (nx != gl2->n.x) {
    gl2->n.x = nx;
    if (gts_vector_norm (&gl2->n.x) == 0.) {
      gl2->n.x = 1.;
      gtk_spin_button_set_value (GTK_SPIN_BUTTON (editable), 1.);
    }
    GFS_GL2D (gl->gl)->n = gl2->n;
    gfs_gl2D_update_plane (GFS_GL2D (gl->gl));
    gfk_gl2D_update_pos_bounds (gl2);
    gfk_gl_expose (gl);
  }
}


void
on_spinbuttony_changed                 (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (editable);
  GfkGl2D * gl2 = GFK_GL2D (gl);
  gdouble ny = gtk_spin_button_get_value (GTK_SPIN_BUTTON (editable));

  if (ny != gl2->n.y) {
    gl2->n.y = ny;
    if (gts_vector_norm (&gl2->n.x) == 0.) {
      gl2->n.y = 1.;
      gtk_spin_button_set_value (GTK_SPIN_BUTTON (editable), 1.);
    }
    GFS_GL2D (gl->gl)->n = gl2->n;
    gfs_gl2D_update_plane (GFS_GL2D (gl->gl));
    gfk_gl2D_update_pos_bounds (gl2);
    gfk_gl_expose (gl);    
  }
}


void
on_spinbuttonz_changed                 (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (editable);
  GfkGl2D * gl2 = GFK_GL2D (gl);
  gdouble nz = gtk_spin_button_get_value (GTK_SPIN_BUTTON (editable));

  if (nz != gl2->n.z) {
    gl2->n.z = nz;
    if (gts_vector_norm (&gl2->n.x) == 0.) {
      gl2->n.z = 1.;
      gtk_spin_button_set_value (GTK_SPIN_BUTTON (editable), 1.);
    }
    GFS_GL2D (gl->gl)->n = gl2->n;
    gfs_gl2D_update_plane (GFS_GL2D (gl->gl));
    gfk_gl2D_update_pos_bounds (gl2);
    gfk_gl_expose (gl);    
  }
}

void
on_spinbuttonpos_changed               (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (editable);
  GfsGl2D * gl2 = GFS_GL2D (gl->gl);
  gdouble pos = gtk_spin_button_get_value (GTK_SPIN_BUTTON (editable));

  if (pos != gl2->pos) {
    gl2->pos = pos;
    gfs_gl2D_update_plane (gl2);  
    gfk_gl_expose (gl);    
  }
}

void
on_spinbuttonlevel_changed             (GtkEditable     *editable,
                                        gpointer         user_data)
{
#if (!FTT_2D)
  GfkGl * gl = lookup_gl (editable);
  GfsGlIsosurface * gli = GFS_GL_ISOSURFACE (gl->gl);
  gdouble level = gtk_spin_button_get_value (GTK_SPIN_BUTTON (editable));

  if (level != gli->level) {
    gli->level = level;
    gfs_gl_isosurface_reset (gli);
    gfk_gl_expose (gl);
  }
#endif /* 3D */
}


void
on_reversed_toggled                    (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
#if (!FTT_2D)
  GfkGl * gl = lookup_gl (togglebutton);
  GfsGlIsosurface * gli = GFS_GL_ISOSURFACE (gl->gl);
  gli->reversed = !gli->reversed;
  gfk_gl_expose (gl);
#endif /* 3D */
}


void
on_solid_reversed_toggled              (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
#if (!FTT_2D)
  GfkGl * gl = lookup_gl (togglebutton);
  GfsGlSolid * gls = GFS_GL_SOLID (gl->gl);
  gls->reversed = !gls->reversed;
  gfk_gl_expose (gl);
#endif /* 3D */
}

gboolean
main_quit                              (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data)
{
  struct stat sb;
  if (fstat (STDIN_FILENO, &sb) != -1 && 
      (sb.st_mode & S_IFMT) == S_IFIFO &&
      !feof (stdin)) {
    GtkWidget * view = lookup_widget (widget, "view");
    GtkWidget * msg = gtk_message_dialog_new (GTK_WINDOW (view),
					      GTK_DIALOG_MODAL,
					      GTK_MESSAGE_INFO,
					      GTK_BUTTONS_OK_CANCEL,
					      "GfsView is connected to another application\n"
					      "Quitting GfsView now will also terminate this\n"
					      "application");
    if (gtk_dialog_run (GTK_DIALOG (msg)) == GTK_RESPONSE_OK) {
      gtk_widget_destroy (msg);
      gtk_main_quit ();
    }
    else
      gtk_widget_destroy (msg);
  }
  else
    gtk_main_quit ();
  return TRUE;
}

void
on_spinbuttoniso_value_changed         (GtkSpinButton   *spinbutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (spinbutton);
  GfsGlIsoline * gli = GFS_GL_ISOLINE (gl->gl);
  gdouble n = gtk_spin_button_get_value (spinbutton);

  if (n != gli->n) {
    gli->n = n;
    gfk_gl_expose (gl);
  }
}


void
on_entryiso_activate                   (GtkEntry        *entry,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (entry);
  GfsGlIsoline * gli = GFS_GL_ISOLINE (gl->gl);

  gfs_gl_isoline_set_levels (gli, gtk_entry_get_text (entry));
  gfk_gl_expose (gl);
}


void
on_ellipse_scale_changed               (GtkSpinButton   *spinbutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (spinbutton);
  GfsGlEllipses * gle = GFS_GL_ELLIPSES (gl->gl);
  gdouble scale = gtk_spin_button_get_value_as_float (spinbutton);

  if (scale != gle->scale) {
    gle->scale = scale;
    gfk_gl_expose (GFK_GL (gl));
  }
}


void
on_symbol_size_value_changed           (GtkSpinButton   *spinbutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (spinbutton);
  GfsGlLocation * gll = GFS_GL_LOCATION (gl->gl);
  gdouble size = gtk_spin_button_get_value_as_float (spinbutton);

  if (size != gll->size) {
    gll->size = size;
    gfk_gl_expose (gl);
  }
}


void
on_substract1_activate                 (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gtk_widget_show (lookup_widget (GTK_WIDGET (menuitem), "filews"));
}


void
on_back_clicked                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget * play = lookup_widget (GTK_WIDGET (button), "play");
  glob_t * files = g_object_get_data (G_OBJECT (play), "glob_t");

  if (files->gl_offs > 0) {
    GtkWidget * parent = g_object_get_data (G_OBJECT (play), "parent");

    gfk_gl_simulation_read (files->gl_pathv[--files->gl_offs], parent, TRUE);
    if (files->gl_offs == 0)
      gtk_widget_set_sensitive (GTK_WIDGET (button), FALSE);
    if (files->gl_offs == files->gl_pathc - 2)
      gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (button), "forward"), TRUE);
  } 
}


void
on_forward_clicked                     (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget * play = lookup_widget (GTK_WIDGET (button), "play");
  glob_t * files = g_object_get_data (G_OBJECT (play), "glob_t");

  if (files->gl_offs < files->gl_pathc - 1) {
    GtkWidget * parent = g_object_get_data (G_OBJECT (play), "parent");
   
    gfk_gl_simulation_read (files->gl_pathv[++files->gl_offs], parent, TRUE);
    if (files->gl_offs == files->gl_pathc - 1)
      gtk_widget_set_sensitive (GTK_WIDGET (button), FALSE);
    if (files->gl_offs == 1)
      gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (button), "back"), TRUE);
  }
}


void
on_gfsview_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2ps_ppm_set_sensitive (GTK_WIDGET (menuitem), FALSE, FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "orientationlabel"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "orientationmenu"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "label25"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "spinbutton1"), FALSE);
  set_extension (menuitem, "gfv");
  gl2psparams (menuitem)->format = GFSGL_GFSVIEW;
}


void
on_gerris_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gl2ps_ppm_set_sensitive (GTK_WIDGET (menuitem), FALSE, FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "orientationlabel"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "orientationmenu"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "label25"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (menuitem), "spinbutton1"), FALSE);
  set_extension (menuitem, "gfs");
  gl2psparams (menuitem)->format = GFSGL_GERRIS;
}


void
on_scripting_toggled                   (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  G_LOCK (gfk_gl_scripting);
  gfk_gl_scripting = !gfk_gl_scripting;
  G_UNLOCK (gfk_gl_scripting);
}

void
mark_edited                           (GtkEditable     *editable,
				       gpointer         user_data)
{
  g_object_set_data (G_OBJECT (editable), "edited", editable);
}

void
unmark_edited                          (GObject   *object,
                                        gpointer         user_data)
{
  g_object_set_data (object, "edited", NULL);
}

void
on_delete_streamline_clicked           (GtkButton       *button,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (button);

  if (gfs_gl_streamlines_remove_selected (GFS_GL_STREAMLINES (gl->gl))) {
    gtk_widget_set_sensitive (lookup_widget_params (GFK_GL_STREAMLINES (gl)->stream, "delete"), 
			      FALSE);
    gfk_gl_expose (gl);
  }
}

void
on_new_streamline_toggled              (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GfkGlStreamlines * gl = lookup_gl (togglebutton);

  if (gtk_toggle_button_get_active (togglebutton)) {
    if (GFS_GL_STREAMLINES (GFK_GL (gl)->gl)->selected) {
      GFS_GL_STREAMLINES (GFK_GL (gl)->gl)->selected = NULL;
      gtk_widget_set_sensitive (lookup_widget_params (gl->stream, "delete"), FALSE);
      gfk_gl_expose (GFK_GL (gl));
    }
    gl->edit = FALSE;
  }
  else
    gl->edit = TRUE;
}


void
on_show_cells_toggled                  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (togglebutton);
  
  GFS_GL_STREAMLINES (gl->gl)->show_cells = gtk_toggle_button_get_active (togglebutton);
  gfk_gl_expose (gl);
}

void
on_dmin_changed                        (GtkSpinButton   *spinbutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (spinbutton);
  GfsGlStreamlines * gls = GFS_GL_STREAMLINES (gl->gl);
  gdouble v = gtk_spin_button_get_value (spinbutton);

  if (gls->dmin != v) {
    GfkGlStreamlines * glk = GFK_GL_STREAMLINES (gl);

    if (v > 0.)
      gtk_widget_set_sensitive (lookup_widget_params (glk->stream, "evenly_spaced"), TRUE);
    else
      gtk_widget_set_sensitive (lookup_widget_params (glk->stream, "evenly_spaced"), FALSE);
    gls->dmin = v;
    gfs_gl_streamlines_reset (gls);
    gfk_gl_expose (gl);
  }
}

static gboolean callback (GfsGlStreamlines * gls, gpointer data)
{
  GfkGlStreamlines * gl = data;

  if (gtk_toggle_button_get_active (lookup_widget_params (gl->stream, "show_execution")))
    gfk_gl_expose (GFK_GL (gl));
  while (gtk_events_pending ())
    gtk_main_iteration ();
  return (g_object_get_data (lookup_widget_params (gl->stream, "stop_evenly"), "stop") != NULL);
}

void
on_evenly_spaced_clicked               (GtkButton       *button,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (button);
  GfkGlStreamlines * glk = GFK_GL_STREAMLINES (gl);
  GfsGlStreamlines * gls = GFS_GL_STREAMLINES (gl->gl);
  GtkWidget * stop = lookup_widget_params (glk->stream, "stop_evenly");
  GtkWidget * view = lookup_widget (gl->list, "view");
  GtkWidget * pref = g_object_get_data (G_OBJECT (view), "preferences");
  GtkToggleButton * scripting_on = GTK_TOGGLE_BUTTON (lookup_widget (pref, "scripting_on"));
  GtkToggleButton * scripting_off = GTK_TOGGLE_BUTTON (lookup_widget (pref, "scripting_off"));
  gboolean sensitive, scripting;

  gtk_widget_set_sensitive (stop, TRUE);
  g_object_set_data (G_OBJECT (stop), "stop", NULL);
  gtk_widget_set_sensitive (GTK_WIDGET (button), FALSE);

  if ((scripting = gtk_toggle_button_get_active (scripting_on)))
    gtk_toggle_button_set_active (scripting_off, TRUE);
  if ((sensitive = GTK_WIDGET_IS_SENSITIVE (lookup_widget (pref, "scripting"))))
    gtk_widget_set_sensitive (lookup_widget (pref, "scripting"), FALSE);

  gfs_gl_streamlines_evenly_spaced (gls, callback, gl);

  gtk_widget_set_sensitive (lookup_widget (pref, "scripting"), sensitive);
  gtk_toggle_button_set_active (scripting_on, scripting);

  gtk_widget_set_sensitive (stop, FALSE);
  gtk_widget_set_sensitive (GTK_WIDGET (button), TRUE);
  gfk_gl_expose (gl);
}


void
on_stop_evenly_clicked                 (GtkButton       *button,
                                        gpointer         user_data)
{
  g_object_set_data (G_OBJECT (button), "stop", button);
}


void
on_radius_changed                      (GtkSpinButton   *spinbutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (spinbutton);
  GfsGlStreamlines * gls = GFS_GL_STREAMLINES (gl->gl);
  gdouble v = gtk_spin_button_get_value (spinbutton);

  if (gls->radius != v) {
    gls->radius = v;
    gfs_gl_streamlines_update_display_lists (gls);
    gfk_gl_expose (gl);
  }
}

gboolean 
gl2ps_update_ppm_size (GtkWidget * widget, GdkEventExpose * event, GtkWidget * size)
{
  if (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON (lookup_widget (size, "screen")))) {
    guint width = widget->allocation.width;
    guint height = widget->allocation.height;
    gdouble * ratio = g_object_get_data (G_OBJECT (size), "ratio");

    if (width % 2) width++;
    if (height % 2) height++;
    gtk_spin_button_set_value (GTK_SPIN_BUTTON (lookup_widget (size, "width")), width);
    gtk_spin_button_set_value (GTK_SPIN_BUTTON (lookup_widget (size, "height")), height);
    *ratio = width/(gdouble) height;
  }
  return TRUE;
}


void
on_screen_ppm_toggled                  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget * size = lookup_widget (GTK_WIDGET (togglebutton), "size");

  if (gtk_toggle_button_get_active (togglebutton)) {
    GtkWidget * view = lookup_widget (GTK_WIDGET (togglebutton), "view");
    GtkWidget * glarea = g_object_get_data (G_OBJECT (view), "glarea");

    gl2ps_update_ppm_size (glarea, NULL, size);
    gtk_widget_set_sensitive (size, FALSE);
  }
  else
    gtk_widget_set_sensitive (size, TRUE);
}


void
on_width_value_changed                 (GtkSpinButton   *spinbutton,
                                        gpointer         user_data)
{
  GfsGl2PSParams * p = gl2psparams (spinbutton);

  p->width = gtk_spin_button_get_value (spinbutton);
  if (!gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON (lookup_widget (GTK_WIDGET (spinbutton), 
								      "screen"))) &&
      gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON (lookup_widget (GTK_WIDGET (spinbutton), 
								      "keep_ratio")))) {
    GtkWidget * size = lookup_widget (GTK_WIDGET (spinbutton), "size");
    GtkSpinButton * heights = GTK_SPIN_BUTTON (lookup_widget (size, "height"));
    gdouble * ratio = g_object_get_data (G_OBJECT (size), "ratio");
    guint height = p->width/(*ratio);

    if (height % 2) height++;
    if (fabs (height - gtk_spin_button_get_value (heights)) > 2.)
      gtk_spin_button_set_value (heights, height);
  }
}


void
on_height_value_changed                (GtkSpinButton   *spinbutton,
                                        gpointer         user_data)
{
  GfsGl2PSParams * p = gl2psparams (spinbutton);

  p->height = gtk_spin_button_get_value (spinbutton);
  if (!gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON (lookup_widget (GTK_WIDGET (spinbutton), 
								      "screen"))) &&
      gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON (lookup_widget (GTK_WIDGET (spinbutton), 
								      "keep_ratio")))) {
    GtkWidget * size = lookup_widget (GTK_WIDGET (spinbutton), "size");
    GtkSpinButton * widths = GTK_SPIN_BUTTON (lookup_widget (size, "width"));
    gdouble * ratio = g_object_get_data (G_OBJECT (size), "ratio");
    guint width = p->height*(*ratio);

    if (width % 2) width++;
    if (fabs (width - gtk_spin_button_get_value (widths)) > 2.)
      gtk_spin_button_set_value (widths, width);
  }
}


void
on_keep_ratio_toggled                  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active (togglebutton)) {
    GtkWidget * size = lookup_widget (GTK_WIDGET (togglebutton), "size");
    gdouble * ratio = g_object_get_data (G_OBJECT (size), "ratio");
    guint width = gtk_spin_button_get_value (GTK_SPIN_BUTTON (lookup_widget (size, "width")));
    guint height = gtk_spin_button_get_value (GTK_SPIN_BUTTON (lookup_widget (size, "height")));
    *ratio = width/(gdouble) height;
  }
}


void
on_information_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GfkGlCells * gl = lookup_gl (togglebutton);

  if (gtk_toggle_button_get_active (togglebutton)) {
    gtk_widget_set_sensitive (lookup_widget_params (gl->cells, "edit_level"), FALSE);
    gl->edit = FALSE;
  }
  else {
    gtk_widget_set_sensitive (lookup_widget_params (gl->cells, "edit_level"), TRUE);
    gl->edit = TRUE;
  }
}


void
on_locate_x_entry_changed              (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkEntry * entry = GTK_ENTRY (editable);
  GfkGl * gl = lookup_gl (editable);
  GfsGlLocate * locate = GFS_GL_LOCATE (gl->gl);
  gdouble v = atof (gtk_entry_get_text (entry));

  if (locate->p.x != v) {
    locate->p.x = v;
    gfk_gl_expose (gl);
  }
}


void
on_locate_y_entry_changed              (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkEntry * entry = GTK_ENTRY (editable);
  GfkGl * gl = lookup_gl (editable);
  GfsGlLocate * locate = GFS_GL_LOCATE (gl->gl);
  gdouble v = atof (gtk_entry_get_text (entry));

  if (locate->p.y != v) {
    locate->p.y = v;
    gfk_gl_expose (gl);
  }  
}


void
on_locate_z_entry_changed              (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkEntry * entry = GTK_ENTRY (editable);
  GfkGl * gl = lookup_gl (editable);
  GfsGlLocate * locate = GFS_GL_LOCATE (gl->gl);
  gdouble v = atof (gtk_entry_get_text (entry));

  if (locate->p.z != v) {
    locate->p.z = v;
    gfk_gl_expose (gl);
  }
}

void
on_x_scaling_changed                   (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkWidget * view = lookup_widget (GTK_WIDGET (editable), "view");
  GtkWidget * glarea = lookup_widget (GTK_WIDGET (view), "glarea");
  GfsGlViewParams * p = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");
  gdouble sx = gtk_spin_button_get_value_as_float (GTK_SPIN_BUTTON (editable));

  if (sx != p->sx) {
    p->sx = sx;
    gdk_window_invalidate_rect (glarea->window, &glarea->allocation, FALSE);
  }
}


void
on_y_scaling_changed                   (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkWidget * view = lookup_widget (GTK_WIDGET (editable), "view");
  GtkWidget * glarea = lookup_widget (GTK_WIDGET (view), "glarea");
  GfsGlViewParams * p = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");
  gdouble sy = gtk_spin_button_get_value_as_float (GTK_SPIN_BUTTON (editable));

  if (sy != p->sy) {
    p->sy = sy;
    gdk_window_invalidate_rect (glarea->window, &glarea->allocation, FALSE);
  }
}


void
on_z_scaling_changed                   (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkWidget * view = lookup_widget (GTK_WIDGET (editable), "view");
  GtkWidget * glarea = lookup_widget (GTK_WIDGET (view), "glarea");
  GfsGlViewParams * p = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");
  gdouble sz = gtk_spin_button_get_value_as_float (GTK_SPIN_BUTTON (editable));

  if (sz != p->sz) {
    p->sz = sz;
    gdk_window_invalidate_rect (glarea->window, &glarea->allocation, FALSE);
  }
}

void
on_vof_reversed_toggled                (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (togglebutton);
  GfsGlVOF * gli = GFS_GL_VOF (gl->gl);
  gli->reversed = !gli->reversed;
  gfk_gl_expose (gl);
}

void
on_vof_visible_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (togglebutton);
  GfsGlVOF * gli = GFS_GL_VOF (gl->gl);
  gli->draw_edges = !gli->draw_edges;
  gfk_gl_expose (gl);
}

void
on_linear_reversed_toggled             (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (togglebutton);
  GfsGlLinear * gli = GFS_GL_LINEAR (gl->gl);
  gli->reversed = !gli->reversed;
  gfk_gl_expose (gl);
}

void
on_location_label_toggled              (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (togglebutton);
  GfsGlLocation * gll = GFS_GL_LOCATION (gl->gl);
  gll->label = !gll->label;
  gfk_gl_expose (gl);
}


void
on_font_size_changed                   (GtkSpinButton   *spinbutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (spinbutton);
  gfs_gl_set_font_size (gl->gl, gtk_spin_button_get_value_as_float (spinbutton));
  gfk_gl_expose (gl);
}


void
on_vector_font_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (togglebutton);
  gfs_gl_set_raster_font (gl->gl, !gl->gl->use_raster_font);
  gfk_gl_expose (gl);
}


void
on_raster_font_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


static void set_labelpos_increment (GfkGl * gl)
{
  double maxlevel = gl->gl->maxlevel;
  if (maxlevel < 0) {
    GtkSpinButton * max = lookup_widget_params (gl->properties, "maxlevel");
    GtkAdjustment * amax = gtk_spin_button_get_adjustment (max);
    maxlevel = amax->upper;
  }
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++) {
    const char * name[] = { "labelx", "labely", "labelz" };
    GtkSpinButton * spos = lookup_widget_params (GFK_GL_LABEL (gl)->label, name[c]);
    GtkAdjustment * pos = gtk_spin_button_get_adjustment (spos);
    pos->page_increment = (pos->upper - pos->lower)/exp (maxlevel*log (2.));
    pos->step_increment = pos->page_increment/10.;
    gtk_spin_button_set_adjustment (spos, pos);
  }
}


void
on_labelx_changed                      (GtkSpinButton   *spinbutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (spinbutton);
  GfsGlLabel * l = GFS_GL_LABEL (gl->gl);
  gdouble x = gtk_spin_button_get_value (spinbutton);
  
  if (x != l->p.x) {
    l->p.x = x;
    set_labelpos_increment (gl);
    gfk_gl_expose (gl);
  }
}


void
on_labely_changed                      (GtkSpinButton   *spinbutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (spinbutton);
  GfsGlLabel * l = GFS_GL_LABEL (gl->gl);
  gdouble y = gtk_spin_button_get_value (spinbutton);
  
  if (y != l->p.y) {
    l->p.y = y;
    set_labelpos_increment (gl);
    gfk_gl_expose (gl);
  }
}


void
on_labelz_changed                      (GtkSpinButton   *spinbutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (spinbutton);
  GfsGlLabel * l = GFS_GL_LABEL (gl->gl);
  gdouble z = gtk_spin_button_get_value (spinbutton);
  
  if (z != l->p.z) {
    l->p.z = z;
    set_labelpos_increment (gl);
    gfk_gl_expose (gl);
  }
}


void
on_labelentry_activate                 (GtkEntry        *entry,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (entry);
  GfsGlLabel * l = GFS_GL_LABEL (gl->gl);

  gfs_gl_label_set_label (l, gtk_entry_get_text (entry), gl->gl->sim);
  
  gfk_gl_update_properties (gl);
  gfk_gl_expose (gl);
}


void
on_label_symbol_check_toggled          (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (togglebutton);
  GFS_GL_LABEL (gl->gl)->symbol = !GFS_GL_LABEL (gl->gl)->symbol;
  gfk_gl_expose (gl);
}


void
on_showscale_toggled                   (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GfkGl * gl = lookup_gl (togglebutton);
  GFS_GL_SCALAR (gl->gl)->show = !GFS_GL_SCALAR (gl->gl)->show;
  gfk_gl_expose (gl);
}
