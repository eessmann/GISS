/* Gerris - The GNU Flow Solver
 * Copyright (C) 2004-2012 National Institute of Water and Atmospheric
 * Research
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
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <glob.h>
#include <math.h>
#include <gts.h>
#include <gtk/gtkgl.h>
#include <gdk/x11/gdkglx.h>
#include <gdk/x11/gdkglglxext.h>
#if defined(__APPLE__)
#  include <OpenGL/glu.h>
#else
#  include <GL/glu.h>
#endif

#include "gfkgl.h"
#include "gl/trackball.h"
#include "glade/interface.h"
#include "glade/callbacks.h"
#include "glade/support.h"

G_LOCK_DEFINE (gfk_gl_scripting);
gboolean gfk_gl_scripting = FALSE;

gpointer lookup_gl               (gpointer widget);
gpointer lookup_widget_params    (gpointer widget, const gchar * widget_name);
void     gl2ps_ppm_set_sensitive (GtkWidget * w, gboolean s, gboolean s1);

static void    gl2D_destroy  (GtsObject * o);

typedef GtsFile * (* GfkFunctionSet) (GfkGl *, const gchar *);

static GtkWidget * gfk_function (GtkWidget * parent, 
				 GfkGl * gl, 
				 const gchar * name,
				 GfkFunctionSet set);

/* GfkGlSymmetry: Object */

static void gl2D_read (GtsObject ** o, GtsFile * fp)
{
  GtkWidget * params = GFK_GL2D (*o)->params;
  GfsGl2D * gl = GFS_GL2D (GFK_GL (*o)->gl);

  (* GTS_OBJECT_CLASS (gfk_gl2D_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GFK_GL2D (*o)->n = gl->n;
  gtk_spin_button_set_value (lookup_widget_params (params, "spinbuttonx"), gl->n.x);
  gtk_spin_button_set_value (lookup_widget_params (params, "spinbuttony"), gl->n.y);
  gtk_spin_button_set_value (lookup_widget_params (params, "spinbuttonz"), gl->n.z);
  gtk_spin_button_set_value (lookup_widget_params (params, "spinbuttonpos"), gl->pos);
}

static void pos_bounds (FttCell * root, GfsGl2D * gl, GtkAdjustment * pos)
{
  static FttVector d[8] = {{1.,1.,1.},{1.,1.,-1.},{1.,-1.,1.},{1.,-1.,-1.},
			   {-1.,1.,1.},{-1.,1.,-1.},{-1.,-1.,1.},{-1.,-1.,-1.}};
  FttVector p, o;
  gdouble h = ftt_cell_size (root)/2.;
  guint i;

  pos->lower = G_MAXDOUBLE;
  pos->upper = -G_MAXDOUBLE;
  ftt_cell_pos (root, &o);
  for (i = 0; i < 8; i++) {
    gdouble e;
    p.x = o.x + h*d[i].x; p.y = o.y + h*d[i].y; p.z = o.z + h*d[i].z;
    e = p.x*gl->n.x + p.y*gl->n.y + p.z*gl->n.z;
    if (e > pos->upper) pos->upper = e;
    if (e < pos->lower) pos->lower = e;
  }
}

static void find_bounds (FttCell * root, GfsGl2D * gl, GtkAdjustment * pos)
{
  GtkAdjustment p;

  pos_bounds (root, gl, &p);
  if (FTT_CELL_IS_LEAF (root)) {
    if (p.upper > pos->upper) pos->upper = p.upper;
    if (p.lower < pos->lower) pos->lower = p.lower;
  }
  else if (p.upper > pos->upper || p.lower < pos->lower) {
    FttCellChildren child;
    guint n;

    ftt_cell_children (root, &child);
    for (n = 0; n < FTT_CELLS; n++)
      if (child.c[n])
	find_bounds (child.c[n], gl, pos);
  }
}

static void find_box_bounds (GfsBox * box, gpointer * data)
{
  find_bounds (box->root, data[0], data[1]);
}

void gfk_gl2D_update_pos_bounds (GfkGl2D * gl)
{
  GtkSpinButton * spos = lookup_widget_params (gl->params, "spinbuttonpos");
  GtkAdjustment * pos = gtk_spin_button_get_adjustment (spos);

  if (GFK_IS_GL_PERIODIC (gl)) {
    pos->lower = -G_MAXDOUBLE;
    pos->upper = G_MAXDOUBLE;
  }
  else {
    gpointer data[2];
    data[0] = GFK_GL (gl)->gl;
    data[1] = pos;
    pos->lower = G_MAXDOUBLE;
    pos->upper = -G_MAXDOUBLE;
    gts_container_foreach (GTS_CONTAINER (GFK_GL (gl)->gl->sim), (GtsFunc) find_box_bounds, data);
  }
  gtk_spin_button_set_adjustment (spos, pos);
  if (pos->value < pos->lower)
    gtk_spin_button_set_value (spos, pos->lower);
  else if (pos->value > pos->upper)
    gtk_spin_button_set_value (spos, pos->upper);
}

static void set_pos_increment (GfkGl * gl)
{
  GtkSpinButton * maxlevel = lookup_widget_params (gl->properties, "maxlevel");
  GtkAdjustment * amax = gtk_spin_button_get_adjustment (maxlevel);
  GtkSpinButton * spos = lookup_widget_params (GFK_GL2D (gl)->params, "spinbuttonpos");
  GtkAdjustment * pos = gtk_spin_button_get_adjustment (spos);
  pos->step_increment = pos->page_increment = 1./exp (amax->upper*log (2.));
  gtk_spin_button_set_adjustment (spos, pos);
}

static void gl2D_set_simulation (GfkGl * gl, GfsSimulation * sim)
{
  (* GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl2D_class ())->parent_class)->set_simulation) (gl, sim);
  gfk_gl2D_update_pos_bounds (GFK_GL2D (gl));
  set_pos_increment (gl);
}

static void gl2D_post_init (GfkGl * gl)
{
  GfkGl2D * gl2 = GFK_GL2D (gl);

  (* GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl2D_class ())->parent_class)->post_init) (gl);

  FttComponent c;
  for (c = 0; c < 3; c++) {
    static gchar * name[] = {"spinbuttonx", "spinbuttony", "spinbuttonz"};
    GtkSpinButton * spin = lookup_widget_params (gl2->params, name[c]);
    gdouble v = gtk_spin_button_get_value (spin);
    if (v != (&gl2->n.x)[c])
      gtk_spin_button_set_value (spin, (&gl2->n.x)[c]);
  }

  GFS_GL2D (gl->gl)->n = gl2->n;
  gfs_gl2D_update_plane (GFS_GL2D (gl->gl));

  if (gl->gl->sim) {
    gfk_gl2D_update_pos_bounds (gl2);
    set_pos_increment (gl);
  }
}

static gchar * gl_symmetry_name (GfkGlClass * klass)
{
  static gchar name[] = "Symmetry";
  return name;
}

static GtkWidget * gl_symmetry_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "symmetry-16x16.png");
}

static void gl_symmetry_class_init (GfkGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gl2D_destroy;
  GTS_OBJECT_CLASS (klass)->read = gl2D_read;
  klass->post_init = gl2D_post_init;
  klass->set_simulation = gl2D_set_simulation;

  klass->gl_class = gfs_gl_symmetry_class ();
  klass->name = gl_symmetry_name;
  klass->icon = gl_symmetry_icon;
}

static void gl_symmetry_init (GfkGl2D * object)
{
  object->params = create_gl2D_params ();
  object->n.x = 1.; object->n.y = 0.; object->n.z = 0.;
  gfk_gl_prepend_params (GFK_GL (object), object->params, gtk_label_new ("2D Plane"));
#if FTT_2D
  gtk_widget_hide (lookup_widget_params (object->params, "spinbuttonz"));
#endif  
  gtk_widget_show (object->params);
}

GfkGlClass * gfk_gl_symmetry_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_symmetry_info = {
      "GfkGlSymmetry",
      sizeof (GfkGl2D),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_symmetry_class_init,
      (GtsObjectInitFunc) gl_symmetry_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_class ()), &gfk_gl_symmetry_info);
  }

  return klass;
}

/* GfkGlPeriodic: Object */

static gchar * gl_periodic_name (GfkGlClass * klass)
{
  static gchar name[] = "Periodic";
  return name;
}

static GtkWidget * gl_periodic_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "periodic-16x16.png");
}

static void gl_periodic_class_init (GfkGlClass * klass)
{
  klass->gl_class = gfs_gl_periodic_class ();
  klass->name = gl_periodic_name;
  klass->icon = gl_periodic_icon;
}

GfkGlClass * gfk_gl_periodic_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_periodic_info = {
      "GfkGlPeriodic",
      sizeof (GfkGl2D),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_periodic_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_symmetry_class ()), 
				  &gfk_gl_periodic_info);
  }

  return klass;
}

#if FTT_2D
# include "gfkgl2D.h"
#else  /* 3D */
# include "gfkgl3D.h"
#endif /* 3D */

/* GfkGl: Object */

static void gl_read (GtsObject ** o, GtsFile * fp)
{
  GtkWidget * properties =  GFK_GL (*o)->properties;
  GfsGl * gl = GFK_GL (*o)->gl;
  GtsObject * object = GTS_OBJECT (gl);
  GtkSpinButton * maxlevel = lookup_widget_params (properties, "maxlevel");
  GtkSpinButton * linewidth = lookup_widget_params (properties, "linewidth");
  GtkAdjustment * amax = gtk_spin_button_get_adjustment (maxlevel);

  (* object->klass->read) (&object, fp);
  if (fp->type == GTS_ERROR)
    return;

  gfk_gl_set_color (GFK_GL (*o), gl->lc);

  gtk_spin_button_set_value (linewidth, gl->line_width);

  if (gl->maxlevel == -1)
    gtk_spin_button_set_value (maxlevel, amax->upper);
  else {
    gtk_spin_button_set_value (maxlevel, gl->maxlevel);
    gl->maxlevel = -1;
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (lookup_widget_params (properties, "finest")), 
				  FALSE);
  }

  GtkWidget * font = GFK_GL (*o)->font;
  gtk_spin_button_set_value (lookup_widget_params (font, "font_size"), gl->font_size);
  if (!gl->use_raster_font) {
    gl->use_raster_font = TRUE;
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (lookup_widget_params (font, "vector_font")), 
				  TRUE);
  }
}

static void gl_write (GtsObject * o, FILE * fp)
{
  GfsGl * gl = GFK_GL (o)->gl;

  (* GTS_OBJECT (gl)->klass->write) (GTS_OBJECT (gl), fp);
}

static void gl_destroy (GtsObject * object)
{
  GfkGl * gl = GFK_GL (object);

  gtk_widget_destroy (gl->params);
  gtk_widget_destroy (gl->color_selector);
  gts_object_destroy (GTS_OBJECT (gl->gl));
  g_free (gl->props);
  (* GTS_OBJECT_CLASS (gfk_gl_class ())->parent_class->destroy) (object);
}

static gchar * gfk_gl_get_name (GfkGl * gl)
{
  return (* GFK_GL_CLASS (GTS_OBJECT (gl)->klass)->name) (GFK_GL_CLASS (GTS_OBJECT (gl)->klass));
}

static GtkWidget * gfk_gl_get_icon (GfkGl * gl)
{
  return (* GFK_GL_CLASS (GTS_OBJECT (gl)->klass)->icon) (GFK_GL_CLASS (GTS_OBJECT (gl)->klass));
}

static gchar * gfk_gl_get_properties (GfkGl * gl)
{
  if (GFK_GL_CLASS (GTS_OBJECT (gl)->klass)->properties)
    return (* GFK_GL_CLASS (GTS_OBJECT (gl)->klass)->properties) (gl);
  return gfk_gl_get_name (gl);
}

static void gl_update_interface (GfkGl * gl)
{
  GtkSpinButton * maxlevel = lookup_widget_params (gl->properties, "maxlevel");
  GtkAdjustment * amax = gtk_spin_button_get_adjustment (maxlevel);
  
  amax->upper = gfs_domain_depth (GFS_DOMAIN (gl->gl->sim));
  gtk_spin_button_set_adjustment (maxlevel, amax);
  if (gl->gl->maxlevel == -1) {
    gtk_spin_button_set_value (maxlevel, amax->upper);
    gl->gl->maxlevel = -1;
  }
  else {
    gtk_spin_button_set_value (maxlevel, gl->gl->maxlevel);
    GtkToggleButton * finest = GTK_TOGGLE_BUTTON (lookup_widget_params (gl->properties, "finest"));
    if (gtk_toggle_button_get_active (finest)) {
      gl->gl->maxlevel = -1;
      gtk_toggle_button_set_active (finest, FALSE);
    }
  }
  gtk_spin_button_set_value (lookup_widget_params (gl->font, "font_size"), gl->gl->font_size);
}

static void gl_set_simulation (GfkGl * gl, GfsSimulation * sim)
{
  gfs_gl_set_simulation (gl->gl, sim);
  gl_update_interface (gl);
}

static void gl_post_init (GfkGl * gl)
{
  gtk_window_set_transient_for (GTK_WINDOW (gl->params),
				GTK_WINDOW (gtk_widget_get_toplevel (gl->glarea)));
  gtk_window_set_destroy_with_parent (GTK_WINDOW (gl->params), FALSE);
  gtk_window_set_position (GTK_WINDOW (gl->params), GTK_WIN_POS_CENTER_ON_PARENT);
  gtk_window_set_title (GTK_WINDOW (gl->params), gfk_gl_get_name (gl));

  gtk_window_set_transient_for (GTK_WINDOW (gl->color_selector),
				GTK_WINDOW (gtk_widget_get_toplevel (gl->glarea)));
  gtk_window_set_destroy_with_parent (GTK_WINDOW (gl->color_selector), FALSE);
  gtk_window_set_position (GTK_WINDOW (gl->color_selector), GTK_WIN_POS_CENTER_ON_PARENT);

  if (gl->gl->sim)
    gl_update_interface (gl);
}

static void gl_add_gl (GtkWidget * glarea, GtkWidget * list, GfkGl * gl)
{
  GtkListStore * store;
  GtkTreeSelection * select;
  GtkTreeIter iter;
  GtkImage * icon;
  
  g_return_if_fail (glarea != NULL);
  g_return_if_fail (list != NULL);
  g_return_if_fail (gl != NULL);

  store = GTK_LIST_STORE (gtk_tree_view_get_model (GTK_TREE_VIEW (list)));
  select = gtk_tree_view_get_selection (GTK_TREE_VIEW (list));
  icon = GTK_IMAGE (gfk_gl_get_icon (gl));
  
  gtk_list_store_append (store, &iter);
  gtk_list_store_set (store, &iter,
		      VISIBLE_COLUMN,    TRUE,
		      ICON_COLUMN,       icon ? gtk_image_get_pixbuf (icon) : NULL,
		      PROPERTIES_COLUMN, gfk_gl_get_properties (gl),
		      GL_COLUMN, gl,
		      SELECTED_COLUMN,   TRUE,
		      -1);
  gtk_tree_selection_select_iter (select, &iter);

  gtk_widget_set_sensitive (lookup_widget (glarea, "save1"), TRUE);
  gtk_widget_set_sensitive (lookup_widget (glarea, "edit1"), TRUE);
  gtk_widget_set_sensitive (lookup_widget (glarea, "view1"), TRUE);
}

static void gl_add (GObject * tool)
{
  GtkWidget * glarea = g_object_get_data (tool, "glarea");
  GtkWidget * list = g_object_get_data (tool, "list");  
  GfkGlClass * klass = g_object_get_data (tool, "klass");
  GfsSimulation * sim = g_object_get_data (G_OBJECT (glarea), "sim");
  GfkGl * gl = gfk_gl_new (klass, glarea, list);

  gfk_gl_set_simulation (gl, sim);
  gl_add_gl (glarea, list, gl);
  gfk_gl_expose (gl);
}

static GtkWidget * gl_icon (GfkGlClass * klass)
{
  return NULL;
}

static gchar * gl_name (GfkGlClass * klass)
{
  return GTS_OBJECT_CLASS (klass)->info.name;
}

static void gl_class_init (GfkGlClass * klass)
{
  klass->gl_class = gfs_gl_class ();
  klass->post_init = gl_post_init;
  klass->set_simulation = gl_set_simulation;
  klass->icon = gl_icon;
  klass->name = gl_name;
  GTS_OBJECT_CLASS (klass)->read = gl_read;
  GTS_OBJECT_CLASS (klass)->write = gl_write;
  GTS_OBJECT_CLASS (klass)->destroy = gl_destroy;
}

static gboolean hide_params (GtkWidget * p)
{
  GfkGl * gl = g_object_get_data (G_OBJECT (p), "GfkGl");

  gtk_widget_hide (p);
  g_object_set_data (G_OBJECT (gl->list), "former", NULL);
  return TRUE;
}

void gfk_gl_prepend_params (GfkGl * gl, GtkWidget * widget, GtkWidget * label)
{
  gtk_notebook_prepend_page (g_object_get_data (G_OBJECT (gl->params), "book"), widget, label);
}

void gfk_gl_set_color (GfkGl * gl, GtsColor c)
{
  GdkColor gc;
  GtkWidget * b;

  g_return_if_fail (gl != NULL);

  gl->gl->lc = c;
  b = lookup_widget_params (gl->properties, "default_color");
  gc.red = c.r*65535.;
  gc.green = c.g*65535.;
  gc.blue = c.b*65535.;
  gtk_widget_modify_bg (b, GTK_STATE_NORMAL, &gc);
  gtk_widget_modify_bg (b, GTK_STATE_PRELIGHT, &gc);
  gtk_widget_modify_bg (b, GTK_STATE_ACTIVE, &gc);
  gtk_color_selection_set_current_color (GTK_COLOR_SELECTION (lookup_widget (gl->color_selector, 
									     "colorselection1")),
					 &gc);
}

static void gl_init (GfkGl * gl)
{
  GtkWidget * b = gtk_notebook_new ();

  gtk_widget_show (b);
  gl->params = gtk_window_new (GTK_WINDOW_TOPLEVEL);

  g_signal_connect (G_OBJECT (gl->params), "delete_event", G_CALLBACK (hide_params), NULL);
  g_object_set_data (G_OBJECT (gl->params), "GfkGl", gl);
  g_object_set_data (G_OBJECT (gl->params), "book", b);
  gtk_container_add (GTK_CONTAINER (gl->params), b);

  gl->properties = create_gl_params ();
  gfk_gl_prepend_params (gl, gl->properties, gtk_label_new ("Properties"));
  gtk_widget_show (gl->properties);

  gl->font = create_font_params ();
  gfk_gl_prepend_params (gl, gl->font, gtk_label_new ("Font"));
  /* dot not show by default */

  gl->color_selector = create_color_selector ();
  g_signal_connect (G_OBJECT (gl->color_selector), "delete_event", 
		    G_CALLBACK (gtk_widget_hide_on_delete), NULL);
  g_object_set_data (G_OBJECT (gl->color_selector), "GfkGl", gl);

  gl->gl = gfs_gl_new (GFK_GL_CLASS (GTS_OBJECT (gl)->klass)->gl_class);
  
  gfk_gl_set_color (gl, gl->gl->lc);

  GtkSpinButton * linewidth = lookup_widget_params (gl->properties, "linewidth");
  GtkAdjustment * a = gtk_spin_button_get_adjustment (linewidth);
  GLint range[2];
  glGetIntegerv (GL_ALIASED_LINE_WIDTH_RANGE, range);
  a->lower = range[0];
  a->upper = range[1];
  gtk_spin_button_set_adjustment (linewidth, a);
}

GfkGlClass * gfk_gl_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_info = {
      "GfkGl",
      sizeof (GfkGl),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_class_init,
      (GtsObjectInitFunc) gl_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &gfk_gl_info);
  }

  return klass;
}

GfkGl * gfk_gl_new (GfkGlClass * klass,
		    GtkWidget * glarea,
		    GtkWidget * list)
{
  GfkGl * object;

  g_return_val_if_fail (klass != NULL, NULL);
  g_return_val_if_fail (glarea != NULL, NULL);
  g_return_val_if_fail (list != NULL, NULL);

  object = GFK_GL (gts_object_new (GTS_OBJECT_CLASS (klass)));
  object->glarea = glarea;
  object->gl->p = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");
  object->list = list;
  (* klass->post_init) (object);

  return object;
}

void gfk_gl_expose (GfkGl * gl)
{
  g_return_if_fail (gl != NULL);

  if (gl->glarea->window)
    gdk_window_invalidate_rect (gl->glarea->window, &gl->glarea->allocation, FALSE);
}

void gfk_gl_set_sensitive (GfkGl * gl, GtkWidget * page, gboolean sensitive)
{
  GtkNotebook * book;

  g_return_if_fail (gl != NULL);
  g_return_if_fail (page != NULL);

  book = g_object_get_data (G_OBJECT (gl->params), "book");
  gtk_widget_set_sensitive (page, sensitive);
  gtk_widget_set_sensitive (gtk_notebook_get_tab_label (book, page), sensitive);
}

void gfk_gl_set_simulation (GfkGl * gl, GfsSimulation * sim)
{
  g_return_if_fail (gl != NULL);
  g_return_if_fail (sim != NULL);

  if (GFK_GL_CLASS (GTS_OBJECT (gl)->klass)->set_simulation)
    (* GFK_GL_CLASS (GTS_OBJECT (gl)->klass)->set_simulation) (gl, sim);
}

void gfk_gl_update_properties (GfkGl * gl)
{
  GtkTreeSelection * select = gtk_tree_view_get_selection (GTK_TREE_VIEW (gl->list));
  GtkTreeModel * model;
  GtkTreeIter iter;

  if (gtk_tree_selection_get_selected (select, &model, &iter)) {
    GfkGl * sgl;

    gtk_tree_model_get (model, &iter, GL_COLUMN, &sgl, -1);
    if (sgl == gl)
      gtk_list_store_set (GTK_LIST_STORE (model), &iter, 
			  PROPERTIES_COLUMN, gfk_gl_get_properties (gl), 
			  -1);
  }
}

void gfk_gl_update_interface (GfkGl * gl)
{
  g_return_if_fail (gl != NULL);

  if (GFK_GL_CLASS (GTS_OBJECT (gl)->klass)->update_interface)
    (* GFK_GL_CLASS (GTS_OBJECT (gl)->klass)->update_interface) (gl);
}

static void gfk_gl_menu_append (GfkGlClass * klass,
				GtkMenu * objects,
				GtkWidget * list,
				GtkWidget * glarea)
{
  GtkWidget * menuitem;

  menuitem = gtk_menu_item_new_with_mnemonic ((* klass->name) (klass));
  g_signal_connect (G_OBJECT (menuitem), "activate", G_CALLBACK (gl_add), NULL);
  g_object_set_data (G_OBJECT (menuitem), "klass", klass);
  g_object_set_data (G_OBJECT (menuitem), "list", list);
  g_object_set_data (G_OBJECT (menuitem), "glarea", glarea);
  gtk_widget_show (menuitem);
  gtk_container_add (GTK_CONTAINER (objects), menuitem);
}

static void gfk_gl_tool_append (GfkGlClass * klass,
				GtkToolbar * toolbar,
				GtkMenu * objects,
				GtkWidget * list,
				GtkWidget * glarea)
{
  GtkWidget * tool;

  gfk_gl_menu_append (klass, objects, list, glarea);
  tool = gtk_toolbar_append_element (toolbar,
				     GTK_TOOLBAR_CHILD_BUTTON,
				     NULL,
				     (* klass->name) (klass),
				     NULL, NULL,
				     (* klass->icon) (klass),
				     (GtkSignalFunc) gl_add, NULL);
  g_object_set_data (G_OBJECT (tool), "klass", klass);
  g_object_set_data (G_OBJECT (tool), "list", list);
  g_object_set_data (G_OBJECT (tool), "glarea", glarea);
  gtk_widget_show (tool);
}

/* GfkGlLabel: Object */

static void gl_label_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfk_gl_label_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsGlLabel * gl = GFS_GL_LABEL (GFK_GL (*o)->gl);
  GtkWidget * label = GFK_GL_LABEL (*o)->label;
  gtk_entry_set_text (lookup_widget_params (label, "labelentry"), gl->label);
  gtk_spin_button_set_value (lookup_widget_params (label, "labelx"), gl->p.x);
  gtk_spin_button_set_value (lookup_widget_params (label, "labely"), gl->p.y);
  gtk_spin_button_set_value (lookup_widget_params (label, "labelz"), gl->p.z);
  if (gl->symbol) {
    gl->symbol = FALSE;
    gtk_toggle_button_set_active (lookup_widget_params (label, "label_symbol_check"), TRUE);
  }
}

typedef struct {
  GfsSimulation * sim;
  GtkAdjustment * pos[3];
} BoundPar;

static void cell_bounds (FttCell * cell, BoundPar * par)
{
  static FttVector d[8] = {{ 1.,1.,1.},{ 1.,1.,-1.},{ 1.,-1.,1.},{ 1.,-1.,-1.},
			   {-1.,1.,1.},{-1.,1.,-1.},{-1.,-1.,1.},{-1.,-1.,-1.}};
  FttVector p, o;
  gdouble h = ftt_cell_size (cell)/2.;
  guint i;

  ftt_cell_pos (cell, &o);
  for (i = 0; i < 8; i++) {
    p.x = o.x + h*d[i].x; p.y = o.y + h*d[i].y; p.z = o.z + h*d[i].z;
    gfs_simulation_map_inverse (par->sim, &p);
    FttComponent c;
    for (c = 0; c < 3; c++) {
      gdouble e = (&p.x)[c];
      if (e > par->pos[c]->upper) par->pos[c]->upper = e;
      if (e < par->pos[c]->lower) par->pos[c]->lower = e;
    }
  }
}

static void gfk_gl_label_update_bounds (GfkGlLabel * gl)
{
  GtkSpinButton * label[3];
  BoundPar p;
  p.sim = GFK_GL (gl)->gl->sim;
  FttComponent c;
  for (c = 0; c < 3; c++) {
    const char * name[] = { "labelx", "labely", "labelz" };
    label[c] = lookup_widget_params (gl->label, name[c]);
    p.pos[c] = gtk_spin_button_get_adjustment (label[c]);
    p.pos[c]->lower = G_MAXDOUBLE;
    p.pos[c]->upper = -G_MAXDOUBLE;
  }
  gfs_domain_cell_traverse (GFS_DOMAIN (GFK_GL (gl)->gl->sim), 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS | FTT_TRAVERSE_LEVEL, 5,
			    (FttCellTraverseFunc) cell_bounds, &p);
  for (c = 0; c < 3; c++) {
    gdouble range = p.pos[c]->upper - p.pos[c]->lower;
    p.pos[c]->lower -= range/2.;
    p.pos[c]->upper += range/2.;
    gtk_spin_button_set_adjustment (label[c], p.pos[c]);
    if (p.pos[c]->value < p.pos[c]->lower)
      gtk_spin_button_set_value (label[c], p.pos[c]->lower);
    else if (p.pos[c]->value > p.pos[c]->upper)
      gtk_spin_button_set_value (label[c], p.pos[c]->upper);
    gdouble max = fabs (p.pos[c]->upper);
    if (fabs (p.pos[c]->lower) > max)
      max = fabs (p.pos[c]->lower);
    if (max == 0.)
      gtk_spin_button_set_digits (label[c], 3);
    else {
      gint n = 6. - log10 (max);
      gtk_spin_button_set_digits (label[c], MAX (n, 0));
    }
  }
}

static void gl_label_set_simulation (GfkGl * gl, GfsSimulation * sim)
{
  (* GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_label_class ())->parent_class)->set_simulation) 
    (gl, sim);
  gfk_gl_label_update_bounds (GFK_GL_LABEL (gl));
}

static gchar * gl_label_name (GfkGlClass * klass)
{
  static gchar name[] = "Label";
  return name;
}

static GtkWidget * gl_label_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "label-16x16.png");
}

#define MAXLENGTH 16

static gchar * gl_label_properties (GfkGl * gl)
{
  g_free (gl->props);
  gl->props = g_strndup (GFS_GL_LABEL (gl->gl)->label, MAXLENGTH + 3);
  if (strlen (gl->props) > MAXLENGTH) {
    guint i;
    for (i = 0; i < 3; i++)
      gl->props[MAXLENGTH + i] = '.';
  }
  return gl->props;
}

static void gl_label_class_init (GfkGlClass * klass)
{
  klass->gl_class = gfs_gl_label_class ();
  klass->name = gl_label_name;
  klass->icon = gl_label_icon;
  klass->properties = gl_label_properties;
  klass->set_simulation = gl_label_set_simulation;
  GTS_OBJECT_CLASS (klass)->read = gl_label_read;
}

static void gl_label_init (GfkGl * gl)
{
  gtk_widget_show (gl->font);

  GFK_GL_LABEL (gl)->label = create_label_params ();
  gfk_gl_prepend_params (gl, GFK_GL_LABEL (gl)->label, gtk_label_new ("Label"));
  gtk_widget_show (GFK_GL_LABEL (gl)->label);
}

GfkGlClass * gfk_gl_label_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_label_info = {
      "GfkGlLabel",
      sizeof (GfkGlLabel),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_label_class_init,
      (GtsObjectInitFunc) gl_label_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_class ()), &gfk_gl_label_info);
  }

  return klass;
}

/* GfkGl2D: Object */

static void gl2D_destroy (GtsObject * o)
{
  g_free (GFK_GL2D (o)->pickinfo);

  (* GTS_OBJECT_CLASS (gfk_gl2D_class ())->parent_class->destroy) (o);
}

static void gl2D_class_init (GfkGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gl2D_destroy;
  GTS_OBJECT_CLASS (klass)->read = gl2D_read;
  klass->post_init = gl2D_post_init;
  klass->set_simulation = gl2D_set_simulation;
  klass->pickinfo = gl2D_pickinfo;

  klass->gl_class = gfs_gl_symmetry_class ();
  klass->name = gl_symmetry_name;
  klass->icon = gl_symmetry_icon;
}

GfkGlClass * gfk_gl2D_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl2D_info = {
      "GfkGl2D",
      sizeof (GfkGl2D),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl2D_class_init,
      (GtsObjectInitFunc) gl2D_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_class ()), &gfk_gl2D_info);
  }

  return klass;
}

/* GfkGlCells: Object */

static gchar * gl_cells_name (GfkGlClass * klass)
{
  static gchar name[] = "Cells";
  return name;
}

static GtkWidget * gl_cells_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "cells-16x16.png");
}

static gboolean coarsenable (FttCell * cell)
{
  return TRUE;
}

static gchar * gl_cells_pickinfo (GfkGl * gl, gboolean motion)
{
  GfkGlCells * glk = GFK_GL_CELLS (gl);

  if (gtk_toggle_button_get_active (lookup_widget_params (glk->cells, "edit"))) {
    FttCell * cell = GFS_GL2D (gl->gl)->picked;
    guint level = gtk_spin_button_get_value (lookup_widget_params (glk->cells, "level"));
    guint l = ftt_cell_level (cell);
    if (FTT_CELL_IS_LEAF (cell) && l < level) {
      GfsDomain * domain = GFS_DOMAIN (gl->gl->sim);
      while (l < level) {
	ftt_cell_refine_corners (cell, (FttCellInitFunc) gfs_cell_fine_init, domain);
	ftt_cell_refine_single (cell, (FttCellInitFunc) gfs_cell_fine_init, domain);
	if (++l < level)
	  cell = ftt_cell_locate (cell, GFS_GL2D (gl->gl)->pickedpos, -1);
      }
      gfk_gl_expose (gl);
    }
    else if ((FTT_CELL_IS_LEAF (cell) && l > level) ||
	     (!FTT_CELL_IS_LEAF (cell) && l >= level)) {
      while (l > level) {
	cell = ftt_cell_parent (cell);
	l--;
      }
      ftt_cell_coarsen (cell, (FttCellCoarsenFunc) coarsenable, NULL,
			(FttCellCleanupFunc) gfs_cell_cleanup, GFS_DOMAIN (gl->gl->sim));
      gfk_gl_expose (gl);
    }
  }
  return gl2D_pickinfo (gl, motion);
}

static void gl_cells_class_init (GfkGlClass * klass)
{
  klass->gl_class = gfs_gl_cells_class ();
  klass->name = gl_cells_name;
  klass->icon = gl_cells_icon;
  klass->pickinfo = gl_cells_pickinfo;
}

static void gl_cells_init (GfkGl * gl)
{
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading"), FALSE);

  GFK_GL_CELLS (gl)->cells = create_cells_params ();
  gfk_gl_prepend_params (gl, GFK_GL_CELLS (gl)->cells, gtk_label_new ("Cells"));
  gtk_widget_show (GFK_GL_CELLS (gl)->cells);
}

GfkGlClass * gfk_gl_cells_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_cells_info = {
      "GfkGlCells",
      sizeof (GfkGlCells),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_cells_class_init,
      (GtsObjectInitFunc) gl_cells_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl2D_class ()), &gfk_gl_cells_info);
  }

  return klass;
}

/* GfkGlFractions: Object */

static gchar * gl_fractions_name (GfkGlClass * klass)
{
  static gchar name[] = "Fractions";
  return name;
}

static GtkWidget * gl_fractions_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "fractions-16x16.png");
}

static void gl_fractions_class_init (GfkGlClass * klass)
{
  klass->gl_class = gfs_gl_fractions_class ();
  klass->name = gl_fractions_name;
  klass->icon = gl_fractions_icon;
}

GfkGlClass * gfk_gl_fractions_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_fractions_info = {
      "GfkGlFractions",
      sizeof (GfkGl),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_fractions_class_init,
      (GtsObjectInitFunc) gl_fractions_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_class ()), &gfk_gl_fractions_info);
  }

  return klass;
}

/* GfkGlBoundaries: Object */

static gchar * gl_boundaries_name (GfkGlClass * klass)
{
  static gchar name[] = "Boundaries";
  return name;
}

static GtkWidget * gl_boundaries_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "boundaries-16x16.png");
}

static void gl_boundaries_class_init (GfkGlClass * klass)
{
  klass->gl_class = gfs_gl_boundaries_class ();
  klass->name = gl_boundaries_name;
  klass->icon = gl_boundaries_icon;
}

static void gl_boundaries_init (GfkGl * gl)
{
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading"), FALSE);
}

GfkGlClass * gfk_gl_boundaries_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_boundaries_info = {
      "GfkGlBoundaries",
      sizeof (GfkGl),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_boundaries_class_init,
      (GtsObjectInitFunc) gl_boundaries_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_class ()), &gfk_gl_boundaries_info);
  }

  return klass;
}

/* GfkGlLevels: Object */

static gchar * gl_levels_name (GfkGlClass * klass)
{
  static gchar name[] = "Levels";
  return name;
}

static GtkWidget * gl_levels_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "levels-16x16.png");
}

static void gl_levels_class_init (GfkGlClass * klass)
{
  klass->gl_class = gfs_gl_levels_class ();
  klass->name = gl_levels_name;
  klass->icon = gl_levels_icon;
}

static void gl_levels_init (GfkGl * gl)
{
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading"), FALSE);
}

GfkGlClass * gfk_gl_levels_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_levels_info = {
      "GfkGlLevels",
      sizeof (GfkGlLevels),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_levels_class_init,
      (GtsObjectInitFunc) gl_levels_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl2D_class ()),
				  &gfk_gl_levels_info);
  }

  return klass;
}

/* GfkFunction */

enum {
  COMPLETION_NAME_COLUMN,
  COMPLETION_DESCRIPTION_COLUMN,
  COMPLETION_POINTER_COLUMN
} CompletionColumns;

#define GLADE_HOOKUP_OBJECT(component,widget,name) \
  g_object_set_data_full (G_OBJECT (component), name, \
    gtk_widget_ref (widget), (GDestroyNotify) gtk_widget_unref)

static void on_function_activate (GtkEntry * entry, GfkFunctionSet set)
{
  GfkGl * gl = lookup_gl (entry);
  const gchar * needle = gtk_entry_get_text (entry);
  GtsFile * fp = (*set) (gl, needle);
  GtkEntryCompletion * c = gtk_entry_get_completion (entry);
  GtkTreeModel * list = gtk_entry_completion_get_model (c);
  GtkTreeIter iter;
  gboolean valid, found = FALSE;
  gchar * haystack;
  
  valid = gtk_tree_model_get_iter_first (list, &iter);
  while (valid && !found) {
    gtk_tree_model_get (list, &iter, 0, &haystack, -1);
    if (!strcmp (haystack, needle)) {
      found = TRUE;
      break;
    }
    valid = gtk_tree_model_iter_next (list, &iter);
  }

  if (fp) {
    GtkWidget * view = lookup_widget (gl->list, "view");
    GtkWidget * msg = gtk_message_dialog_new (GTK_WINDOW (view),
					      GTK_DIALOG_DESTROY_WITH_PARENT,
					      GTK_MESSAGE_WARNING,
					      GTK_BUTTONS_CLOSE,
					      "%s",
					      fp->error);
    gts_file_destroy (fp);
    gtk_dialog_run (GTK_DIALOG (msg));
    gtk_widget_destroy (msg);
    if (found)
      gtk_list_store_remove (GTK_LIST_STORE (list), &iter);
  }
  else {
    if (!found) {
      gtk_list_store_append (GTK_LIST_STORE (list), &iter);
      gtk_list_store_set (GTK_LIST_STORE (list), &iter, 
			  COMPLETION_NAME_COLUMN, needle, 
			  COMPLETION_DESCRIPTION_COLUMN, NULL,
			  COMPLETION_POINTER_COLUMN, NULL, 
			  -1);
    }
    gfk_gl_update_interface (gl);
    g_object_set_data (G_OBJECT (entry), "edited", NULL);
  }
}

static void variable_selected (GtkButton * button, GtkEntry * entry)
{
  GtkWidget * tree = lookup_widget (GTK_WIDGET (button), "tree");
  GtkTreeSelection * select = gtk_tree_view_get_selection (GTK_TREE_VIEW (tree));
  GtkTreeModel * model;
  GtkTreeIter iter;

  if (gtk_tree_selection_get_selected (select, &model, &iter)) {
    gchar * name;
    gtk_tree_model_get (model, &iter, 0, &name, -1);
    gtk_entry_set_text (entry, name);
    on_function_activate (entry, g_object_get_data (G_OBJECT (entry), "set"));
  }
}

static GtkWidget * variables_view (GtkEntry * entry, GtkTreeModel * model, GtkWindow * parent)
{
  GtkWidget * variables = create_variables ();

  gtk_window_set_transient_for (GTK_WINDOW (variables), parent);
  gtk_window_set_position (GTK_WINDOW (variables), GTK_WIN_POS_CENTER_ON_PARENT);
  g_signal_connect (G_OBJECT (variables), "delete_event", 
		    G_CALLBACK (gtk_widget_hide_on_delete), NULL);
  g_signal_connect_swapped (G_OBJECT (lookup_widget (variables, "cancel")), "clicked", 
			    G_CALLBACK (gtk_widget_hide), variables);
  g_signal_connect_swapped (G_OBJECT (lookup_widget (variables, "OK")), "clicked", 
			    G_CALLBACK (gtk_widget_hide), variables);
  g_signal_connect (G_OBJECT (lookup_widget (variables, "OK")), "clicked",
		    G_CALLBACK (variable_selected), entry);

  GtkWidget * tree = lookup_widget (variables, "tree");
  gtk_tree_view_set_model (GTK_TREE_VIEW (tree), model);
  GtkTreeSelection * select = gtk_tree_view_get_selection (GTK_TREE_VIEW (tree));
  gtk_tree_selection_set_mode (select, GTK_SELECTION_SINGLE);

  GtkCellRenderer * renderer = gtk_cell_renderer_text_new ();
  GtkTreeViewColumn * column = 
    gtk_tree_view_column_new_with_attributes ("Name", renderer,
					      "text", COMPLETION_NAME_COLUMN,
					      NULL);
  gtk_tree_view_column_set_resizable (column, TRUE);
  gtk_tree_view_column_set_sizing (column, GTK_TREE_VIEW_COLUMN_FIXED);
  gtk_tree_view_column_set_fixed_width (column, 80);
  gtk_tree_view_append_column (GTK_TREE_VIEW (tree), column);

  column = gtk_tree_view_column_new_with_attributes ("Description", renderer,
						     "text", COMPLETION_DESCRIPTION_COLUMN,
						     NULL);
  gtk_tree_view_column_set_resizable (column, TRUE);
  gtk_tree_view_column_set_sizing (column, GTK_TREE_VIEW_COLUMN_FIXED);
  gtk_tree_view_column_set_fixed_width (column, 80);
  gtk_tree_view_append_column (GTK_TREE_VIEW (tree), column);

  return variables;
}

static gboolean match_selected (GtkEntryCompletion * widget,
				GtkTreeModel * model,
				GtkTreeIter * iter)
{
  GtkWidget * entry = gtk_entry_completion_get_entry (widget);
  gchar * text;
  gtk_tree_model_get (model, iter, COMPLETION_NAME_COLUMN, &text, -1);
  gtk_entry_set_text (GTK_ENTRY (entry), text);
  g_signal_emit_by_name (entry, "activate", entry, NULL);
  return TRUE;
}

static GtkWidget * gfk_function (GtkWidget * parent,
				 GfkGl * gl, 
				 const gchar * name,
				 GfkFunctionSet set)
{
  GtkWidget *gfk_function;
  GtkWidget *scalar;
  GtkWidget *select;

  gfk_function = gtk_hbox_new (FALSE, 0);

  scalar = gtk_entry_new ();
  gtk_widget_show (scalar);
  gtk_box_pack_start (GTK_BOX (gfk_function), scalar, TRUE, TRUE, 0);
  gtk_entry_set_invisible_char (GTK_ENTRY (scalar), 9679);
  g_object_set_data (G_OBJECT (scalar), "set", set);

  GtkEntryCompletion * c = gtk_entry_completion_new ();
  GtkTreeModel * completion = g_object_get_data (G_OBJECT (lookup_widget (gl->list, "view")),
						 "completion");
  gtk_entry_completion_set_model (c, completion);
  gtk_entry_completion_set_text_column (c, COMPLETION_NAME_COLUMN);
  g_signal_connect (G_OBJECT (c), "match-selected", G_CALLBACK (match_selected), NULL);
  gtk_entry_set_completion (GTK_ENTRY (scalar), c);

  select = gtk_button_new_with_mnemonic ("...");
  gtk_widget_show (select);
  gtk_box_pack_start (GTK_BOX (gfk_function), select, FALSE, FALSE, 0);

  GtkWidget * variables = variables_view (GTK_ENTRY (scalar), completion, GTK_WINDOW (gl->params));
  GLADE_HOOKUP_OBJECT (select, variables, "variables");
  g_signal_connect_swapped ((gpointer) select, "clicked",
			    G_CALLBACK (gtk_widget_show), variables);

  g_signal_connect ((gpointer) scalar, "changed",
                    G_CALLBACK (mark_edited),
                    NULL);
  g_signal_connect ((gpointer) scalar, "move_cursor",
                    G_CALLBACK (mark_edited),
                    NULL);
  g_signal_connect ((gpointer) scalar, "activate",
                    G_CALLBACK (on_function_activate),
                    set);

  /* Store pointers to all widgets, for use by lookup_widget(). */
  GLADE_HOOKUP_OBJECT (parent, gfk_function, "gfk_function");
  GLADE_HOOKUP_OBJECT (parent, scalar, name);
  GLADE_HOOKUP_OBJECT (parent, select, "select");

  return gfk_function;
}

/* GfkGlScalar: Object */

static void gl_scalar_read (GtsObject ** o, GtsFile * fp)
{
  GtkWidget * scalar = GFK_GL_SCALAR (*o)->scalar;
  GfsGlScalar * gl = GFS_GL_SCALAR (GFK_GL (*o)->gl);
  GtkWidget * amin = lookup_widget_params (scalar, "amin");
  GtkWidget * amax = lookup_widget_params (scalar, "amax");
  GtkSpinButton * spinbuttonmin = lookup_widget_params (scalar, "spinbuttonmin");
  GtkSpinButton * spinbuttonmax = lookup_widget_params (scalar, "spinbuttonmax");
  gdouble maxv = gtk_spin_button_get_value (spinbuttonmax);

  (* GTS_OBJECT_CLASS (gfk_gl_scalar_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (!gl->amin) {
    gl->amin = TRUE;
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (amin), FALSE);
    if (gl->amax && gl->min > maxv)
      gtk_spin_button_set_value (spinbuttonmax, gl->min);
    gtk_spin_button_set_value (spinbuttonmin, gl->min);
  }

  if (!gl->amax) {
    gl->amax = TRUE;
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (amax), FALSE);
    gtk_spin_button_set_value (spinbuttonmax, gl->max);
  }

  if (gl->show) {
    gl->show = FALSE;
    gtk_toggle_button_set_active 
      (GTK_TOGGLE_BUTTON (lookup_widget_params (scalar, "show")), TRUE);
  }

  GtkOptionMenu * colormap = lookup_widget_params (scalar, "colormap");
  if (colormap) {
    if (!strcmp (gl->cmap->name, "Jet"))
      gtk_option_menu_set_history (colormap, 0);
    else if (!strcmp (gl->cmap->name, "Cool"))
      gtk_option_menu_set_history (colormap, 1);
    else if (!strcmp (gl->cmap->name, "Gray"))
      gtk_option_menu_set_history (colormap, 2);
  }
}

guint gfk_decimal_digits (double x, guint significant)
{
  if (x == 0.)
    return significant;
  else {
    gdouble n = significant - 1. - floor (log10 (fabs (x)));
    return n > 0 ? n : 0;
  }
}

static void set_value (GtkSpinButton * s, gdouble v)
{
  GtkAdjustment * as = gtk_spin_button_get_adjustment (s);
  if (as->lower > v)
    as->lower = v;
  if (as->upper < v)
    as->upper = v;
  gtk_spin_button_set_adjustment (s, as);
  gtk_spin_button_set_value (s, v);
}

static void gl_scalar_update_interface (GfkGl * g)
{
  GfkGlScalar * gl = GFK_GL_SCALAR (g);
  GfsGlScalar * gls = GFS_GL_SCALAR (GFK_GL (gl)->gl);
  GtkSpinButton * min = lookup_widget_params (gl->scalar, "spinbuttonmin");
  GtkSpinButton * max = lookup_widget_params (gl->scalar, "spinbuttonmax");
  GtkEntry * scalar = lookup_widget_params (gl->scalar, "scalar");

  if (!GFK_IS_EDITED (scalar))
    gtk_entry_set_text (scalar, gls->expr->str);

  if (gls->amin)
    set_value (min, gls->aminv);
  else if (!GFK_IS_EDITED (min))
    set_value (min, gls->min);
  if (gls->amax)
    set_value (max, gls->amaxv);
  else if (!GFK_IS_EDITED (max))
    set_value (max, gls->max);

  if (gls->max > gls->min) {
    gdouble step = (gls->max - gls->min)/1000.;
    guint digits = gfk_decimal_digits (step, 2);
    GtkAdjustment * amin = gtk_spin_button_get_adjustment (min);
    GtkAdjustment * amax = gtk_spin_button_get_adjustment (max);
    amin->step_increment = amax->step_increment = step;
    amin->page_increment = amax->page_increment = step*10.;
    if (!GFK_IS_EDITED (min))
      gtk_spin_button_configure (min, amin, 2.*step, digits);
    if (!GFK_IS_EDITED (max))
      gtk_spin_button_configure (max, amax, 2.*step, digits);
  }
  gfk_gl_update_properties (GFK_GL (gl));
  gfk_gl_expose (GFK_GL (gl));  
}

static void gl_scalar_set_simulation (GfkGl * object, GfsSimulation * sim)
{
  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_scalar_class ())->parent_class)->set_simulation) 
    (object, sim);
  
  gl_scalar_update_interface (object);
}

static GtsFile * gl_scalar_set (GfkGl * gl, const gchar * s)
{
  GtsFile * fp = gfs_gl_scalar_set (GFS_GL_SCALAR (gl->gl), s);
  if (!fp)
    gl_scalar_update_interface (gl);
  return fp;
}

static void gl_scalar_post_init (GfkGl * object)
{
  GfkGlScalar * gls = GFK_GL_SCALAR (object);

  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_scalar_class ())->parent_class)->post_init) (object);

  GtkWidget * scalar = gfk_function (gls->scalar, object, "scalar", gl_scalar_set);
  gtk_widget_show (scalar);
  gtk_table_attach (GTK_TABLE (lookup_widget_params (gls->scalar, "table1")),
		    scalar, 1, 3, 0, 1,
                    (GtkAttachOptions) (GTK_FILL),
                    (GtkAttachOptions) (0), 0, 0);

  gl_scalar_update_interface (object);
}

static gchar * gl_scalar_properties (GfkGl * gl)
{
  g_free (gl->props);
  gl->props = g_strndup (GFS_GL_SCALAR (gl->gl)->expr->str, MAXLENGTH + 3);
  if (strlen (gl->props) > MAXLENGTH) {
    guint i;
    for (i = 0; i < 3; i++)
      gl->props[MAXLENGTH + i] = '.';
  }
  return gl->props;
}

static gchar * gl_scalar_pickinfo (GfkGl * gl, gboolean motion)
{
  gchar * s = 
    (* GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_scalar_class ())->parent_class)->pickinfo) (gl, motion);
  gdouble val = GFS_VALUE (GFS_GL2D (gl->gl)->picked, 
			   GFS_GL_SCALAR (gl->gl)->v);
  gchar * s1 = (val != GFS_NODATA ? 
		g_strdup_printf ("%s %s %g", s, gl_scalar_properties (gl), val) :
		g_strdup_printf ("%s %s NODATA", s, gl_scalar_properties (gl))
		);
  g_free (GFK_GL2D (gl)->pickinfo);
  GFK_GL2D (gl)->pickinfo = s1;
  return s1;
}

static void gl_scalar_class_init (GfkGlClass * klass)
{
  klass->post_init = gl_scalar_post_init;
  klass->set_simulation = gl_scalar_set_simulation;
  klass->update_interface = gl_scalar_update_interface;
  klass->properties = gl_scalar_properties;
  klass->pickinfo = gl_scalar_pickinfo;
  GTS_OBJECT_CLASS (klass)->read = gl_scalar_read;
}

static void set_cmap_cool_warm (GtkWidget * cmap, GfkGl * gl)
{
  gfs_colormap_destroy (GFS_GL_SCALAR (gl->gl)->cmap);
  GFS_GL_SCALAR (gl->gl)->cmap = gfs_colormap_cool_warm ();
  gfk_gl_expose (gl);
}

static void set_cmap_gray (GtkWidget * cmap, GfkGl * gl)
{
  gfs_colormap_destroy (GFS_GL_SCALAR (gl->gl)->cmap);
  GFS_GL_SCALAR (gl->gl)->cmap = gfs_colormap_gray ();
  gfk_gl_expose (gl);
}

static void set_cmap_jet (GtkWidget * cmap, GfkGl * gl)
{
  gfs_colormap_destroy (GFS_GL_SCALAR (gl->gl)->cmap);
  GFS_GL_SCALAR (gl->gl)->cmap = gfs_colormap_jet ();
  gfk_gl_expose (gl);
}

static void gl_scalar_init (GfkGlScalar * object)
{
  GtkOptionMenu * colormap;

  object->scalar = create_scalar_params ();
  if ((colormap = lookup_widget_params (object->scalar, "colormap"))) {
    GtkWidget * m = gtk_menu_new (), * i;

    i = gtk_menu_item_new_with_label ("Jet");
    g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_cmap_jet), object);
    gtk_menu_append (m, i);
    gtk_widget_show (i);
    gtk_option_menu_set_menu (colormap, m);
    gtk_widget_show (m);

    i = gtk_menu_item_new_with_label ("Cool");
    g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_cmap_cool_warm), object);
    gtk_menu_append (m, i);
    gtk_widget_show (i);
    gtk_option_menu_set_menu (colormap, m);
    gtk_widget_show (m);

    i = gtk_menu_item_new_with_label ("Gray");
    g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_cmap_gray), object);
    gtk_menu_append (m, i);
    gtk_widget_show (i);
    gtk_option_menu_set_menu (colormap, m);
    gtk_widget_show (m);
  }
  gfk_gl_prepend_params (GFK_GL (object), object->scalar, gtk_label_new ("Scalar"));
  gtk_widget_show (object->scalar);
}

GfkGlClass * gfk_gl_scalar_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_scalar_info = {
      "GfkGlScalar",
      sizeof (GfkGlScalar),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_scalar_class_init,
      (GtsObjectInitFunc) gl_scalar_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl2D_class ()),
				  &gfk_gl_scalar_info);
  }

  return klass;
}

/* GfkGlSquares: Object */

static gchar * gl_squares_name (GfkGlClass * klass)
{
  static gchar name[] = "Squares";
  return name;
}

static GtkWidget * gl_squares_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "squares-16x16.png");
}

static void gl_squares_class_init (GfkGlClass * klass)
{
  klass->gl_class = gfs_gl_squares_class ();
  klass->name = gl_squares_name;
  klass->icon = gl_squares_icon;
}

static void gl_squares_init (GfkGl * gl)
{
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color"), FALSE);
}

GfkGlClass * gfk_gl_squares_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_squares_info = {
      "GfkGlSquares",
      sizeof (GfkGlSquares),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_squares_class_init,
      (GtsObjectInitFunc) gl_squares_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_scalar_class ()),
				  &gfk_gl_squares_info);
  }

  return klass;
}

/* GfkGlLinear: Object */

static void update_linear_color (GfkGl * gl)
{
  gfk_gl_set_sensitive (GFK_GL (gl), GFK_GL_SCALAR (gl)->scalar, FALSE);
  gfk_gl_update_properties (GFK_GL (gl));
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color_label"), 
			    TRUE);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color"), TRUE);
}

static void gl_linear_read (GtsObject ** o, GtsFile * fp)
{
  GtkWidget * scalar = GFK_GL_LINEAR (*o)->scalar;
  GfsGlLinear * gl = GFS_GL_LINEAR (GFK_GL (*o)->gl);

  (* GTS_OBJECT_CLASS (gfk_gl_linear_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (gl->reversed) {
    gl->reversed = FALSE;
    gtk_toggle_button_set_active (lookup_widget_params (scalar, "reversed"), TRUE);
  }

  if (!gl->use_scalar) {
    gtk_option_menu_set_history (lookup_widget_params (scalar, "color"), 0);
    update_linear_color (GFK_GL (*o));
  }
}

static gchar * gl_linear_properties (GfkGl * gl)
{
  GfsGlLinear * gll = GFS_GL_LINEAR (gl->gl);
  gchar * s;
  if (gll->vf->f && gfs_function_get_constant_value (gll->vf->f) != 0.)
    s = g_strjoin (" ",
		   gll->expr->str,
		   gll->use_scalar ?
		   ((* gfk_gl_scalar_class ()->properties) (gl)) : NULL, 
		   NULL);
  else
    s = g_strjoin (" ",
		   gll->use_scalar ?
		   ((* gfk_gl_scalar_class ()->properties) (gl)) : NULL, 
		   NULL);
  g_free (gl->props);
  gl->props = s;
  return gl->props;
}

static gchar * gl_linear_name (GfkGlClass * klass)
{
  static gchar name[] = "Linear";
  return name;
}

static GtkWidget * gl_linear_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "linear-16x16.png");
}

static gchar * gl_linear_pickinfo (GfkGl * gl, gboolean motion)
{
  gchar * s = 
    (* GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_scalar_class ())->parent_class)->pickinfo) 
    (gl, motion);
  gchar * s1 = GFS_VALUE (GFS_GL2D (gl->gl)->picked, GFS_GL_SCALAR (gl->gl)->v) != GFS_NODATA ?
    g_strdup_printf ("%s %s %g", s, gl_scalar_properties (gl),
		     gfs_interpolate (GFS_GL2D (gl->gl)->picked,
				      GFS_GL2D (gl->gl)->pickedpos,
				      GFS_GL_SCALAR (gl->gl)->v)) :
    g_strdup_printf ("%s %s NODATA", s, gl_scalar_properties (gl));
  g_free (GFK_GL2D (gl)->pickinfo);
  GFK_GL2D (gl)->pickinfo = s1;
  return s1;
}

static void gl_linear_update_interface (GfkGl * g)
{
  GfkGlLinear * gl = GFK_GL_LINEAR (g);
  GfsGlLinear * gli = GFS_GL_LINEAR (GFK_GL (gl)->gl);
  GtkEntry * scalar = lookup_widget_params (gl->scalar, "scalar");

  if (!GFK_IS_EDITED (scalar))
    gtk_entry_set_text (scalar, gli->expr->str);

  gfk_gl_update_properties (GFK_GL (gl));
  gfk_gl_expose (GFK_GL (gl));
}

static void gl_linear_set_simulation (GfkGl * object, GfsSimulation * sim)
{
  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_linear_class ())->parent_class)->set_simulation)
    (object, sim);
  
  gl_linear_update_interface (object);
}

static void set_linear_scalar (GtkWidget * scalar, GfkGl * gl)
{
  GFS_GL_LINEAR (gl->gl)->use_scalar = GFS_GL_SCALAR (gl->gl)->v;
  gfk_gl_set_sensitive (gl, GFK_GL_SCALAR (gl)->scalar, TRUE);
  gfk_gl_update_properties (gl);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color"), FALSE);
  gfk_gl_expose (gl);
}

static void set_linear_color (GtkWidget * color, GfkGlLinear * gl)
{
  GFS_GL_LINEAR (GFK_GL (gl)->gl)->use_scalar = NULL;
  update_linear_color (GFK_GL (gl));
  gfk_gl_expose (GFK_GL (gl));
}

static GtsFile * gl_linear_set (GfkGl * gl, const char * s)
{
  return gfs_gl_linear_set (GFS_GL_LINEAR (gl->gl), s);
}

static void gl_linear_post_init (GfkGl * object)
{
  GfkGlLinear * gli = GFK_GL_LINEAR (object);

  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_linear_class ())->parent_class)->post_init)
    (object);

  GtkWidget * scalar = gfk_function (gli->scalar, object, "scalar", gl_linear_set);
  gtk_widget_show (scalar);
  gtk_table_attach (GTK_TABLE (lookup_widget_params (gli->scalar, "table1")),
		    scalar, 1, 3, 0, 1,
                    (GtkAttachOptions) (GTK_FILL),
                    (GtkAttachOptions) (0), 0, 0);

  gl_linear_update_interface (object);
}

static void gl_linear_class_init (GfkGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_linear_read;
  klass->gl_class = gfs_gl_linear_class ();
  klass->set_simulation = gl_linear_set_simulation;
  klass->update_interface = gl_linear_update_interface;
  klass->post_init = gl_linear_post_init;
  klass->properties = gl_linear_properties;
  klass->name = gl_linear_name;
  klass->icon = gl_linear_icon;
  klass->pickinfo = gl_linear_pickinfo;
}

static void gl_linear_init (GfkGl * gl)
{
  GfkGlLinear * gls = GFK_GL_LINEAR (gl);
  gls->scalar = create_linear_params ();

  GtkWidget * m = gtk_menu_new ();
  GtkWidget * i = gtk_menu_item_new_with_label ("Default");
  g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_linear_color), gl);
  gtk_menu_append (m, i);
  gtk_widget_show (i);

  i = gtk_menu_item_new_with_label ("Scalar");
  g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_linear_scalar), gl);
  gtk_menu_append (m, i);
  gtk_widget_show (i);
  
  gtk_option_menu_set_menu (lookup_widget_params (gls->scalar, "color"), m);
  gtk_option_menu_set_history (lookup_widget_params (gls->scalar, "color"), 1);
  gtk_widget_show (m);  
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color"), FALSE);  

  gfk_gl_prepend_params (gl, gls->scalar, gtk_label_new ("Linear"));
#if FTT_2D
  gtk_widget_show (gls->scalar);
#endif
}

GfkGlClass * gfk_gl_linear_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_linear_info = {
      "GfkGlLinear",
      sizeof (GfkGlLinear),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_linear_class_init,
      (GtsObjectInitFunc) gl_linear_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_scalar_class ()),
				  &gfk_gl_linear_info);
  }

  return klass;
}

/* GfkGlIsoline: Object */

static void gl_isoline_read (GtsObject ** o, GtsFile * fp)
{
  GtkWidget * scalar = GFK_GL_SCALAR (*o)->scalar;
  GfsGlIsoline * gl = GFS_GL_ISOLINE (GFK_GL (*o)->gl);

  (* GTS_OBJECT_CLASS (gfk_gl_isoline_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  gtk_spin_button_set_value (lookup_widget_params (scalar, "spinbuttoniso"), gl->n);
  if (gl->ls)
    gtk_entry_set_text (lookup_widget_params (scalar, "entryiso"), gl->ls);
}

static gchar * gl_isoline_name (GfkGlClass * klass)
{
  static gchar name[] = "Isoline";
  return name;
}

static GtkWidget * gl_isoline_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "isoline-16x16.png");
}

static void gl_isoline_class_init (GfkGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_isoline_read;
  klass->gl_class = gfs_gl_isoline_class ();
  klass->name = gl_isoline_name;
  klass->icon = gl_isoline_icon;
}

static void gl_isoline_init (GfkGl * gl)
{
  GfkGlScalar * gls = GFK_GL_SCALAR (gl);
  gtk_widget_destroy (gls->scalar);
  gls->scalar = create_isoline_params ();
  gfk_gl_prepend_params (gl, gls->scalar, gtk_label_new ("Isoline"));
  gtk_widget_show (gls->scalar);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading"), FALSE);

  gtk_option_menu_set_history (lookup_widget_params (GFK_GL_LINEAR (gl)->scalar, "color"), 0);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL_LINEAR (gl)->scalar, "normals"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL_LINEAR (gl)->scalar, "reversed"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color_label"), 
			    TRUE);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color"), TRUE);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL_LINEAR (gl)->scalar, "color"), FALSE);
}

GfkGlClass * gfk_gl_isoline_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_isoline_info = {
      "GfkGlIsoline",
      sizeof (GfkGlIsoline),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_isoline_class_init,
      (GtsObjectInitFunc) gl_isoline_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_linear_class ()),
				  &gfk_gl_isoline_info);
  }

  return klass;
}

/* GfkGlVOF: Object */

static void update_vof_scalar (GfkGl * gl)
{
  gfk_gl_set_sensitive (gl, GFK_GL_SCALAR (gl)->scalar, TRUE);
  gfk_gl_update_properties (gl);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color"), FALSE);
}

static void gl_vof_read (GtsObject ** o, GtsFile * fp)
{
  GtkWidget * scalar = GFK_GL_VOF (*o)->scalar;
  GfsGlVOF * gl = GFS_GL_VOF (GFK_GL (*o)->gl);

  (* GTS_OBJECT_CLASS (gfk_gl_vof_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (gl->reversed) {
    gl->reversed = FALSE;
    gtk_toggle_button_set_active (lookup_widget_params (scalar, "reversed"), TRUE);
  }

  if (gl->draw_edges) {
    gl->draw_edges = FALSE;
    gtk_toggle_button_set_active (lookup_widget_params (scalar, "visible"), TRUE);
  }

  if (gl->use_scalar) {
    gtk_option_menu_set_history (lookup_widget_params (scalar, "color"), gl->interpolate ? 2 : 1);
    update_vof_scalar (GFK_GL (*o));
  }
}

static gchar * gl_vof_name (GfkGlClass * klass)
{
  static gchar name[] = "VOF";
  return name;
}

static GtkWidget * gl_vof_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "vof-16x16.png");
}

static void gl_vof_update_interface (GfkGl * g)
{
  GfkGlVOF * gl = GFK_GL_VOF (g);
  GfsGlVOF * gli = GFS_GL_VOF (GFK_GL (gl)->gl);
  GtkEntry * scalar = lookup_widget_params (gl->scalar, "scalar");

  if (!GFK_IS_EDITED (scalar))
    gtk_entry_set_text (scalar, gli->expr->str);

  gfk_gl_update_properties (GFK_GL (gl));
  gfk_gl_expose (GFK_GL (gl));
}

static void gl_vof_set_simulation (GfkGl * object, GfsSimulation * sim)
{
  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_vof_class ())->parent_class)->set_simulation)
    (object, sim);
  
  gl_vof_update_interface (object);
}

static void set_vof_scalar (GtkWidget * scalar, GfkGl * gl)
{
  GFS_GL_VOF (gl->gl)->use_scalar = GFS_GL_SCALAR (gl->gl)->v;
  GFS_GL_VOF (gl->gl)->interpolate = FALSE;
  update_vof_scalar (gl);
  gfk_gl_expose (gl);
}

static void set_vof_scalar_interpolated (GtkWidget * scalar, GfkGl * gl)
{
  GFS_GL_VOF (gl->gl)->use_scalar = GFS_GL_SCALAR (gl->gl)->v;
  GFS_GL_VOF (gl->gl)->interpolate = TRUE;
  update_vof_scalar (gl);
  gfk_gl_expose (gl);
}

static void set_vof_color (GtkWidget * color, GfkGlVOF * gl)
{
  GFS_GL_VOF (GFK_GL (gl)->gl)->use_scalar = NULL;
  gfk_gl_set_sensitive (GFK_GL (gl), GFK_GL_SCALAR (gl)->scalar, FALSE);
  gfk_gl_update_properties (GFK_GL (gl));
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color_label"), 
			    TRUE);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color"), TRUE);
  gfk_gl_expose (GFK_GL (gl));
}

static GtsFile * gl_vof_set (GfkGl * gl, const gchar * s)
{
  return gfs_gl_vof_set (GFS_GL_VOF (gl->gl), s);
}

static void gl_vof_post_init (GfkGl * object)
{
  GfkGlVOF * gli = GFK_GL_VOF (object);

  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_vof_class ())->parent_class)->post_init)
    (object);

  GtkWidget * scalar = gfk_function (gli->scalar, object, "scalar", gl_vof_set);
  gtk_widget_show (scalar);
  gtk_table_attach (GTK_TABLE (lookup_widget_params (gli->scalar, "table1")),
		    scalar, 1, 3, 0, 1,
                    (GtkAttachOptions) (GTK_FILL),
                    (GtkAttachOptions) (0), 0, 0);
  
  gl_vof_update_interface (object);
}

static gchar * gl_vof_properties (GfkGl * gl)
{
  gchar * s = g_strjoin (" ",
			 GFS_GL_VOF (gl->gl)->expr->str,
			 GFS_GL_VOF (gl->gl)->use_scalar ?
			 ((* gfk_gl_scalar_class ()->properties) (gl)) : NULL, 
			 NULL);
  g_free (gl->props);
  gl->props = s;
  return gl->props;
}

static void gl_vof_class_init (GfkGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_vof_read;
  klass->gl_class = gfs_gl_vof_class ();
  klass->set_simulation = gl_vof_set_simulation;
  klass->update_interface = gl_vof_update_interface;
  klass->post_init = gl_vof_post_init;
  klass->properties = gl_vof_properties;
  klass->name = gl_vof_name;
  klass->icon = gl_vof_icon;
}

static void gl_vof_init (GfkGl * gl)
{
  GfkGlVOF * gli = GFK_GL_VOF (gl);

  gli->scalar = create_vof_params ();

  gfk_gl_set_sensitive (gl, GFK_GL_SCALAR (gl)->scalar, FALSE);

  GtkWidget * m = gtk_menu_new ();
  GtkWidget * i = gtk_menu_item_new_with_label ("Default");
  g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_vof_color), gl);
  gtk_menu_append (m, i);
  gtk_widget_show (i);

  i = gtk_menu_item_new_with_label ("Scalar");
  g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_vof_scalar), gl);
  gtk_menu_append (m, i);
  gtk_widget_show (i);
  
  i = gtk_menu_item_new_with_label ("Scalar (interpolated)");
  g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_vof_scalar_interpolated), gl);
  gtk_menu_append (m, i);
  gtk_widget_show (i);
  
  gtk_option_menu_set_menu (lookup_widget_params (gli->scalar, "color"), m);
  gtk_widget_show (m);

#if FTT_2D
  gtk_widget_hide (lookup_widget_params (gli->scalar, "normals"));
  gtk_widget_hide (lookup_widget_params (gli->scalar, "reversed"));
  gtk_widget_hide (lookup_widget_params (gli->scalar, "edges"));
  gtk_widget_hide (lookup_widget_params (gli->scalar, "visible"));
#else
  gtk_widget_hide (GFK_GL2D (gl)->params);
#endif

  gfk_gl_prepend_params (gl, gli->scalar, gtk_label_new ("VOF"));
  gtk_widget_show (gli->scalar);
}

GfkGlClass * gfk_gl_vof_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_vof_info = {
      "GfkGlVOF",
      sizeof (GfkGlVOF),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_vof_class_init,
      (GtsObjectInitFunc) gl_vof_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_scalar_class ()),
				  &gfk_gl_vof_info);
  }

  return klass;
}

/* GfkGlVectors: Object */

static void set_vector_color (GtkWidget * color, GfkGlVectors * gl)
{
  GFS_GL_VECTORS (GFK_GL (gl)->gl)->use_scalar = FALSE;
  gfk_gl_set_sensitive (GFK_GL (gl), GFK_GL_SCALAR (gl)->scalar, FALSE);
  gfk_gl_update_properties (GFK_GL (gl));
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color_label"), 
			    TRUE);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color"), TRUE);
  gfk_gl_expose (GFK_GL (gl));
}

static void update_vector_scalar (GfkGlVectors * gl)
{
  gfk_gl_set_sensitive (GFK_GL (gl), GFK_GL_SCALAR (gl)->scalar, TRUE);
  gfk_gl_update_properties (GFK_GL (gl));
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color_label"), 
			    FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color"), 
			    FALSE);
}

static void set_vector_scalar (GtkWidget * scalar, GfkGlVectors * gl)
{
  GFS_GL_VECTORS (GFK_GL (gl)->gl)->use_scalar = TRUE;
  update_vector_scalar (gl);
  gfk_gl_expose (GFK_GL (gl));
}

static void gl_vectors_read (GtsObject ** o, GtsFile * fp)
{
  GtkWidget * vector = GFK_GL_VECTORS (*o)->vector;
  GfsGlVectors * gl = GFS_GL_VECTORS (GFK_GL (*o)->gl);

  (* GTS_OBJECT_CLASS (gfk_gl_vectors_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  gtk_spin_button_set_value (lookup_widget_params (vector, "scale"), gl->scale);

  if (gl->use_scalar) {
    gtk_option_menu_set_history (lookup_widget_params (vector, "color"), 1);
    update_vector_scalar (GFK_GL_VECTORS (*o));
  }
}

static void gl_vectors_update_interface (GfkGlVectors * gl)
{
  GfsGlVectors * glv = GFS_GL_VECTORS (GFK_GL (gl)->gl);
  GtkSpinButton * scale = lookup_widget_params (gl->vector, "scale");
  GtkAdjustment * ascale = gtk_spin_button_get_adjustment (scale);
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++) {
    static gchar * name[] = {"Vx", "Vy", "Vz"};
    GtkEntry * entry = lookup_widget_params (gl->vector, name[c]);
    if (!GFK_IS_EDITED (entry))
      gtk_entry_set_text (entry, glv->expr[c]->str);
  }

  ascale->upper = glv->max > 0. ? 10000.*glv->h/glv->max : 1.;
  ascale->step_increment = glv->max > 0. ? glv->h/glv->max/10. : 0.1;
  ascale->page_increment = 10.*ascale->step_increment;
  ascale->value = glv->scale;
  if (!GFK_IS_EDITED (scale))
    gtk_spin_button_configure (scale, ascale, 2.*ascale->step_increment,
			       gfk_decimal_digits (ascale->step_increment, 2));
  gfk_gl_update_properties (GFK_GL (gl));
  gfk_gl_expose (GFK_GL (gl));
}

static void gl_vectors_set_simulation (GfkGl * object, GfsSimulation * sim)
{
  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_vectors_class ())->parent_class)->set_simulation)
    (object, sim);
  
  gl_vectors_update_interface (GFK_GL_VECTORS (object));
}

static GtsFile * set_x (GfkGl * gl, const gchar * func)
{
  return gfs_gl_vectors_set (GFS_GL_VECTORS (gl->gl), 0, func);
}

static GtsFile * set_y (GfkGl * gl, const gchar * func)
{
  return gfs_gl_vectors_set (GFS_GL_VECTORS (gl->gl), 1, func);
}

static GtsFile * set_z (GfkGl * gl, const gchar * func)
{
  return gfs_gl_vectors_set (GFS_GL_VECTORS (gl->gl), 2, func);
}

static void gl_vectors_post_init (GfkGl * object)
{
  GfkGlVectors * gls = GFK_GL_VECTORS (object);
  FttComponent c;

  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_vectors_class ())->parent_class)->post_init) (object);

  GtkTable * table = GTK_TABLE (lookup_widget_params (gls->vector, "table2"));
  for (c = 0; c < FTT_DIMENSION; c++) {
    static gchar * name[] = {"Vx", "Vy", "Vz"};
    static GfkFunctionSet func[] = {set_x, set_y, set_z};
    
    GtkWidget * scalar = gfk_function (gls->vector, object, name[c], func[c]);
    gtk_widget_show (scalar);
    gtk_table_attach (table, scalar, 1, 2, c, c + 1,
		      (GtkAttachOptions) (GTK_FILL),
		      (GtkAttachOptions) (0), 0, 0);

    GtkWidget * label = gtk_label_new (name[c]);
    gtk_widget_show (label);
    gtk_table_attach (table, label, 0, 1, c, c + 1,
		      (GtkAttachOptions) (GTK_FILL),
		      (GtkAttachOptions) (0), 0, 0);
  }

  gl_vectors_update_interface (gls);
}

static gchar * gl_vectors_name (GfkGlClass * klass)
{
  static gchar name[] = "Vectors";
  return name;
}

static GtkWidget * gl_vectors_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "vectors-16x16.png");
}

static gchar * gl_vectors_properties (GfkGl * gl)
{
  gchar * f[3] = { NULL, NULL, NULL };
  FttComponent c;

  g_free (gl->props);
  for (c = 0; c < FTT_DIMENSION; c++)
    f[c] = GFS_GL_VECTORS (gl->gl)->expr[c]->str;
  gl->props = g_strjoin (",", f[0], f[1], f[2], NULL);
  return gl->props;
}

static void gl_vectors_class_init (GfkGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_vectors_read;
  klass->gl_class = gfs_gl_vectors_class ();
  klass->post_init = gl_vectors_post_init;
  klass->set_simulation = gl_vectors_set_simulation;
  klass->name = gl_vectors_name;
  klass->icon = gl_vectors_icon;
  klass->properties = gl_vectors_properties;
}

static void gl_vectors_init (GfkGlVectors * object)
{
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (object)->properties, "shading_label"), 
			    FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (object)->properties, "shading"), FALSE);
  object->vector = create_vector_params ();
  gfk_gl_set_sensitive (GFK_GL (object), GFK_GL_SCALAR (object)->scalar, FALSE);

  GtkWidget * m = gtk_menu_new ();
  GtkWidget * i = gtk_menu_item_new_with_label ("Default");
  g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_vector_color), object);
  gtk_menu_append (m, i);
  gtk_widget_show (i);

  i = gtk_menu_item_new_with_label ("Scalar");
  g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_vector_scalar), object);
  gtk_menu_append (m, i);
  gtk_widget_show (i);
  
  gtk_option_menu_set_menu (lookup_widget_params (object->vector, "color"), m);
  gtk_widget_show (m);

  gfk_gl_prepend_params (GFK_GL (object), object->vector, gtk_label_new ("Vector"));
  gtk_widget_show (object->vector);
}

GfkGlClass * gfk_gl_vectors_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_vectors_info = {
      "GfkGlVectors",
      sizeof (GfkGlVectors),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_vectors_class_init,
      (GtsObjectInitFunc) gl_vectors_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_scalar_class ()),
				  &gfk_gl_vectors_info);
  }

  return klass;
}

/* GfkGlStreamlines: Object */

static void gl_streamlines_read (GtsObject ** o, GtsFile * fp)
{
  GtkWidget * stream = GFK_GL_STREAMLINES (*o)->stream;
  GfsGlStreamlines * gl = GFS_GL_STREAMLINES (GFK_GL (*o)->gl);

  (* GTS_OBJECT_CLASS (gfk_gl_streamlines_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (!gl->show_cells)
    gtk_toggle_button_set_active (lookup_widget_params (stream, "showcells"), FALSE);
  gtk_spin_button_set_value (lookup_widget_params (stream, "dmin"), gl->dmin);
  if (gl->dmin > 0.)
    gtk_widget_set_sensitive (lookup_widget_params (stream, "evenly_spaced"), TRUE);
  gtk_spin_button_set_value (lookup_widget_params (stream, "radius"), gl->radius);
}

static gchar * gl_streamlines_name (GfkGlClass * klass)
{
  static gchar name[] = "Stream";
  return name;
}

static GtkWidget * gl_streamlines_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "stream-16x16.png");
}

static void snap_to_spacing (GfsGlStreamlines * gls)
{
  if (gls->dmin > 0.) {
    FttVector * picked = &GFS_GL2D (gls)->pickedpos;
    GtsPoint p;
    
    if (gfs_gl_streamlines_closest (gls, picked, &p) < G_MAXDOUBLE) {
      GtsVector n;
      gts_vector_init (n, &p, picked);
      gts_vector_normalize (n);
      picked->x = p.x + gls->dmin*n[0];
      picked->y = p.y + gls->dmin*n[1];
      picked->z = p.z + gls->dmin*n[2];
    }
  }
}

static gchar * gl_streamlines_pickinfo (GfkGl * gl, gboolean motion)
{
  GfkGlStreamlines * glk = GFK_GL_STREAMLINES (gl);
  GfsGlStreamlines * gls = GFS_GL_STREAMLINES (gl->gl);
  gboolean snapspacing = FALSE;

  if (gtk_toggle_button_get_active (lookup_widget_params (glk->stream, "snapcenters"))) {
    FttCell * cell = GFS_GL2D (gls)->picked;
    
    while (gl->gl->maxlevel >= 0 && gl->gl->maxlevel < ftt_cell_level (cell))
      cell = ftt_cell_parent (cell);
    ftt_cell_pos (cell, &GFS_GL2D (gls)->pickedpos);
  }
  else
    snapspacing = gtk_toggle_button_get_active (lookup_widget_params (glk->stream, "snapspacing"));

  if (motion) {
    if (gls->selected) {
      GfsGlStreamline * s = gls->selected->data;
      
      gfs_gl_streamlines_reset_selected (gls);
      if (snapspacing)
	snap_to_spacing (gls);
      s->c = GFS_GL2D (gls)->pickedpos;
      gfk_gl_expose (gl);
    }
  }
  else if (!glk->edit) {
    if (!gls->selected)
      gtk_widget_set_sensitive (lookup_widget_params (glk->stream, "delete"), TRUE);
    if (snapspacing)
      snap_to_spacing (gls);
    if (gfs_gl_streamlines_add (gls, GFS_GL2D (gls)->pickedpos))
      gfk_gl_expose (gl);
  }
  else if (gls->candidate) { 
    if (gls->candidate != gls->selected) {
      if (!gls->selected)
	gtk_widget_set_sensitive (lookup_widget_params (glk->stream, "delete"), TRUE);
      gls->selected = gls->candidate;
      gfk_gl_expose (gl);
    }
  }
  else if (gls->selected) {
    gtk_widget_set_sensitive (lookup_widget_params (glk->stream, "delete"), FALSE);
    gls->selected = NULL;
    gfk_gl_expose (gl);
  }
  return gl2D_pickinfo (gl, motion);
}

static void streamlines_changed (GtkEntry * entry, GfsGlStreamlines * gl)
{
  gfs_gl_streamlines_reset (gl);
}

static void gl_streamlines_post_init (GfkGl * object)
{
  FttComponent c;

  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_streamlines_class ())->parent_class)->post_init) 
    (object);

  for (c = 0; c < FTT_DIMENSION; c++) {
    static gchar * name[] = {"Vx", "Vy", "Vz"};
    GtkWidget * entry = lookup_widget_params (GFK_GL_VECTORS (object)->vector, name[c]);

    g_signal_connect (G_OBJECT (entry), "activate", G_CALLBACK (streamlines_changed), object->gl);
  }
}

static void gl_streamlines_class_init (GfkGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_streamlines_read;
  klass->gl_class = gfs_gl_streamlines_class ();
  klass->name = gl_streamlines_name;
  klass->icon = gl_streamlines_icon;
  klass->pickinfo = gl_streamlines_pickinfo;
  klass->post_init = gl_streamlines_post_init;
}

static void gl_streamlines_init (GfkGlVectors * gl)
{
  GfkGlStreamlines * gls = GFK_GL_STREAMLINES (gl);

  gtk_widget_set_sensitive (lookup_widget_params (gl->vector, "label11"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->vector, "scale"), FALSE);

  gls->stream = create_streamlines_params ();
  gtk_widget_set_sensitive (lookup_widget_params (gls->stream, "delete"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gls->stream, "evenly_spaced"), FALSE);
  if (FTT_DIMENSION < 3)
    gtk_widget_hide (lookup_widget_params (gls->stream, "params_3D"));
  gfk_gl_prepend_params (GFK_GL (gl), gls->stream, gtk_label_new ("Streamlines"));
  gtk_widget_show (gls->stream);
}

GfkGlClass * gfk_gl_streamlines_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_streamlines_info = {
      "GfkGlStreamlines",
      sizeof (GfkGlStreamlines),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_streamlines_class_init,
      (GtsObjectInitFunc) gl_streamlines_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_vectors_class ()),
				  &gfk_gl_streamlines_info);
  }

  return klass;
}

/* GfkGlEllipses: Object */

static void update_ellipse_scalar (GfkGlEllipses * gl)
{
  gfk_gl_set_sensitive (GFK_GL (gl), GFK_GL_SCALAR (gl)->scalar, TRUE);
  gfk_gl_update_properties (GFK_GL (gl));
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color_label"), 
			    FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color"), 
			    FALSE);
}

static void gl_ellipses_read (GtsObject ** o, GtsFile * fp)
{
  GtkWidget * ellipse = GFK_GL_ELLIPSES (*o)->ellipse;
  GfsGlEllipses * gl = GFS_GL_ELLIPSES (GFK_GL (*o)->gl);

  (* GTS_OBJECT_CLASS (gfk_gl_ellipses_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  gtk_spin_button_set_value (lookup_widget_params (ellipse, "scale"), gl->scale);

  if (gl->use_scalar) {
    gtk_option_menu_set_history (lookup_widget_params (ellipse, "color"), 1);
    update_ellipse_scalar (GFK_GL_ELLIPSES (*o));
  }
}

static void gl_ellipses_update_interface (GfkGl * g)
{
  GfkGlEllipses * gl = GFK_GL_ELLIPSES (g);
  GfsGlEllipses * gle = GFS_GL_ELLIPSES (GFK_GL (gl)->gl);
  GtkSpinButton * scale = lookup_widget_params (gl->ellipse, "scale");
  GtkAdjustment * ascale = gtk_spin_button_get_adjustment (scale);
  guint j;
  
  for (j = 0; j < 4; j++) {
    static gchar * name[] = {"Ax", "Ay", "Bx", "By"};
    GtkEntry * entry = lookup_widget_params (gl->ellipse, name[j]);
    if (!GFK_IS_EDITED (entry))
      gtk_entry_set_text (entry, gle->expr[j]->str);
  }

  ascale->upper = gle->max > 0. ? 10000.*gle->h/gle->max : 1.;
  ascale->step_increment = gle->max > 0. ? gle->h/gle->max/10. : 0.1;
  ascale->page_increment = 10.*ascale->step_increment;
  ascale->value = gle->scale;
  if (!GFK_IS_EDITED (scale))
    gtk_spin_button_configure (scale, ascale, 2.*ascale->step_increment,
			       gfk_decimal_digits (ascale->step_increment, 2));
  gfk_gl_update_properties (GFK_GL (gl));
  gfk_gl_expose (GFK_GL (gl));
}

static void set_ellipse_color (GtkWidget * color, GfkGlEllipses * gl)
{
  GFS_GL_ELLIPSES (GFK_GL (gl)->gl)->use_scalar = FALSE;
  gfk_gl_set_sensitive (GFK_GL (gl), GFK_GL_SCALAR (gl)->scalar, FALSE);
  gfk_gl_update_properties (GFK_GL (gl));
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color_label"), 
			    TRUE);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color"), TRUE);
  gfk_gl_expose (GFK_GL (gl));
}

static void set_ellipse_scalar (GtkWidget * scalar, GfkGlEllipses * gl)
{
  GFS_GL_ELLIPSES (GFK_GL (gl)->gl)->use_scalar = TRUE;
  update_ellipse_scalar (gl);
  gfk_gl_expose (GFK_GL (gl));
}

static void gl_ellipses_set_simulation (GfkGl * object, GfsSimulation * sim)
{
  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_ellipses_class ())->parent_class)->set_simulation)
    (object, sim);

  gl_ellipses_update_interface (object);
}

static GtsFile * set_ax (GfkGl * gl, const gchar * func)
{
  return gfs_gl_ellipses_set (GFS_GL_ELLIPSES (gl->gl), 0, func);
}

static GtsFile * set_ay (GfkGl * gl, const gchar * func)
{
  return gfs_gl_ellipses_set (GFS_GL_ELLIPSES (gl->gl), 1, func);
}

static GtsFile * set_bx (GfkGl * gl, const gchar * func)
{
  return gfs_gl_ellipses_set (GFS_GL_ELLIPSES (gl->gl), 2, func);
}

static GtsFile * set_by (GfkGl * gl, const gchar * func)
{
  return gfs_gl_ellipses_set (GFS_GL_ELLIPSES (gl->gl), 3, func);
}

static void gl_ellipses_post_init (GfkGl * object)
{
  GfkGlEllipses * gls = GFK_GL_ELLIPSES (object);
  guint i;

  if (GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_ellipses_class ())->parent_class)->post_init)
    (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_ellipses_class ())->parent_class)->post_init)
      (object);

  GtkTable * table = GTK_TABLE (lookup_widget_params (gls->ellipse, "table12"));
  for (i = 0; i < 4; i++) {
    static gchar * name[] = {"Ax", "Ay", "Bx", "By"};
    static GfkFunctionSet func[] = {set_ax, set_ay, set_bx, set_by};

    GtkWidget * scalar = gfk_function (gls->ellipse, object, name[i], func[i]);
    gtk_widget_show (scalar);
    gtk_table_attach (table, scalar, 1, 2, i, i + 1,
		      (GtkAttachOptions) (GTK_FILL),
		      (GtkAttachOptions) (0), 0, 0);
  }

  gl_ellipses_update_interface (object);
}

static gchar * gl_ellipses_name (GfkGlClass * klass)
{
  static gchar name[] = "Ellipses";
  return name;
}

static GtkWidget * gl_ellipses_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "ellipses-16x16.png");
}

static gchar * gl_ellipses_properties (GfkGl * gl)
{
  GfsGlEllipses * gle = GFS_GL_ELLIPSES (gl->gl);
  
  g_free (gl->props);
  gl->props = g_strjoin (",", 
			 gle->expr[0]->str, gle->expr[1]->str, 
			 gle->expr[2]->str, gle->expr[3]->str, NULL);
  return gl->props;
}

static void gl_ellipses_class_init (GfkGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_ellipses_read;
  klass->gl_class = gfs_gl_ellipses_class ();
  klass->post_init = gl_ellipses_post_init;
  klass->set_simulation = gl_ellipses_set_simulation;
  klass->update_interface = gl_ellipses_update_interface;
  klass->name = gl_ellipses_name;
  klass->icon = gl_ellipses_icon;
  klass->properties = gl_ellipses_properties;
}

static void gl_ellipses_init (GfkGlEllipses * object)
{
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (object)->properties, "shading_label"), 
			    FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (object)->properties, "shading"), FALSE);

  object->ellipse = create_ellipse_params ();

  gfk_gl_set_sensitive (GFK_GL (object), GFK_GL_SCALAR (object)->scalar, FALSE);
  GtkWidget * m = gtk_menu_new ();
  GtkWidget * i = gtk_menu_item_new_with_label ("Default");
  g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_ellipse_color), object);
  gtk_menu_append (m, i);
  gtk_widget_show (i);

  i = gtk_menu_item_new_with_label ("Scalar");
  g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_ellipse_scalar), object);
  gtk_menu_append (m, i);
  gtk_widget_show (i);
  
  gtk_option_menu_set_menu (lookup_widget_params (object->ellipse, "color"), m);
  gtk_widget_show (m);

  gfk_gl_prepend_params (GFK_GL (object), object->ellipse, gtk_label_new ("Ellipse"));
  gtk_widget_show (object->ellipse);
}

GfkGlClass * gfk_gl_ellipses_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_ellipses_info = {
      "GfkGlEllipses",
      sizeof (GfkGlEllipses),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_ellipses_class_init,
      (GtsObjectInitFunc) gl_ellipses_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_scalar_class ()),
				  &gfk_gl_ellipses_info);
  }

  return klass;
}

/* GfkGlLocation: Object */

static void gl_location_read (GtsObject ** o, GtsFile * fp)
{
  GtkWidget * location = GFK_GL_LOCATION (*o)->location;
  GfsGlLocation * gl = GFS_GL_LOCATION (GFK_GL (*o)->gl);

  (* GTS_OBJECT_CLASS (gfk_gl_location_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  
  gtk_spin_button_set_value (lookup_widget_params (location, "size"), gl->size);
#if HAVE_FTGL
  if (gl->label) {
    gl->label = FALSE;
    gtk_toggle_button_set_active (lookup_widget_params (location, "location_label_check"), TRUE);
  }
#endif /* HAVE_FTGL */
}

static void gl_location_init (GfkGl * gl)
{
  GtkWidget * location = create_location_params ();

  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "maxlevel_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "maxlevel"), FALSE);
#if HAVE_FTGL
  gtk_widget_show (gl->font);
#else /* not HAVE_FTGL */
  gtk_widget_set_sensitive (lookup_widget_params (location, "location_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (location, "location_label_check"), FALSE);
#endif /* not HAVE_FTGL */
  gfk_gl_prepend_params (gl, location, gtk_label_new ("Location"));
  gtk_spin_button_set_value (lookup_widget_params (location, "size"), 
			     GFS_GL_LOCATION (gl->gl)->size);
  gtk_widget_show (location);
  GFK_GL_LOCATION (gl)->location = location;
}

static gchar * gl_location_name (GfkGlClass * klass)
{
  static gchar name[] = "Location";
  return name;
}

static GtkWidget * gl_location_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "location-16x16.png");
}

static void gl_location_class_init (GfkGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_location_read;
  klass->gl_class = gfs_gl_location_class ();
  klass->name = gl_location_name;
  klass->icon = gl_location_icon;
}

GfkGlClass * gfk_gl_location_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_location_info = {
      "GfkGlLocation",
      sizeof (GfkGlLocation),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_location_class_init,
      (GtsObjectInitFunc) gl_location_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_class ()), &gfk_gl_location_info);
  }

  return klass;
}

/* GfkGlHeight: Object */

static gchar * gl_height_name (GfkGlClass * klass)
{
  static gchar name[] = "Height";
  return name;
}

static GtkWidget * gl_height_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "height-16x16.png");
}

static void gl_height_class_init (GfkGlClass * klass)
{
  klass->gl_class = gfs_gl_height_class ();
  klass->name = gl_height_name;
  klass->icon = gl_height_icon;
}

GfkGlClass * gfk_gl_height_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfkGlHeight",
      sizeof (GfkGlLocation),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_height_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_location_class ()), &info);
  }

  return klass;
}

/* GfkGlLocate: Object */

static void gl_locate_read (GtsObject ** o, GtsFile * fp)
{
  GtkWidget * locate = GFK_GL_LOCATE (*o)->locate;
  GfsGlLocate * gl = GFS_GL_LOCATE (GFK_GL (*o)->gl);

  (* GTS_OBJECT_CLASS (gfk_gl_locate_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  FttComponent c;
  for (c = 0; c < 3; c++) {
    static gchar * name[] = {"locate_x_entry", "locate_y_entry", "locate_z_entry"};
    gchar * s = g_strdup_printf ("%g", (&gl->p.x)[c]);
    gtk_entry_set_text (lookup_widget_params (locate, name[c]), s);
    g_free (s);
  }
}

static gchar * gl_locate_name (GfkGlClass * klass)
{
  static gchar name[] = "Locate";
  return name;
}

static GtkWidget * gl_locate_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "locate-16x16.png");
}

static void gl_locate_update_interface (GfkGl * g)
{
  GtkWidget * locate = GFK_GL_LOCATE (g)->locate;
  GfsGlLocate * gl = GFS_GL_LOCATE (g->gl);

  FttComponent c;
  for (c = 0; c < 3; c++) {
    static gchar * name[] = {"locate_x_entry", "locate_y_entry", "locate_z_entry"};
    GtkEntry * entry = lookup_widget_params (locate, name[c]);
    gchar * s = g_strdup_printf ("%g", (&gl->p.x)[c]);
    gtk_entry_set_text (entry, s);
    g_free (s);
  }
  gfk_gl_update_properties (g);
  gfk_gl_expose (g);
}

static void gl_locate_set_simulation (GfkGl * object, GfsSimulation * sim)
{
  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_locate_class ())->parent_class)->set_simulation)
    (object, sim);

  gl_locate_update_interface (object);
}

static void gl_locate_class_init (GfkGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_locate_read;
  klass->gl_class = gfs_gl_locate_class ();
  klass->name = gl_locate_name;
  klass->icon = gl_locate_icon;
  klass->update_interface = gl_locate_update_interface;
  klass->set_simulation = gl_locate_set_simulation;
}

static void gl_locate_init (GfkGl * gl)
{
  GtkWidget * locate = create_locate_params ();

  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading"), FALSE);
  gfk_gl_prepend_params (gl, locate, gtk_label_new ("Locate"));
  gtk_widget_show (locate);

  GFK_GL_LOCATE (gl)->locate = locate;
}

GfkGlClass * gfk_gl_locate_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_locate_info = {
      "GfkGlLocate",
      sizeof (GfkGl),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_locate_class_init,
      (GtsObjectInitFunc) gl_locate_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_class ()), &gfk_gl_locate_info);
  }

  return klass;
}

/* GfkGlPipes: Object */

static void gl_pipes_init (GfkGl * gl)
{
  gtk_widget_hide (lookup_widget_params (gl->properties, "maxlevel_label"));
  gtk_widget_hide (lookup_widget_params (gl->properties, "maxlevel"));
  gtk_widget_hide (lookup_widget_params (gl->properties, "finest"));
  gtk_widget_hide (lookup_widget_params (gl->properties, "shading_label"));
  gtk_widget_hide (lookup_widget_params (gl->properties, "shading"));
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "label117"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "linewidth"), FALSE);
#if HAVE_FTGL
  gtk_widget_show (gl->font);
#endif
}

static gchar * gl_pipes_name (GfkGlClass * klass)
{
  static gchar name[] = "Pipes";
  return name;
}

static GtkWidget * gl_pipes_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "pipes-16x16.png");
}

static void gl_pipes_class_init (GfkGlClass * klass)
{
  klass->gl_class = gfs_gl_pipes_class ();
  klass->name = gl_pipes_name;
  klass->icon = gl_pipes_icon;
}

GfkGlClass * gfk_gl_pipes_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfkGlPipes",
      sizeof (GfkGl),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_pipes_class_init,
      (GtsObjectInitFunc) gl_pipes_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_class ()), &info);
  }

  return klass;
}

/* GfkGlClipPlane: Object */

static gchar * gl_clip_plane_name (GfkGlClass * klass)
{
  static gchar name[] = "Clipping";
  return name;
}

static GtkWidget * gl_clip_plane_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "clipping-16x16.png");
}

static void gl_clip_plane_class_init (GfkGlClass * klass)
{
  klass->gl_class = gfs_gl_clip_plane_class ();
  klass->name = gl_clip_plane_name;
  klass->icon = gl_clip_plane_icon;
  klass->pickinfo = NULL;
}

static void gl_clip_plane_init (GfkGl * gl)
{
  gtk_widget_show (GFK_GL2D (gl)->params);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "gl_params"), FALSE);
}

GfkGlClass * gfk_gl_clip_plane_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_clip_plane_info = {
      "GfkGlClipPlane",
      sizeof (GfkGl2D),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_clip_plane_class_init,
      (GtsObjectInitFunc) gl_clip_plane_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl2D_class ()), &gfk_gl_clip_plane_info);
  }

  return klass;
}

/* GfkGlInfo: Object */

typedef struct {
  guint cells, leaves;
  GArray * level;
} CellCount;

static void count_cells (FttCell * cell, CellCount * c)
{
  guint i = ftt_cell_level (cell), n = 0;
  if (i >= c->level->len)
    g_array_append_val (c->level, n);
  c->cells++;
  if (FTT_CELL_IS_LEAF (cell)) {
    c->leaves++;
    g_array_index (c->level, guint, i)++;
  }
}

static void gl_info_set_simulation (GfkGl * object, GfsSimulation * sim)
{
  CellCount count = {0, 0, NULL};

  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_info_class ())->parent_class)->set_simulation) 
    (object, sim);
  
  count.level = g_array_new (FALSE, FALSE, sizeof (guint));
  gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) count_cells, &count);
  GtkLabel * label = lookup_widget_params (GFK_GL_INFO (object)->info, "number_of_levels");
  gchar * s = g_strdup_printf ("%d", count.level->len);
  gtk_label_set_text (label, s);
  g_free (s);

  label = lookup_widget_params (GFK_GL_INFO (object)->info, "number_of_cells");
  s = g_strdup_printf ("%d", count.cells);
  gtk_label_set_text (label, s);
  g_free (s);

  label = lookup_widget_params (GFK_GL_INFO (object)->info, "number_of_leaf_cells");
  s = g_strdup_printf ("%d", count.leaves);
  gtk_label_set_text (label, s);
  g_free (s);

  GtkWidget * table = lookup_widget_params (GFK_GL_INFO (object)->info, "leaf_cells_per_level");
  gtk_container_foreach (GTK_CONTAINER (table), (GtkCallback) gtk_widget_destroy, NULL);

  guint i, l = 0;
  for (i = 0; i < count.level->len; i++)
    if (g_array_index (count.level, guint, i)) {
      gchar * s = g_strdup_printf ("Level %d", i);
      GtkWidget * label = gtk_label_new (s);
      gtk_widget_show (label);
      g_free (s);
      gtk_table_attach (GTK_TABLE (table), label, 0, 1, l, l + 1,
			(GtkAttachOptions) (GTK_FILL),
			(GtkAttachOptions) (0), 16, 0);
      gtk_misc_set_alignment (GTK_MISC (label), 0, 0.5);
      
      s = g_strdup_printf ("%d", g_array_index (count.level, guint, i));
      label = gtk_label_new (s);
      gtk_widget_show (label);
      g_free (s);
      gtk_table_attach (GTK_TABLE (table), label, 1, 2, l, l + 1,
			(GtkAttachOptions) (GTK_FILL),
			(GtkAttachOptions) (0), 0, 0);
      gtk_label_set_justify (GTK_LABEL (label), GTK_JUSTIFY_RIGHT);
      gtk_misc_set_alignment (GTK_MISC (label), 1, 0.5);
      
      l++;
    }

  g_array_free (count.level, TRUE);
}

static gchar * gl_info_name (GfkGlClass * klass)
{
  static gchar name[] = "Info";
  return name;
}

static GtkWidget * gl_info_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "info-16x16.png");
}

static void gl_info_class_init (GfkGlClass * klass)
{
  klass->set_simulation = gl_info_set_simulation;
  klass->name = gl_info_name;
  klass->icon = gl_info_icon;
}

static void gl_info_init (GfkGl * gl)
{
  GFK_GL_INFO (gl)->info = create_info_params ();
  gfk_gl_set_sensitive (gl, gl->properties, FALSE);
  gfk_gl_prepend_params (gl, GFK_GL_INFO (gl)->info, gtk_label_new ("Info"));
  gtk_widget_show (GFK_GL_INFO (gl)->info);
}

GfkGlClass * gfk_gl_info_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_info_info = {
      "GfkGlInfo",
      sizeof (GfkGlInfo),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_info_class_init,
      (GtsObjectInitFunc) gl_info_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_class ()),
				  &gfk_gl_info_info);
  }

  return klass;
}

/* GfkGlView: Object */

static void write_ppm_pixbuf (GdkPixbuf * pixbuf, FILE * fp)
{
  guint width, height = gdk_pixbuf_get_height (pixbuf);
  guchar * pixels = gdk_pixbuf_get_pixels (pixbuf);
  guint j, rowstride = gdk_pixbuf_get_rowstride (pixbuf);

  g_assert (gdk_pixbuf_get_colorspace (pixbuf) == GDK_COLORSPACE_RGB);
  g_assert (gdk_pixbuf_get_bits_per_sample (pixbuf) == 8);
  g_assert (!gdk_pixbuf_get_has_alpha (pixbuf));
  g_assert (gdk_pixbuf_get_n_channels (pixbuf) == 3);

  width = rowstride/3;
  if (width % 2)
    width -= 1;
  if (height % 2)
    height -= 1;
  fprintf (fp, "P6 %d %d 255\n", width, height);
  for (j = 0; j < height; j++, pixels += rowstride)
    fwrite (pixels, 3*width, sizeof (guchar), fp);  
}

static void write_ppm (guint width, guint height, FILE * fp)
{
  gchar * p, * p1;
  guint j;

  p = g_malloc0 (3*width*height*sizeof (gchar));
  glPixelStorei (GL_PACK_ALIGNMENT, 1);
  glReadPixels (0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, p);
  fprintf (fp, "P6 %d %d 255\n", width, height);
  p1 = (gchar *)(p + sizeof (gchar)*(3*width*(height - 1)));
  for (j = height; j--; p1 = (gchar *)(p1 - sizeof (gchar)*3*width))
    fwrite (p1, 3*width, sizeof (gchar), fp);
  g_free (p);
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

void gfs_gl2ps (GfsGl2PSParams * p, 
		FILE * fp, 
		const gchar * fname, 
		GtkWidget * view)
{
  GtkWidget * glarea;
  GfsGlViewParams * viewp;

  g_return_if_fail (p != NULL);
  g_return_if_fail (fp != NULL);
  g_return_if_fail (fname != NULL);
  g_return_if_fail (view != NULL);

  glarea = g_object_get_data (G_OBJECT (view), "glarea");
  viewp = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");

  switch (p->format) {
  case GFSGL_PPM_OFFSCREEN: {
    gboolean written = FALSE;
    GdkGLConfig * glconfig = gdk_gl_config_new_by_mode (GDK_GL_MODE_RGB    |
							GDK_GL_MODE_DEPTH  |
							GDK_GL_MODE_SINGLE);
    if (!glconfig)
      g_warning ("Cannot get OpenGL config\n");
    else {
      guint width  = p->width ?  p->width :  glarea->allocation.width;
      guint height = p->height ? p->height : glarea->allocation.height;
      GdkPixmap * pixmap = gdk_pixmap_new (glarea->window, width, height, -1);
      GdkGLDrawable * gldrawable = GDK_GL_DRAWABLE (gdk_pixmap_set_gl_capability (pixmap,
										  glconfig,
										  NULL));
      GdkGLContext * glcontext = gdk_gl_context_new (gldrawable,
						     NULL,
						     FALSE,
						     GDK_GL_RGBA_TYPE);
    
      if (!glcontext)
	g_warning ("Cannot create the OpenGL rendering context\n");
      else {
	GfsGlViewParams * info = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");
	GLfloat m[4][4];
	GdkPixbuf * pixbuf;
	gdouble max;

	gdk_gl_drawable_gl_begin (gldrawable, glcontext);

	glViewport (0, 0, width, height);
	gfs_gl_init_gl ();
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	GtkTreeModel * list = 
	  gtk_tree_view_get_model (GTK_TREE_VIEW (lookup_widget (view, "gl_list")));
	GList * symmetries = get_symmetries (list);
	max = gfs_gl_domain_extent (g_object_get_data (G_OBJECT (glarea), "sim"), symmetries);
	g_list_free (symmetries);
	gluPerspective (info->fov, width/(float)height, 1., 1. + 2.*max);
	glMatrixMode (GL_MODELVIEW);
      
	glClearColor (info->bg.r, info->bg.g, info->bg.b, 1);
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity ();
	glTranslatef (info->tx, info->ty, - (1. + max));
	gfs_gl_build_rotmatrix (m, info->quat);
	glMultMatrixf (&m[0][0]);
	glScalef (info->sx, info->sy, info->sz);

	gfk_gl_view_draw (view, GFSGL_PPM_OFFSCREEN);

	glFlush ();
	gdk_gl_drawable_gl_end (gldrawable);

	pixbuf = gdk_pixbuf_get_from_drawable (NULL, GDK_DRAWABLE (pixmap), NULL, 
					       0, 0, 0, 0, -1, -1);

	if (pixbuf == NULL)
	  g_warning ("Cannot get pixbuf\n");
	else {
	  write_ppm_pixbuf (pixbuf, fp);
	  written = TRUE;
	  g_object_unref (G_OBJECT (pixbuf));
	}
	gdk_gl_context_destroy (glcontext);
	info->do_init = TRUE; /* this is necessary to reinitialise textures for fonts */
      }
      gdk_pixmap_unset_gl_capability (pixmap);
      g_object_unref (G_OBJECT (pixmap));
    }
    if (written)
      break;
  }
  case GFSGL_PPM_SCREEN:
    write_ppm (glarea->allocation.width, glarea->allocation.height, fp);
    break;
  case GFSGL_GNUPLOT: case GFSGL_OBJ: case GFSGL_KML: {
    guint buffsize = 0;
    gboolean done = FALSE;
    gfloat base_res = viewp->base_res;
    viewp->base_res = 0.;
    while (!done) {
      GfsGlFeedback * f;

      buffsize += 2048*2048;
      f = gfs_gl_feedback_begin (buffsize);
      gfk_gl_view_draw (view, GFSGL_SCREEN);
      done = gfs_gl_feedback_end (f, g_object_get_data (G_OBJECT (glarea), "sim"), fp, p->format);
    }
    viewp->base_res = base_res;
    break;
  }
  default: {
    GLint buffsize = 0, state = GL2PS_OVERFLOW;

    while (state == GL2PS_OVERFLOW) {
      buffsize += 2048*2048;
      gl2psBeginPage ("", "GfsView",
		      NULL,
		      p->format, p->sort, p->options, 
		      GL_RGBA, 0, NULL, 
		      0, 0, 0,
		      buffsize, fp, fname);
      GtkWidget * glarea = g_object_get_data (G_OBJECT (view), "glarea");
      GfsGlViewParams * v = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");
      v->lw = p->lw;
      gfk_gl_view_draw (view, p->format);
      state = gl2psEndPage();
    }
  }
  }
}

static void visible_toggled (GtkCellRendererToggle * cellrenderertoggle,
			     gchar * path_string,
			     GtkTreeModel * model)
{
  GtkTreeIter iter;
  gboolean visible;
  GfkGl * gl;

  g_assert (gtk_tree_model_get_iter_from_string (model, &iter, path_string));
  gtk_tree_model_get (model, &iter, GL_COLUMN, &gl, VISIBLE_COLUMN, &visible, -1);
  gtk_list_store_set (GTK_LIST_STORE (model), &iter, VISIBLE_COLUMN, !visible, -1);
  gfk_gl_expose (gl);
}

static gboolean tree_selection_func (GtkTreeSelection * select,
				     GtkTreeModel * model,
				     GtkTreePath * path,
				     gboolean path_currently_selected,
				     gpointer data)
{
  if (path_currently_selected)
    return TRUE;
  else {
    GtkTreeIter iter;
    gboolean visible;

    g_assert (gtk_tree_model_get_iter (model, &iter, path));
    gtk_tree_model_get (model, &iter, VISIBLE_COLUMN, &visible, -1);
    if (visible && gtk_tree_selection_get_selected (select, &model, &iter))
      gtk_list_store_set (GTK_LIST_STORE (model), &iter, SELECTED_COLUMN, FALSE, -1);
    return visible;
  }
}

static void tree_selection_changed_cb (GtkTreeSelection * selection, GObject * list)
{
  GtkTreeIter iter;
  GtkTreeModel * model;
  
  if (gtk_tree_selection_get_selected (selection, &model, &iter)) {
    GfkGl * gl, * former = g_object_get_data (list, "former");

    gtk_tree_model_get (model, &iter, GL_COLUMN, &gl, -1);
    gtk_list_store_set (GTK_LIST_STORE (model), &iter, SELECTED_COLUMN, TRUE, -1);
    if (former && former != gl) {
      gint x, y;
      gtk_window_get_position (GTK_WINDOW (former->params), &x, &y);
      gtk_window_move (GTK_WINDOW (gl->params), x, y);
      gtk_widget_show (gl->params);
      gtk_widget_hide (former->params);
      g_object_set_data (list, "former", gl);
    }
    gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (list), "properties1"), TRUE);
    gtk_widget_set_sensitive (lookup_widget (GTK_WIDGET (list), "delete3"), TRUE);
  }
}

static gchar * format_glob (GSList * list)
{
  gchar * s = g_strdup ("");

  while (list) {
    GfsFormat * f = list->data;
    gchar * s2 = NULL;

    switch (f->t) {
    case GFS_NONE_FORMAT:
      s2 = g_strconcat (s, f->s, NULL);
      break;
    case GFS_ITER_FORMAT: case GFS_TIME_FORMAT: 
      s2 = g_strconcat (s, "*", NULL);
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

static void free_glob (glob_t * files)
{
  guint i;

  for (i = 0; i < files->gl_pathc; i++)
    g_free (files->gl_pathv[i]);
  g_free (files->gl_pathv);
  g_free (files);
}

static void read_formats (const gchar * fname,
			  GtkWidget * view)
{
  GtkWindow * parent = GTK_WINDOW (view);
  gboolean dynamic = FALSE, parallel = FALSE;
  GSList * formats = gfs_format_new (fname, NULL, &dynamic, &parallel);

  if (formats == NULL || parallel) {
    GtkWidget * msg = gtk_message_dialog_new (parent,
					      GTK_DIALOG_DESTROY_WITH_PARENT,
					      GTK_MESSAGE_WARNING,
					      GTK_BUTTONS_CLOSE,
					      "Cannot open file `%s`",
					      fname);
    gtk_dialog_run (GTK_DIALOG (msg));
    gtk_widget_destroy (msg);
  }
  else {
    gchar * g = format_glob (formats);
    glob_t * files;

    files = g_malloc (sizeof (glob_t));
    glob (g, GLOB_NOSORT, NULL, files);
    if (!files->gl_pathc) {
      GtkWidget * msg = gtk_message_dialog_new (parent,
						GTK_DIALOG_DESTROY_WITH_PARENT,
						GTK_MESSAGE_WARNING,
						GTK_BUTTONS_CLOSE,
						"Pattern `%s` did not match any files",
						fname);
      gtk_dialog_run (GTK_DIALOG (msg));
      gtk_widget_destroy (msg);
    }
    else {
      GtkWidget * play;
      guint i, j;

      for (i = 0; i < files->gl_pathc - 1; i++)
	for (j = 0; j < files->gl_pathc - 1 - i; j++)
	  if (gfs_format_time_value (formats, files->gl_pathv[j + 1]) < 
	      gfs_format_time_value (formats, files->gl_pathv[j])) {
	    gchar * tmp = files->gl_pathv[j];
	    files->gl_pathv[j] = files->gl_pathv[j + 1];
	    files->gl_pathv[j + 1] = tmp;
	  }

      i = 0;
      while (i < files->gl_pathc && 
	     gfs_format_time_value (formats, files->gl_pathv[i]) < G_MAXDOUBLE)
	i++;
      j = i;
      while (j < files->gl_pathc) g_free (files->gl_pathv[j++]);
      files->gl_pathc = i;
      if (files->gl_pathc == 0) {
	GtkWidget * msg = gtk_message_dialog_new (parent,
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_MESSAGE_WARNING,
						  GTK_BUTTONS_CLOSE,
						  "Pattern `%s` did not match any files",
						  fname);
	gtk_dialog_run (GTK_DIALOG (msg));
	gtk_widget_destroy (msg);
	free_glob (files);
	gfs_format_destroy (formats);
	return;
      }
      
      files->gl_offs = 0;

      play = create_play ();
      g_object_set_data_full (G_OBJECT (play), "glob_t", files, (GDestroyNotify) free_glob);
      g_object_set_data (G_OBJECT (play), "parent", parent);
      gtk_window_set_transient_for (GTK_WINDOW (play), parent);
      gtk_window_set_position (GTK_WINDOW (play), GTK_WIN_POS_CENTER_ON_PARENT);
      gtk_widget_show (play);

      gfk_gl_simulation_read (files->gl_pathv[0], view, TRUE);
    }
    g_free (g);
  }
  gfs_format_destroy (formats);
}

static GfsSimulation * gl_simulation_read (GtsFile * fp,
					   const gchar * fname,
					   GtkWindow * parent)
{
  GfsSimulation * sim;

  if ((sim = gfs_simulation_read (fp)) == NULL) {
    gchar * basename = g_path_get_basename (fname);
    GtkWidget * msg = gtk_message_dialog_new (parent,
					      GTK_DIALOG_DESTROY_WITH_PARENT,
					      GTK_MESSAGE_WARNING,
					      GTK_BUTTONS_CLOSE,
					      "File `%s' is not a valid Gerris file\n"
					      "%s:%d:%d: %s\n",
					      basename, basename,
					      fp->line, fp->pos, fp->error);
    g_free (basename);
    gtk_dialog_run (GTK_DIALOG (msg));
    gtk_widget_destroy (msg);
    return NULL;
  }
  else
    gfs_simulation_init (sim);
  return sim;
}

gboolean gfk_gl_view_read_parameters (GtkWidget * view, GtsFile * fp, gboolean discard)
{
  GfkGl * gl = NULL;
  GtkWidget * pref, * glarea, * list;
  GfsGlViewParams * params;
  GfsSimulation * sim;

  g_return_val_if_fail (view != NULL, FALSE);
  g_return_val_if_fail (fp != NULL, FALSE);

  pref = lookup_widget (view, "preferences");
  glarea = lookup_widget (view, "glarea");
  params = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");
  list = lookup_widget (view, "gl_list");
  sim = g_object_get_data (G_OBJECT (glarea), "sim");

  while (fp->type == GTS_STRING) {
    if (!strcmp (fp->token->str, "View")) {
      GfsGlViewParams p = *params;

      gfs_gl_view_params_read (&p, fp);
      if (!discard) {
	GtkWidget * bgcolor = lookup_widget (pref, "bgcolor");
	GtkWidget * colorsel = g_object_get_data (G_OBJECT (bgcolor), "colorsel");
	GdkColor c;

	*params = p;
	c.red = params->bg.r*65535; c.green = params->bg.g*65535; c.blue = params->bg.b*65535;
	gtk_widget_modify_bg (bgcolor, GTK_STATE_NORMAL, &c);
	gtk_widget_modify_bg (bgcolor, GTK_STATE_PRELIGHT, &c);
	gtk_widget_modify_bg (bgcolor, GTK_STATE_ACTIVE, &c);
	gtk_color_selection_set_current_color (GTK_COLOR_SELECTION (lookup_widget 
								    (colorsel, "colorselection2")), 
					       &c);
	gtk_spin_button_set_value (lookup_widget_params (pref, "resolution"), params->base_res);
	gtk_spin_button_set_value (lookup_widget_params (pref, "lc"), params->lc);
	gtk_spin_button_set_value (lookup_widget_params (pref, "reactivity"), params->reactivity);
	gtk_spin_button_set_value (lookup_widget_params (pref, "sx"), params->sx);
	gtk_spin_button_set_value (lookup_widget_params (pref, "sy"), params->sy);
	gtk_spin_button_set_value (lookup_widget_params (pref, "sz"), params->sz);
      }
    }
    else {
      if (discard) {
	GfsGl * g;

	if ((g = gfs_gl_new_from_file (fp)))
	  gts_object_destroy (GTS_OBJECT (g));
	else
	  break;
      }
      else {
	GtsObjectClass * klass = gts_object_class_from_name (fp->token->str);
	GtsObject * o;
      
	if (klass == NULL) {
	  gchar * ename = g_strconcat ("GfkGl", fp->token->str, NULL);
	  klass = gts_object_class_from_name (ename);
	  g_free (ename);
	}
	if (klass == NULL || !gts_object_class_is_from_class (klass, gfk_gl_class ()))
	  break;
	gl = gfk_gl_new (GFK_GL_CLASS (klass), glarea, list);
	o = GTS_OBJECT (gl);
	(* klass->read) (&o, fp);
	if (fp->type == GTS_ERROR) {
	  gts_object_destroy (o);
	  return FALSE;
	}
	if (sim)
	  gfk_gl_set_simulation (gl, sim);
	gl_add_gl (glarea, list, gl);
      }
    }
    if (fp->type == GTS_ERROR)
      return FALSE;
    gts_file_next_token (fp);
  }
  if (gl != NULL)
    gfk_gl_expose (gl);
  return TRUE;
}

GfsSimulation * gfk_gl_simulation_read (const gchar * fname,
					GtkWidget * view,
					gboolean set)
{
  GtsFile * fp;
  GfsSimulation * sim = NULL;
  FILE * fptr;

  g_return_val_if_fail (fname != NULL, NULL);
  g_return_val_if_fail (view != NULL, NULL);

  if (g_file_test (fname, G_FILE_TEST_IS_REGULAR))
    fptr = gfs_gl_popen (fname);
  else {
    read_formats (fname, view);
    return NULL;
  }
  if (fptr == NULL)
    return NULL;
  fp = gts_file_new (fptr);
  while (fp->type == GTS_STRING || fp->type == GTS_INT) {
    if (!gfk_gl_view_read_parameters (view, fp, FALSE)) {
      gchar * basename = g_path_get_basename (fname);
      GtkWidget * msg = gtk_message_dialog_new (GTK_WINDOW (view),
						GTK_DIALOG_DESTROY_WITH_PARENT,
						GTK_MESSAGE_WARNING,
						GTK_BUTTONS_CLOSE,
						"File `%s' is not a valid GfsView file\n"
						"%s:%d:%d: %s\n",
						basename, basename,
						fp->line, fp->pos, fp->error);
      g_free (basename);
      gtk_dialog_run (GTK_DIALOG (msg));
      gtk_widget_destroy (msg);
      gts_file_destroy (fp);
      pclose (fptr);
      return sim;
    }
    if (fp->type == GTS_INT || !strcmp (fp->token->str, "GModule")) {
      GfsSimulation * sim1 = gl_simulation_read (fp, fname, GTK_WINDOW (view));
      if (sim1) {
	if (sim && !set)
	  gts_object_destroy (GTS_OBJECT (sim));
	sim = sim1;
	if (set)
	  gfk_gl_view_set_simulation (view, sim, fname);
      }
    }
  }
  gts_file_destroy (fp);
  pclose (fptr);
  return sim;
}

static gint filew_ok (GtkWidget * widget, GtkWidget * filew)
{
  GtkWidget * view = g_object_get_data (G_OBJECT (filew), "view");
  const gchar * fname = gtk_file_selection_get_filename (GTK_FILE_SELECTION (filew));

  gtk_widget_hide (filew);
  gfk_gl_simulation_read (fname, view, TRUE);
  return TRUE;
}

static void substract (FttCell * from, FttCell * cell, gpointer * data)
{
  GfsDomain * domain = data[0];
  GSList * i = data[1], * j = domain->variables;

  while (j) {
    GfsVariable * v = j->data;
    if (v->name) {
      gint index = GPOINTER_TO_INT (i->data);
      if (index >= 0)
	GFS_VALUE (from, v) -= GFS_VALUEI (cell, index);
      i = i->next;
    }
    j = j->next;
  }
}

static void outside (FttCell * cell, gpointer * data)
{
  GfsDomain * domain = data[0];
  GSList * i = data[1], * j = domain->variables;
  
  while (j) {
    GfsVariable * v = j->data;
    if (v->name) {
      gint index = GPOINTER_TO_INT (i->data);
      if (index >= 0)
	GFS_VALUE (cell, v) = 0.;
      i = i->next;
    }
    j = j->next;
  }
}

static gint filews_ok (GtkWidget * widget, GtkWidget * filew)
{
  GtkWidget * view = g_object_get_data (G_OBJECT (filew), "view");
  const gchar * fname = gtk_file_selection_get_filename (GTK_FILE_SELECTION (filew));
  GfsSimulation * sim;
  gpointer data[2];

  gtk_widget_hide (filew);
  if ((sim = gfk_gl_simulation_read (fname, view, FALSE)) != NULL) {
    GtkWidget * glarea = lookup_widget (view, "glarea");
    GfsDomain * domain = g_object_get_data (G_OBJECT (glarea), "sim");
    GSList * i = domain->variables;
    gchar * missing = g_strdup (""), * title, * basename, * s;
    GSList * indices = NULL;

    while (i) {
      GfsVariable * v = i->data;
      if (v->name) {
	GfsVariable * v1 = gfs_variable_from_name (GFS_DOMAIN (sim)->variables, v->name);
	
	if (v1 != NULL)
	  indices = g_slist_append (indices, GINT_TO_POINTER ((gint) v1->i));
	else {
	  gchar * s = missing;
	  missing = s[0] == '\0' ? g_strdup (v->name) : g_strjoin (",", missing, v->name, NULL);
	  g_free (s);
	  indices = g_slist_append (indices, GINT_TO_POINTER (-1));
	}
      }
      i = i->next;
    }

    if (missing[0] != '\0') {
      GtkWidget * msg = gtk_message_dialog_new (GTK_WINDOW (view),
						GTK_DIALOG_DESTROY_WITH_PARENT,
						GTK_MESSAGE_WARNING,
						GTK_BUTTONS_CLOSE,
						"The following variables are absent\n"
						"from the simulation file:\n%s",
						missing);
      gtk_dialog_run (GTK_DIALOG (msg));
      gtk_widget_destroy (msg);
    }
    g_free (missing);

    gfs_domain_cell_traverse (GFS_DOMAIN (sim),
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, sim);
    data[0] = domain;
    data[1] = indices;
    gfs_domain_combine_traverse (domain, GFS_DOMAIN (sim),
				 (FttCellCombineTraverseFunc) substract, data,
				 (FttCellTraverseFunc) outside, data);
    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    i = domain->variables;
    while (i) {
      gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, i->data);
      i = i->next;
    }

    g_slist_free (indices);
    gts_object_destroy (GTS_OBJECT (sim));

    basename = g_path_get_basename (fname);
    title = g_strjoin (" - ", gtk_window_get_title (GTK_WINDOW (view)), basename, NULL);
    if ((s = strstr (title, "GfsView: ")))
      s += strlen ("GfsView: ");
    else
      s = title;
    gfk_gl_view_set_simulation (view, GFS_SIMULATION (domain), s);
    g_free (title);
    g_free (basename);
  }
  return TRUE;
}

G_LOCK_DEFINE (scripting_pending);

gboolean gfk_receive_scripting_message (gpointer data)
{
  GfkScriptingMessage * msg = data;
  switch (msg->event) {
  case GFS_SAVE_EVENT: case GFS_APPEND_EVENT: {
    GfsGl2PSParams * p = msg->data;
    gfs_gl2ps (p, p->fp, "", msg->view);
    if (p->fp == stdout || p->fp == stderr || (msg->event == GFS_APPEND_EVENT))
      fflush (p->fp);
    else
      fclose (p->fp);
    break;
  }
  case GFS_ECHO_EVENT:
    puts (msg->data);
    fflush (stdout);
    break;
  }

  g_free (msg->data);
  g_free (msg);
  G_UNLOCK (scripting_pending);
  return FALSE;
}

static void on_gl_list_destroy (GtkTreeView * tree)
{
  GtkTreeModel * model = gtk_tree_view_get_model (tree);
  GtkTreeIter iter;
  gboolean valid;
  
  valid = gtk_tree_model_get_iter_first (model, &iter);
  while (valid) {
    GfkGl * gl;

    gtk_tree_model_get (model, &iter, GL_COLUMN, &gl, -1);
    gts_object_destroy (GTS_OBJECT (gl));
    valid = gtk_tree_model_iter_next (model, &iter);
  }
}

void gfk_gl_view_set_scripting (GtkWidget * view, gboolean active)
{
  g_return_if_fail (view != NULL);

  GtkWidget * pref = lookup_widget (view, "preferences");
  if (active) {
    gtk_widget_set_sensitive (lookup_widget (pref, "scripting_label"), TRUE);
    gtk_widget_set_sensitive (lookup_widget (pref, "scripting_on"), TRUE);
    gtk_widget_set_sensitive (lookup_widget (pref, "scripting_off"), TRUE);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (lookup_widget (pref, "scripting_on")), TRUE);
  }
  else {
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (lookup_widget (pref, "scripting_off")), TRUE);
    gtk_widget_set_sensitive (lookup_widget (pref, "scripting_label"), FALSE);
    gtk_widget_set_sensitive (lookup_widget (pref, "scripting_on"), FALSE);
    gtk_widget_set_sensitive (lookup_widget (pref, "scripting_off"), FALSE);
  }
  G_LOCK (gfk_gl_scripting);
  gfk_gl_scripting = active;
  G_UNLOCK (gfk_gl_scripting);
}

GtkWidget * gfk_gl_view (GtkWidget * glarea)
{
  GtkWidget * view, * tree, * pref, * colorsel, * colorb, * gl2ps, * about, * filew, * size;
  GtkListStore * store;
  GtkTreeViewColumn * column;
  GtkToolbar * toolbar;
  GtkTreeSelection * select;
  GtkCellRenderer * renderer;
  GfsGlViewParams * p;
  GfsGl2PSParams * q;
  gdouble * ratio;
  GdkColor c;

  GtkMenu * vector, * scalar, * mesh, * special;

  g_return_val_if_fail (glarea != NULL, NULL);

  p = gfs_gl_view_params_new ();
  g_object_set_data_full (G_OBJECT (glarea), "GfsGlViewParams", p, g_free);

  view = create_view ();
  gtk_container_add (GTK_CONTAINER (lookup_widget (view, "frame1")), glarea);
  g_object_set_data (G_OBJECT (view), "glarea", glarea);

  tree = lookup_widget (view, "gl_list");
  g_signal_connect (G_OBJECT (tree), "destroy", (GCallback) on_gl_list_destroy, NULL);
  store = gtk_list_store_new (N_COLUMNS,
			      G_TYPE_BOOLEAN,
			      GDK_TYPE_PIXBUF,
			      G_TYPE_STRING,
			      G_TYPE_POINTER,
			      G_TYPE_BOOLEAN);
  g_object_set_data (G_OBJECT (glarea), "list", store);
  gtk_tree_view_set_model (GTK_TREE_VIEW (tree), GTK_TREE_MODEL (store));
  select = gtk_tree_view_get_selection (GTK_TREE_VIEW (tree));
  gtk_tree_selection_set_mode (select, GTK_SELECTION_SINGLE);
  g_signal_connect (tree, "row-activated", (GCallback) on_properties1_activate, NULL);
  gtk_tree_selection_set_select_function (select, tree_selection_func, NULL, NULL);
  g_signal_connect (G_OBJECT (select), "changed", G_CALLBACK (tree_selection_changed_cb), tree);

  renderer = gtk_cell_renderer_pixbuf_new ();
  column = gtk_tree_view_column_new_with_attributes ("Icon", renderer,
  						     "pixbuf", ICON_COLUMN, 
						     NULL);
  gtk_tree_view_append_column (GTK_TREE_VIEW (tree), column);

  renderer = gtk_cell_renderer_toggle_new ();
  g_signal_connect (G_OBJECT (renderer), "toggled", G_CALLBACK (visible_toggled), store);
  column = gtk_tree_view_column_new_with_attributes ("Visible", renderer,
  						     "active", VISIBLE_COLUMN,
						     NULL);
  gtk_tree_view_append_column (GTK_TREE_VIEW (tree), column);

  renderer = gtk_cell_renderer_text_new ();
  column = gtk_tree_view_column_new_with_attributes ("Properties", renderer,
						     "text", PROPERTIES_COLUMN, 
						     NULL);
  gtk_tree_view_append_column (GTK_TREE_VIEW (tree), column);

  vector = GTK_MENU (gtk_menu_new ());
  gtk_menu_item_set_submenu (GTK_MENU_ITEM (lookup_widget (view, "vector1")), 
			     GTK_WIDGET (vector));
  scalar = GTK_MENU (gtk_menu_new ());
  gtk_menu_item_set_submenu (GTK_MENU_ITEM (lookup_widget (view, "scalar1")), 
			     GTK_WIDGET (scalar));
  mesh = GTK_MENU (gtk_menu_new ());
  gtk_menu_item_set_submenu (GTK_MENU_ITEM (lookup_widget (view, "mesh1")), 
			     GTK_WIDGET (mesh));
  special = GTK_MENU (gtk_menu_new ());
  gtk_menu_item_set_submenu (GTK_MENU_ITEM (lookup_widget (view, "special1")), 
			     GTK_WIDGET (special));

  toolbar = GTK_TOOLBAR (lookup_widget (view, "toolbar"));
  gfk_gl_tool_append (gfk_gl_vectors_class (), toolbar, vector, tree, glarea);
  gfk_gl_tool_append (gfk_gl_streamlines_class (), toolbar, vector, tree, glarea);
  gfk_gl_tool_append (gfk_gl_squares_class (), toolbar, scalar, tree, glarea);
  gfk_gl_tool_append (gfk_gl_linear_class (), toolbar, scalar, tree, glarea);
  gfk_gl_tool_append (gfk_gl_isoline_class (), toolbar, scalar, tree, glarea);
  gfk_gl_tool_append (gfk_gl_cells_class (), toolbar, mesh, tree, glarea);
  gfk_gl_tool_append (gfk_gl_solid_class (), toolbar, mesh, tree, glarea);
  gfk_gl_tool_append (gfk_gl_vof_class (), toolbar, scalar, tree, glarea);

  gfk_gl_menu_append (gfk_gl_boundaries_class (), mesh, tree, glarea);
  gfk_gl_menu_append (gfk_gl_levels_class (), mesh, tree, glarea);
  gfk_gl_menu_append (gfk_gl_fractions_class (), mesh, tree, glarea);
  gfk_gl_menu_append (gfk_gl_symmetry_class (), mesh, tree, glarea);
  gfk_gl_menu_append (gfk_gl_periodic_class (), mesh, tree, glarea);
  gfk_gl_menu_append (gfk_gl_height_class (), scalar, tree, glarea);
  gfk_gl_menu_append (gfk_gl_location_class (), special, tree, glarea);
  gfk_gl_menu_append (gfk_gl_ellipses_class (), special, tree, glarea);
  gfk_gl_menu_append (gfk_gl_info_class (), special, tree, glarea);
  gfk_gl_menu_append (gfk_gl_label_class (), special, tree, glarea);
  gfk_gl_menu_append (gfk_gl_locate_class (), special, tree, glarea);
  gfk_gl_menu_append (gfk_gl_clip_plane_class (), special, tree, glarea);
  gfk_gl_menu_append (gfk_gl_pipes_class (), special, tree, glarea);

#if (!FTT_2D)
  gfk_gl_tool_append (gfk_gl_isosurface_class (), toolbar, scalar, tree, glarea);
  gfk_gl_menu_append (gfk_gl_cut_plane_class (), special, tree, glarea);
#endif /* 3D */

  about = create_about ();
  gtk_window_set_transient_for (GTK_WINDOW (about), GTK_WINDOW (view));
  gtk_window_set_destroy_with_parent (GTK_WINDOW (about), TRUE);
  gtk_window_set_position (GTK_WINDOW (about), GTK_WIN_POS_CENTER_ON_PARENT);
  g_signal_connect (G_OBJECT (about), "delete_event", 
		    G_CALLBACK (gtk_widget_hide_on_delete), NULL);
  g_object_set_data (G_OBJECT (view), "about", about);

  filew = gtk_file_selection_new ("Select a Gerris simulation");
  gtk_window_set_transient_for (GTK_WINDOW (filew), GTK_WINDOW (view));
  gtk_window_set_destroy_with_parent (GTK_WINDOW (filew), TRUE);
  gtk_window_set_position (GTK_WINDOW (filew), GTK_WIN_POS_CENTER_ON_PARENT);
  g_signal_connect (G_OBJECT (filew), "delete_event", 
		    G_CALLBACK (gtk_widget_hide_on_delete), NULL);
  g_object_set_data (G_OBJECT (filew), "view", view);
  g_signal_connect (G_OBJECT (GTK_FILE_SELECTION (filew)->ok_button), "clicked",
		    G_CALLBACK (filew_ok), filew);
  g_signal_connect_swapped (G_OBJECT (GTK_FILE_SELECTION (filew)->cancel_button), "clicked", 
			    G_CALLBACK (gtk_widget_hide), filew);
  g_object_set_data (G_OBJECT (view), "filew", filew);

  filew = gtk_file_selection_new ("Select a Gerris simulation");
  gtk_window_set_transient_for (GTK_WINDOW (filew), GTK_WINDOW (view));
  gtk_window_set_destroy_with_parent (GTK_WINDOW (filew), TRUE);
  gtk_window_set_position (GTK_WINDOW (filew), GTK_WIN_POS_CENTER_ON_PARENT);
  g_signal_connect (G_OBJECT (filew), "delete_event", 
		    G_CALLBACK (gtk_widget_hide_on_delete), NULL);
  g_object_set_data (G_OBJECT (filew), "view", view);
  g_signal_connect (G_OBJECT (GTK_FILE_SELECTION (filew)->ok_button), "clicked",
		    G_CALLBACK (filews_ok), filew);
  g_signal_connect_swapped (G_OBJECT (GTK_FILE_SELECTION (filew)->cancel_button), "clicked", 
			    G_CALLBACK (gtk_widget_hide), filew);
  g_object_set_data (G_OBJECT (view), "filews", filew);

  pref = create_preferences ();
  gtk_window_set_transient_for (GTK_WINDOW (pref), GTK_WINDOW (view));
  gtk_window_set_destroy_with_parent (GTK_WINDOW (pref), TRUE);
  gtk_window_set_position (GTK_WINDOW (pref), GTK_WIN_POS_CENTER_ON_PARENT);
  g_signal_connect (G_OBJECT (pref), "delete_event", G_CALLBACK (gtk_widget_hide_on_delete), NULL);
  g_object_set_data (G_OBJECT (view), "preferences", pref);
  g_object_set_data (G_OBJECT (pref), "view", view);
  gfk_gl_view_set_scripting (view, FALSE);
#if FTT_2D
  gtk_widget_hide (lookup_widget (pref, "sz_label"));
  gtk_widget_hide (lookup_widget (pref, "sz"));
#endif /* 2D */

  colorsel = create_bg_color_selector ();
  gtk_window_set_transient_for (GTK_WINDOW (colorsel), GTK_WINDOW (pref));
  gtk_window_set_destroy_with_parent (GTK_WINDOW (colorsel), TRUE);
  gtk_window_set_position (GTK_WINDOW (colorsel), GTK_WIN_POS_CENTER_ON_PARENT);
  g_signal_connect (G_OBJECT (colorsel), "delete_event", 
		    G_CALLBACK (gtk_widget_hide_on_delete), NULL);
  c.red = p->bg.r*65535; c.green = p->bg.g*65535; c.blue = p->bg.b*65535;
  gtk_color_selection_set_current_color (GTK_COLOR_SELECTION (lookup_widget (colorsel,
									     "colorselection2")), 
					 &c);

  colorb = lookup_widget_params (pref, "bgcolor");
  g_object_set_data (G_OBJECT (colorb), "colorsel", colorsel);
  gtk_widget_modify_bg (colorb, GTK_STATE_NORMAL, &c);
  gtk_widget_modify_bg (colorb, GTK_STATE_PRELIGHT, &c);
  gtk_widget_modify_bg (colorb, GTK_STATE_ACTIVE, &c);
  g_object_set_data (G_OBJECT (colorsel), "bgcolor", colorb);

  gl2ps = create_gl2ps ();
  gtk_window_set_transient_for (GTK_WINDOW (gl2ps), GTK_WINDOW (view));
  gtk_window_set_destroy_with_parent (GTK_WINDOW (gl2ps), TRUE);
  gtk_window_set_position (GTK_WINDOW (gl2ps), GTK_WIN_POS_CENTER_ON_PARENT);
  g_signal_connect (G_OBJECT (gl2ps), "delete_event",
		    G_CALLBACK (gtk_widget_hide_on_delete), NULL);
  g_object_set_data (G_OBJECT (view), "gl2ps", gl2ps);
  g_object_set_data (G_OBJECT (gl2ps), "view", view);

  size = lookup_widget (gl2ps, "size");
  gtk_widget_set_sensitive (size, FALSE);
  g_signal_connect (G_OBJECT (glarea), "expose_event", G_CALLBACK (gl2ps_update_ppm_size), size);
  ratio = g_malloc (sizeof (gdouble));
  g_object_set_data_full (G_OBJECT (size), "ratio", ratio, g_free);
  
  q = g_malloc (sizeof (GfsGl2PSParams));
  q->width = q->height = 0;
  q->format = GFSGL_PPM;
  q->lw = 1.;
  gl2ps_ppm_set_sensitive (gl2ps, FALSE, TRUE);
  q->sort = GL2PS_SIMPLE_SORT;
  gtk_option_menu_set_history (GTK_OPTION_MENU (lookup_widget (gl2ps, "sort")), 1);
  q->options = (GL2PS_SIMPLE_LINE_OFFSET |
		GL2PS_SILENT |
		GL2PS_BEST_ROOT |
		GL2PS_OCCLUSION_CULL |
		GL2PS_USE_CURRENT_VIEWPORT |
		GL2PS_TIGHT_BOUNDING_BOX);
  g_object_set_data_full (G_OBJECT (gl2ps), "GfsGl2PSParams", q, g_free);

  gtk_widget_set_sensitive (lookup_widget (view, "toolbar"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (view, "objects1"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (view, "save1"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (view, "edit1"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (view, "properties1"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (view, "delete3"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (view, "view1"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (view, "tools1"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (view, "gl_list"), FALSE);

  GtkListStore * completion = gtk_list_store_new (3, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_POINTER);
  g_object_set_data_full (G_OBJECT (view), "completion",
			  g_object_ref (G_OBJECT (completion)), (GDestroyNotify) g_object_unref);

  return view;
}

void gfk_gl_view_set_simulation (GtkWidget * view, GfsSimulation * sim, const gchar * fname)
{
  GtkWidget * glarea;
  GfsSimulation * prev;
  gchar * basename, * s;
  static GStaticMutex mutex = G_STATIC_MUTEX_INIT;
  GtkTreeModel * list;
  GtkTreeIter iter;
  gboolean valid;

  g_return_if_fail (view != NULL);
  g_return_if_fail (sim != NULL);
  g_return_if_fail (fname != NULL);

  g_static_mutex_lock (&mutex);

  glarea = lookup_widget (view, "glarea");
  prev = g_object_get_data (G_OBJECT (glarea), "sim");
  g_object_set_data (G_OBJECT (glarea), "sim", sim);
  list = gtk_tree_view_get_model (GTK_TREE_VIEW (lookup_widget (view, "gl_list")));
  valid = gtk_tree_model_get_iter_first (list, &iter);
  while (valid) {
    GfkGl * gl;
    
    gtk_tree_model_get (list, &iter, GL_COLUMN, &gl, -1);
    gfk_gl_set_simulation (gl, sim);
    valid = gtk_tree_model_iter_next (list, &iter);
  }

  /* Remove obsolete variables from completion */
  list = GTK_TREE_MODEL (lookup_widget (view, "completion"));
  valid = gtk_tree_model_get_iter_first (list, &iter);
  while (valid) {
    gchar * name;
    gpointer p;
    gtk_tree_model_get (list, &iter, 
			COMPLETION_NAME_COLUMN, &name,
			COMPLETION_POINTER_COLUMN, &p, 
			-1);
    if (p) {
      p = gfs_variable_from_name (GFS_DOMAIN (sim)->variables, name);
      if (p) {
	gtk_list_store_set (GTK_LIST_STORE (list), &iter,
			    COMPLETION_NAME_COLUMN, GFS_VARIABLE (p)->name,
			    COMPLETION_DESCRIPTION_COLUMN, GFS_VARIABLE (p)->description,
			    COMPLETION_POINTER_COLUMN, p,
			    -1);
	gfs_object_simulation_set (p, NULL);
	valid = gtk_tree_model_iter_next (list, &iter);
      }
      else if ((p = gfs_derived_variable_from_name (GFS_DOMAIN (sim)->derived_variables, name))) {
	gtk_list_store_set (GTK_LIST_STORE (list), &iter,
			    COMPLETION_NAME_COLUMN, GFS_DERIVED_VARIABLE (p)->name,
			    COMPLETION_DESCRIPTION_COLUMN, GFS_DERIVED_VARIABLE (p)->description,
			    COMPLETION_POINTER_COLUMN, p,
			    -1);
       	gfs_object_simulation_set (p, NULL);
	valid = gtk_tree_model_iter_next (list, &iter);
      }
      else
	valid = gtk_list_store_remove (GTK_LIST_STORE (list), &iter);
    }
    else
      valid = gtk_tree_model_iter_next (list, &iter);
  }

  /* Add new variables to completion*/
  GSList * i = GFS_DOMAIN (sim)->variables;
  while (i) {
    GfsVariable * v = i->data;
    if (v->name) {
      if (gfs_object_simulation (v)) {
	gtk_list_store_append (GTK_LIST_STORE (list), &iter);
	gtk_list_store_set (GTK_LIST_STORE (list), &iter, 
			    COMPLETION_NAME_COLUMN, v->name,
			    COMPLETION_DESCRIPTION_COLUMN, v->description,
			    COMPLETION_POINTER_COLUMN, v,
			    -1);
      }
      else
	gfs_object_simulation_set (v, sim);
    }
    i = i->next;
  }

  /* Add new derived variables to completion */
  i = GFS_DOMAIN (sim)->derived_variables;
  while (i) {
    GfsDerivedVariable * v = i->data;
    if (v->name) {
      if (gfs_object_simulation (v)) {
	gtk_list_store_append (GTK_LIST_STORE (list), &iter);
	gtk_list_store_set (GTK_LIST_STORE (list), &iter, 
			    COMPLETION_NAME_COLUMN, v->name,
			    COMPLETION_DESCRIPTION_COLUMN, v->description,
			    COMPLETION_POINTER_COLUMN, v,
			    -1);
      }
      else
	gfs_object_simulation_set (v, sim);
    }
    i = i->next;
  }

  if (prev == NULL) {
    gtk_widget_set_sensitive (lookup_widget (view, "toolbar"), TRUE);
    gtk_widget_set_sensitive (lookup_widget (view, "objects1"), TRUE);
    gtk_widget_set_sensitive (lookup_widget (view, "tools1"), TRUE);
    gtk_widget_set_sensitive (lookup_widget (view, "gl_list"), TRUE);
  }
  else if (prev != sim)
    gts_object_destroy (GTS_OBJECT (prev));
  
  basename = g_path_get_basename (fname);
  s = g_strdup_printf ("GfsView: %s: t = %g", basename, sim->time.t);
  gtk_window_set_title (GTK_WINDOW (view), s);
  g_free (s);
  g_free (basename);

  g_static_mutex_unlock (&mutex);
}

void gfk_gl_view_draw (GtkWidget * view, guint format)
{
  GtkWidget * glarea;
  GfsGlViewParams * p;
  GtkTreeModel * list;
  GtkTreeIter iter;
  gboolean valid;
  guint size = 0;
  GTimer * timer;

  g_return_if_fail (view != NULL);

  timer = g_timer_new ();
  g_timer_start (timer);

  glarea = g_object_get_data (G_OBJECT (view), "glarea");
  p = g_object_get_data (G_OBJECT (glarea), "GfsGlViewParams");

  if (p->motion) {
    if (p->timing > p->reactivity) {
      if (p->res == 0.) p->res = 1.;
      p->res *= exp (ceil (log (p->timing/p->reactivity)/log(4.))*log (2.));
      if (p->res > 50.)
	p->res = 50.;
    }
  }
  else
    p->res = p->base_res;

  list = gtk_tree_view_get_model (GTK_TREE_VIEW (lookup_widget (view, "gl_list")));
  GList * symmetries = get_symmetries (list);
  GfsFrustum frustum;
  gfs_gl_get_frustum (p, symmetries, &frustum);
  GLuint display_list = glGenLists (1);
  GList * gl_list = NULL;
  valid = gtk_tree_model_get_iter_first (list, &iter);
  while (valid) {
    gboolean visible;
    GfkGl * gl;

    gtk_tree_model_get (list, &iter, GL_COLUMN, &gl, VISIBLE_COLUMN, &visible, -1);
    if (GFK_IS_GL_CLIP_PLANE (gl)) {
      gfs_gl_clip_plane_disable (GFS_GL_CLIP_PLANE (gl->gl));
      GFS_GL_CLIP_PLANE (gl->gl)->disabled = !visible;
    }
    else if (visible || !GFS_IS_GL_CUT_PLANE (gl->gl))
      gl_list = g_list_append (gl_list, gl->gl);
    valid = gtk_tree_model_iter_next (list, &iter);
  }

  glNewList (display_list, GL_COMPILE);
  GList * i = gl_list;
  while (i) {
    GfsGl * gl = i->data;
    if (GFS_IS_GL_CUT_PLANE (gl)) {
      GFS_GL_CUT_PLANE (gl)->list = gl_list;
      gl->format = format;
      gfs_gl_draw (gl, &frustum);
      if (gl->size > size)
	size = gl->size;
      GFS_GL_CUT_PLANE (gl)->list = NULL;
    }
    i = i->next;
  }
  g_list_free (gl_list);

  GSList * clip = NULL;
  gboolean firstclip = TRUE;
  valid = gtk_tree_model_get_iter_first (list, &iter);
  while (valid) {
    gboolean visible;
    GfkGl * gl;

    gtk_tree_model_get (list, &iter, GL_COLUMN, &gl, VISIBLE_COLUMN, &visible, -1);
    if (visible) {
      gl->gl->format = format;
      if (GFK_IS_GL_CLIP_PLANE (gl)) {
	if (firstclip) {
	  g_slist_foreach (clip, (GFunc) gfs_gl_clip_plane_disable, NULL);
	  g_slist_free (clip); clip = NULL;
	  firstclip = FALSE;	  
	}
	gfs_gl_draw (gl->gl, &frustum);
	clip = g_slist_prepend (clip, gl->gl);
      }
      else {
	gfs_gl_draw (gl->gl, &frustum);
	if (gl->gl->size > size)
	  size = gl->gl->size;
	firstclip = TRUE;
      }
    }
    else if (!GFK_IS_GL_CLIP_PLANE (gl))
      firstclip = TRUE;
    valid = gtk_tree_model_iter_next (list, &iter);
  }
  g_slist_free (clip);
  glEndList();

  gfs_gl_symmetry_apply (symmetries, display_list);
  gfs_gl_frustum_free (&frustum);
  g_list_free (symmetries);
  glDeleteLists (display_list, 1);  

  g_timer_stop (timer);
  p->timing = g_timer_elapsed (timer, NULL);
  g_timer_destroy (timer);

  if (size > 0) {
    GtkStatusbar * status = GTK_STATUSBAR (lookup_widget (view, "statusbar1"));
    guint id = gtk_statusbar_get_context_id (status, "glarea");
    gchar * s = g_strdup_printf (" %d items, %g s", size, p->timing);

    gtk_statusbar_pop (status, id);
    gtk_statusbar_push (status, id, s);
    g_free (s);
  }
}

void gfk_gl_view_pick (GtkWidget * view, GfsGlRay * ray, gboolean motion)
{
  GtkWidget * tree;
  GtkTreeModel * list;
  GtkTreeSelection * select;
  GtkTreeIter iter, picked;
  gboolean valid;
  GfkGl * glmin = NULL;
  gdouble zmin = G_MAXDOUBLE/2.;

  GtkStatusbar * status;
  guint id;

  g_return_if_fail (view != NULL);
  g_return_if_fail (ray != NULL);

  tree = lookup_widget (view, "gl_list");
  list = gtk_tree_view_get_model (GTK_TREE_VIEW (tree));
  select = gtk_tree_view_get_selection (GTK_TREE_VIEW (tree));
  valid = gtk_tree_model_get_iter_first (list, &iter);
  while (valid) {
    gboolean visible;
    GfkGl * gl;

    gtk_tree_model_get (list, &iter, GL_COLUMN, &gl, VISIBLE_COLUMN, &visible, -1);
    if (visible && GFS_GL_CLASS (GTS_OBJECT (gl->gl)->klass)->pick) {
      gdouble z = (* GFS_GL_CLASS (GTS_OBJECT (gl->gl)->klass)->pick) (gl->gl, ray);

      if (z < zmin || (z == zmin && gtk_tree_selection_iter_is_selected (select, &iter))) {
	glmin = gl;
	zmin = z;
	picked = iter;
      }
    }
    valid = gtk_tree_model_iter_next (list, &iter);
  }

  status = GTK_STATUSBAR (lookup_widget (view, "statusbar1"));
  id = gtk_statusbar_get_context_id (status, "pick");
  gtk_statusbar_pop (status, id);
  if (glmin) {
    gtk_tree_selection_select_iter (select, &picked);
    g_assert (GFK_GL_CLASS (GTS_OBJECT (glmin)->klass)->pickinfo);
    gtk_statusbar_push (status, id, 
			(* GFK_GL_CLASS (GTS_OBJECT (glmin)->klass)->pickinfo) (glmin, motion));
  }
}

void gfk_gl_view_clear (GtkWidget * view)
{
  GtkWidget * list, * glarea;
  GtkTreeModel * model;

  g_return_if_fail (view != NULL);

  list = lookup_widget (view, "gl_list");
  glarea = lookup_widget (view, "glarea");
  model = gtk_tree_view_get_model (GTK_TREE_VIEW (list));

  on_gl_list_destroy (GTK_TREE_VIEW (list));
  gtk_list_store_clear (GTK_LIST_STORE (model));
  g_object_set_data (G_OBJECT (list), "former", NULL);

  gtk_widget_set_sensitive (lookup_widget (glarea, "save1"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (glarea, "edit1"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (glarea, "properties1"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (glarea, "delete3"), FALSE);
  gtk_widget_set_sensitive (lookup_widget (glarea, "view1"), FALSE);

  gdk_window_invalidate_rect (glarea->window, &glarea->allocation, FALSE);
}
