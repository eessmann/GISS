/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2004 National Institute of Water and Atmospheric
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

/* GfkGl2D: Object */

static gchar * gl2D_pickinfo (GfkGl * gl, gboolean motion)
{
  FttVector p = GFS_GL2D (gl->gl)->pickedpos;
  GfkGl2D * gl2 = GFK_GL2D (gl);

  g_free (gl2->pickinfo);
  if (fabs (p.x) < 1e-10) p.x = 0.;
  if (fabs (p.y) < 1e-10) p.y = 0.;
  if (fabs (p.z) < 1e-10) p.z = 0.;
  gfs_simulation_map_inverse (gl->gl->sim, &p);
  gl2->pickinfo = g_strdup_printf ("(%.8f,%.8f,%.8f)", p.x, p.y, p.z);
  return gl2->pickinfo;
}

static void gl2D_init (GfkGl2D * object)
{
  object->params = create_gl2D_params ();
  object->n.x = 0.; object->n.y = 0.; object->n.z = 1.;
  gfk_gl_prepend_params (GFK_GL (object), object->params, gtk_label_new ("2D Plane"));
  gtk_widget_show (object->params);
}

/* GfkGlSolid: Object */

static void update_solid_scalar (GfkGl * gl)
{
  gfk_gl_set_sensitive (gl, GFK_GL_SCALAR (gl)->scalar, TRUE);
  gfk_gl_update_properties (gl);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color_label"), 
			    FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color"), 
			    FALSE);
}

static void gl_solid_read (GtsObject ** o, GtsFile * fp)
{
  GtkWidget * solid = GFK_GL_SOLID (*o)->solid;
  GfsGlSolid * gl = GFS_GL_SOLID (GFK_GL (*o)->gl);

  (* GTS_OBJECT_CLASS (gfk_gl_solid_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (gl->reversed) {
    gl->reversed = FALSE;
    gtk_toggle_button_set_active (lookup_widget_params (solid, "reversed"), TRUE);
  }

  if (gl->use_scalar) {
    gtk_option_menu_set_history (lookup_widget_params (solid, "color"), 1);
    update_solid_scalar (GFK_GL (*o));
  }
}

static gchar * gl_solid_name (GfkGlClass * klass)
{
  static gchar name[] = "Solid";
  return name;
}

static GtkWidget * gl_solid_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "solid-16x16.png");
}

static gchar * gl_solid_properties (GfkGl * gl)
{
  gchar * s = GFS_GL_SOLID (gl->gl)->use_scalar ? 
    g_strdup ((* gfk_gl_scalar_class ()->properties) (gl)) :
    g_strdup ("");

  g_free (gl->props);
  gl->props = s;
  return gl->props;
}

static void set_solid_color (GtkWidget * color, GfkGl * gl)
{
  gfk_gl_set_sensitive (gl, GFK_GL_SCALAR (gl)->scalar, FALSE);
  gfk_gl_update_properties (gl);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color_label"), 
			    TRUE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color"), TRUE);
  GFS_GL_SOLID (gl->gl)->use_scalar = NULL;
  gfk_gl_expose (gl);
}

static void set_solid_scalar (GtkWidget * scalar, GfkGl * gl)
{
  GFS_GL_SOLID (gl->gl)->use_scalar = GFS_GL_SCALAR (gl->gl)->v;
  update_solid_scalar (gl);
  gfs_gl_solid_reset (GFS_GL_SOLID (gl->gl));
  gfk_gl_expose (gl);
}

static void gl_solid_post_init (GfkGl * object)
{
  GtkWidget * m = gtk_menu_new (), * i;
  GfkGlSolid * gls = GFK_GL_SOLID (object);
  
  if (GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_solid_class ())->parent_class)->post_init)
    (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_solid_class ())->parent_class)->post_init)
      (object);
  
  gfk_gl_set_sensitive (GFK_GL (gls), GFK_GL_SCALAR (gls)->scalar, FALSE);
  m = gtk_menu_new ();
  i = gtk_menu_item_new_with_label ("Default");
  g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_solid_color), object);
  gtk_menu_append (m, i);
  gtk_widget_show (i);

  i = gtk_menu_item_new_with_label ("Scalar");
  g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_solid_scalar), object);
  gtk_menu_append (m, i);
  gtk_widget_show (i);
  
  gtk_option_menu_set_menu (lookup_widget_params (gls->solid, "color"), m);
  gtk_widget_show (m);

  gtk_widget_hide (GFK_GL2D (object)->params);
}

static void gl_solid_class_init (GfkGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_solid_read;
  klass->gl_class = gfs_gl_solid_class ();
  klass->post_init = gl_solid_post_init;
  klass->name = gl_solid_name;
  klass->icon = gl_solid_icon;
  klass->properties = gl_solid_properties;
}

static void gl_solid_init (GfkGl * gl)
{
  GfkGlSolid * gli = GFK_GL_SOLID (gl);

  gli->solid = create_solid_params ();
  gfk_gl_prepend_params (gl, gli->solid, gtk_label_new ("Solid"));
  gtk_widget_show (gli->solid);
}

GfkGlClass * gfk_gl_solid_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_solid_info = {
      "GfkGlSolid",
      sizeof (GfkGlSolid),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_solid_class_init,
      (GtsObjectInitFunc) gl_solid_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_scalar_class ()),
				  &gfk_gl_solid_info);
  }

  return klass;
}

/* GfkGlFractions: Object */

static void gl_fractions_init (GfkGl * gl)
{
}

/* GfkGlIsosurface: Object */

static void update_isosurface_scalar (GfkGl * gl)
{
  gfk_gl_set_sensitive (gl, GFK_GL_SCALAR (gl)->scalar, TRUE);
  gfk_gl_update_properties (gl);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "default_color"), FALSE);
}

static void gl_isosurface_read (GtsObject ** o, GtsFile * fp)
{
  GtkWidget * scalar = GFK_GL_ISOSURFACE (*o)->scalar;
  GfsGlIsosurface * gl = GFS_GL_ISOSURFACE (GFK_GL (*o)->gl);
  GtkSpinButton * level = lookup_widget_params (scalar, "spinbuttonlevel");
  GtkAdjustment * alevel = gtk_spin_button_get_adjustment (level);

  (* GTS_OBJECT_CLASS (gfk_gl_isosurface_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  alevel->value = gl->level;
  if (alevel->value > alevel->upper)
    alevel->upper = alevel->value;
  if (alevel->value < alevel->lower)
    alevel->lower = alevel->value;
  gtk_spin_button_set_adjustment (level, alevel);

  if (gl->reversed) {
    gl->reversed = FALSE;
    gtk_toggle_button_set_active (lookup_widget_params (scalar, "reversed"), TRUE);
  }

  if (gl->use_scalar) {
    gtk_option_menu_set_history (lookup_widget_params (scalar, "color"), 1);
    update_isosurface_scalar (GFK_GL (*o));
  }
}

static gchar * gl_isosurface_name (GfkGlClass * klass)
{
  static gchar name[] = "Isosurface";
  return name;
}

static GtkWidget * gl_isosurface_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "isosurface-16x16.png");
}

static void gl_isosurface_update_interface (GfkGl * g)
{
  GfkGlIsosurface * gl = GFK_GL_ISOSURFACE (g);
  GtkSpinButton * level = lookup_widget_params (gl->scalar, "spinbuttonlevel");
  GfsGlIsosurface * gli = GFS_GL_ISOSURFACE (GFK_GL (gl)->gl);
  GtkEntry * scalar = lookup_widget_params (gl->scalar, "scalar");

  if (!GFK_IS_EDITED (scalar))
    gtk_entry_set_text (scalar, gli->expr->str);

  if (gli->maxv > gli->minv) {
    gdouble step = (gli->maxv - gli->minv)/500.;
    guint digits = gfk_decimal_digits (step, 2);
    GtkAdjustment * alevel = gtk_spin_button_get_adjustment (level);

    alevel->upper = gli->maxv;
    alevel->lower = gli->minv;
    alevel->step_increment = step;
    alevel->page_increment = step*10.;
    if (alevel->value > alevel->upper || alevel->value < alevel->lower)
      alevel->value = (alevel->upper + alevel->lower)/2.;
    if (!GFK_IS_EDITED (level))
      gtk_spin_button_configure (level, alevel, 2.*step, digits);
  }
  gfk_gl_update_properties (GFK_GL (gl));
  gfk_gl_expose (GFK_GL (gl));
}

static void gl_isosurface_set_simulation (GfkGl * object, GfsSimulation * sim)
{
  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_isosurface_class ())->parent_class)->set_simulation)
    (object, sim);
  
  gl_isosurface_update_interface (object);
}

static void set_isosurface_scalar (GtkWidget * scalar, GfkGl * gl)
{
  GFS_GL_ISOSURFACE (gl->gl)->use_scalar = GFS_GL_SCALAR (gl->gl)->v;
  update_isosurface_scalar (gl);
  gfs_gl_isosurface_reset (GFS_GL_ISOSURFACE (gl->gl));
  gfk_gl_expose (gl);
}

static void set_isosurface_color (GtkWidget * color, GfkGlIsosurface * gl)
{
  GFS_GL_ISOSURFACE (GFK_GL (gl)->gl)->use_scalar = NULL;
  gfk_gl_set_sensitive (GFK_GL (gl), GFK_GL_SCALAR (gl)->scalar, FALSE);
  gfk_gl_update_properties (GFK_GL (gl));
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color_label"), 
			    TRUE);
  gtk_widget_set_sensitive (lookup_widget_params (GFK_GL (gl)->properties, "default_color"), TRUE);
  gfk_gl_expose (GFK_GL (gl));
}

static GtsFile * gl_isosurface_set (GfkGl * gl, const char * s)
{
  return gfs_gl_isosurface_set (GFS_GL_ISOSURFACE (gl->gl), s);
}

static void gl_isosurface_post_init (GfkGl * object)
{
  GfkGlIsosurface * gli = GFK_GL_ISOSURFACE (object);

  (*GFK_GL_CLASS (GTS_OBJECT_CLASS (gfk_gl_isosurface_class ())->parent_class)->post_init)
    (object);

  GtkWidget * scalar = gfk_function (gli->scalar, object, "scalar", gl_isosurface_set);
  gtk_widget_show (scalar);
  gtk_table_attach (GTK_TABLE (lookup_widget_params (gli->scalar, "table1")),
		    scalar, 1, 3, 0, 1,
                    (GtkAttachOptions) (GTK_FILL),
                    (GtkAttachOptions) (0), 0, 0);
  
  gl_isosurface_update_interface (object);
}

static gchar * gl_isosurface_properties (GfkGl * gl)
{
  gchar * s = g_strjoin (" ",
			 GFS_GL_ISOSURFACE (gl->gl)->expr->str,
			 GFS_GL_ISOSURFACE (gl->gl)->use_scalar ?
			 ((* gfk_gl_scalar_class ()->properties) (gl)) : NULL, 
			 NULL);
  g_free (gl->props);
  gl->props = s;
  return gl->props;
}

static void gl_isosurface_class_init (GfkGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_isosurface_read;
  klass->gl_class = gfs_gl_isosurface_class ();
  klass->set_simulation = gl_isosurface_set_simulation;
  klass->update_interface = gl_isosurface_update_interface;
  klass->post_init = gl_isosurface_post_init;
  klass->properties = gl_isosurface_properties;
  klass->name = gl_isosurface_name;
  klass->icon = gl_isosurface_icon;
}

static void gl_isosurface_init (GfkGl * gl)
{
  GfkGlIsosurface * gli = GFK_GL_ISOSURFACE (gl);

  gli->scalar = create_isosurface_params ();

  gfk_gl_set_sensitive (gl, GFK_GL_SCALAR (gl)->scalar, FALSE);

  GtkWidget * m = gtk_menu_new ();
  GtkWidget * i = gtk_menu_item_new_with_label ("Default");
  g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_isosurface_color), gl);
  gtk_menu_append (m, i);
  gtk_widget_show (i);

  i = gtk_menu_item_new_with_label ("Scalar");
  g_signal_connect (G_OBJECT (i), "activate", GTK_SIGNAL_FUNC (set_isosurface_scalar), gl);
  gtk_menu_append (m, i);
  gtk_widget_show (i);
  
  gtk_option_menu_set_menu (lookup_widget_params (gli->scalar, "color"), m);
  gtk_widget_show (m);

  gtk_widget_hide (GFK_GL2D (gl)->params);

  gfk_gl_prepend_params (gl, gli->scalar, gtk_label_new ("Isosurface"));
  gtk_widget_show (gli->scalar);
}

GfkGlClass * gfk_gl_isosurface_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_isosurface_info = {
      "GfkGlIsosurface",
      sizeof (GfkGlIsosurface),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_isosurface_class_init,
      (GtsObjectInitFunc) gl_isosurface_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_scalar_class ()),
				  &gfk_gl_isosurface_info);
  }

  return klass;
}

/* GfkGlCutPlane: Object */

static gchar * gl_cut_plane_name (GfkGlClass * klass)
{
  static gchar name[] = "Cut";
  return name;
}

static GtkWidget * gl_cut_plane_icon (GfkGlClass * klass)
{
  return create_pixmap (NULL, "cut-16x16.png");
}

static void gl_cut_plane_class_init (GfkGlClass * klass)
{
  klass->gl_class = gfs_gl_cut_plane_class ();
  klass->name = gl_cut_plane_name;
  klass->icon = gl_cut_plane_icon;
  klass->pickinfo = NULL;
}

static void gl_cut_plane_init (GfkGl * gl)
{
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading"), FALSE);
}

GfkGlClass * gfk_gl_cut_plane_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_cut_plane_info = {
      "GfkGlCutPlane",
      sizeof (GfkGl2D),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_cut_plane_class_init,
      (GtsObjectInitFunc) gl_cut_plane_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl2D_class ()), &gfk_gl_cut_plane_info);
  }

  return klass;
}
