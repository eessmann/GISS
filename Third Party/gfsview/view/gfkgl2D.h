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
  gfs_simulation_map_inverse (gl->gl->sim, &p);
  gl2->pickinfo = g_strdup_printf ("(%.8f,%.8f)", p.x, p.y);
  return gl2->pickinfo;
}

static void gl2D_init (GfkGl2D * object)
{
  object->params = create_gl2D_params ();
  object->n.x = 0.; object->n.y = 0.; object->n.z = 1.;
  gfk_gl_prepend_params (GFK_GL (object), object->params, gtk_label_new ("2D Plane"));
}

/* GfkGlSolid: Object */

static void gl_solid_init (GfkGl * gl)
{
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading"), FALSE);
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

static void gl_solid_class_init (GfkGlClass * klass)
{
  klass->gl_class = gfs_gl_solid_class ();
  klass->name = gl_solid_name;
  klass->icon = gl_solid_icon;
}

GfkGlClass * gfk_gl_solid_class (void)
{
  static GfkGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfk_gl_solid_info = {
      "GfkGlSolid",
      sizeof (GfkGl),
      sizeof (GfkGlClass),
      (GtsObjectClassInitFunc) gl_solid_class_init,
      (GtsObjectInitFunc) gl_solid_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfk_gl_class ()), &gfk_gl_solid_info);
  }

  return klass;
}

/* GfkGlFractions: Object */

static void gl_fractions_init (GfkGl * gl)
{
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading_label"), FALSE);
  gtk_widget_set_sensitive (lookup_widget_params (gl->properties, "shading"), FALSE);
}

