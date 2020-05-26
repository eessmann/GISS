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

#include "gfsgl.h"
#include <gerris/isocube.h>
#include <gerris/river.h>
#include "trackball.h"
#include "config.h"
#if HAVE_FTGL
# include <FTGL/ftgl.h>
#endif

#define DEFAULT_FONT (PACKAGE_DATA_DIR "/fonts/Garuda.ttf")
#define DEFAULT_FONT_SIZE 72.

static gpointer default_font = NULL, default_raster_font = NULL;

#ifndef GL_CONSTANT
# define GL_CONSTANT GL_CONSTANT_EXT
#endif

static void gl_isoline_update_levels (GfsGl * gl);

typedef struct {
  GfsVariable * min;
  GfsVariable * max;
  gdouble level;
  GfsGl * gl;
} IsoParams;

#if FTT_2D
#  include "gfsgl2D.h"
#else  /* 3D */
#  include "gfsgl3D.h"
#endif /* 3D */
#if defined(__APPLE__)
#  include <OpenGL/glu.h>
#else
#  include <GL/glu.h>
#endif

static void draw_face (GtsTriangle * t)
{
  GtsVertex * v1, * v2, * v3;
  gts_triangle_vertices (t, &v1, &v2, &v3);
  glNormal3d (GTS_POINT (v1)->x, GTS_POINT (v1)->y, GTS_POINT (v1)->z);
  glVertex3d (GTS_POINT (v1)->x, GTS_POINT (v1)->y, GTS_POINT (v1)->z);
  glNormal3d (GTS_POINT (v2)->x, GTS_POINT (v2)->y, GTS_POINT (v2)->z);
  glVertex3d (GTS_POINT (v2)->x, GTS_POINT (v2)->y, GTS_POINT (v2)->z);
  glNormal3d (GTS_POINT (v3)->x, GTS_POINT (v3)->y, GTS_POINT (v3)->z);
  glVertex3d (GTS_POINT (v3)->x, GTS_POINT (v3)->y, GTS_POINT (v3)->z);
}

static void draw_sphere (void)
{
  glShadeModel (GL_SMOOTH);
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix ();
  glScalef (0.5, 0.5, 0.5);
  glBegin (GL_TRIANGLES);
  static GtsSurface * sphere = NULL;
  if (!sphere) {
    sphere = gts_surface_new (gts_surface_class (), gts_face_class (), 
			      gts_edge_class (), gts_vertex_class ());
    gts_surface_generate_sphere (sphere, 3);
  }
  gts_surface_foreach_face (sphere, (GtsFunc) draw_face, NULL);
  glEnd ();
  glPopMatrix ();
}

FILE * gfs_gl_popen (const gchar * fname)
{
  g_return_val_if_fail (fname != NULL, NULL);

  FILE * fptr = fopen (fname, "r");
  if (fptr == NULL)
    return NULL;
  fclose (fptr);

  gchar * script = g_strconcat ("if gunzip -t \"", fname, "\" 2> /dev/null; then "
				"  gunzip -c \"", fname, "\" 2> /dev/null; else ",
				"  cat \"", fname, "\"; fi", NULL);
  fptr = popen (script, "r");
  g_free (script);
  return fptr;
}

static GtsColor * color_new (gdouble r, gdouble g, gdouble b)
{
  GtsColor * c = g_malloc (sizeof (GtsColor));
  c->r = r; c->g = g; c->b = b;
  return c;
}

static void color_destroy (GtsColor * color)
{
  g_return_if_fail (color != NULL);

  g_free (color);
}

static void colormap_set_texture (GfsColormap * cmap)
{
  guint i;

  for (i = 0; i < GFS_COLORMAP_TEXTURE_SAMPLES; i++) {
    GtsColor c = gfs_colormap_color (cmap, i/(gdouble) (GFS_COLORMAP_TEXTURE_SAMPLES - 1));
    cmap->texture[3*i] = c.r;
    cmap->texture[3*i + 1] = c.g;
    cmap->texture[3*i + 2] = c.b;
  }
}

GfsColormap * gfs_colormap_jet (void)
{
  GfsColormap * cmap = g_malloc (sizeof (GfsColormap));
  gint i;

  cmap->reversed = FALSE;
  cmap->colors = g_ptr_array_new ();
  for (i = 0; i < 127; i++) {
    gdouble r = 
      i <= 46 ? 0. : 
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. : 
      0.03125*(i - 46);
    gdouble g = 
      i <= 14 || i >= 111 ? 0. : 
      i >= 79 ? -0.03125*(i - 111) : 
      i <= 46 ? 0.03125*(i - 14) : 
      1.;
    gdouble b =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;

    g_ptr_array_add (cmap->colors, color_new (r, g, b));
  }
  colormap_set_texture (cmap);
  cmap->name = g_strdup ("Jet");
  return cmap;
}

GfsColormap * gfs_colormap_cool_warm (void)
{
  GfsColormap * cmap = g_malloc (sizeof (GfsColormap));
  gint i;

  /* diverging cool-warm from:
   *  http://www.sandia.gov/~kmorel/documents/ColorMaps/CoolWarmFloat33.csv
   * see also:
   *  Diverging Color Maps for Scientific Visualization (Expanded)
   *  Kenneth Moreland
   */
  static double basemap[33][3] = {	
    {0.2298057,   0.298717966, 0.753683153},
    {0.26623388,  0.353094838, 0.801466763},
    {0.30386891,  0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334,  0.50941904,  0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708,  0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021,  0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803,  0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856,  0.387970225},
    {0.89904617,  0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379,  0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616,  0.150232812}	
  };
  
  cmap->reversed = FALSE;
  cmap->colors = g_ptr_array_new ();
  for (i = 0; i < 33; i++)
    g_ptr_array_add (cmap->colors, color_new (basemap[i][0], basemap[i][1], basemap[i][2]));
  colormap_set_texture (cmap);
  cmap->name = g_strdup ("Cool");
  return cmap;
}

GfsColormap * gfs_colormap_gray (void)
{
  GfsColormap * cmap = g_malloc (sizeof (GfsColormap));
  cmap->reversed = FALSE;
  cmap->colors = g_ptr_array_new ();
  g_ptr_array_add (cmap->colors, color_new (0,0,0));
  g_ptr_array_add (cmap->colors, color_new (1,1,1));
  colormap_set_texture (cmap);
  cmap->name = g_strdup ("Gray");
  return cmap;
}

void gfs_colormap_destroy (GfsColormap * colormap)
{
  guint i;

  g_return_if_fail (colormap != NULL);

  for (i = 0; i < colormap->colors->len; i++)
    color_destroy (colormap->colors->pdata[i]);
  g_ptr_array_free (colormap->colors, TRUE);
  g_free (colormap->name);
  g_free (colormap);
}

GtsColor gfs_colormap_color (GfsColormap * cmap, gdouble val)
{
  GtsColor c = {1., 1., 1.}, * c1, * c2;
  guint i, n;
  gdouble coef;

  g_return_val_if_fail (cmap != NULL, c);

  if (val > 1.0) val = 1.0;
  else if (val < 0.0) val = 0.0;
  if (cmap->reversed)
    val = 1.0 - val;

  n = cmap->colors->len;
  if (n == 0)
    return c;
  if (n == 1)
    return *((GtsColor *)cmap->colors->pdata[0]);

  i = floor ((gdouble)val*(gdouble)(n - 1));
  if (i == n - 1)
    return *((GtsColor *)cmap->colors->pdata[cmap->colors->len - 1]);
  coef = val*(gdouble)(n - 1) - (gdouble)i;
  c1 = cmap->colors->pdata[i];
  c2 = cmap->colors->pdata[i+1];
  c.r = c1->r + coef*(c2->r - c1->r);
  c.g = c1->g + coef*(c2->g - c1->g);
  c.b = c1->b + coef*(c2->b - c1->b);

  return c;
}

void gfs_colormap_texture (GfsColormap * cmap)
{
  g_return_if_fail (cmap != NULL);

  glTexImage1D (GL_TEXTURE_1D, 0, GL_RGB, GFS_COLORMAP_TEXTURE_SAMPLES, 0, GL_RGB, GL_FLOAT,
		cmap->texture);
  glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}

#define RC(r,c) m[(r)+(c)*4]
#define RCM(m,r,c) (m)[(r)+(c)*4]

static void matrix_multiply (float * m, float * n)
{
  float o[16];
  guint i;
  
  for (i = 0; i < 16; i++) o[i] = m[i];
  RC(0,0)=RCM(o,0,0)*RCM(n,0,0)+RCM(o,0,1)*RCM(n,1,0)+
          RCM(o,0,2)*RCM(n,2,0)+RCM(o,0,3)*RCM(n,3,0);
  RC(0,1)=RCM(o,0,0)*RCM(n,0,1)+RCM(o,0,1)*RCM(n,1,1)+
          RCM(o,0,2)*RCM(n,2,1)+RCM(o,0,3)*RCM(n,3,1);
  RC(0,2)=RCM(o,0,0)*RCM(n,0,2)+RCM(o,0,1)*RCM(n,1,2)+
          RCM(o,0,2)*RCM(n,2,2)+RCM(o,0,3)*RCM(n,3,2);
  RC(0,3)=RCM(o,0,0)*RCM(n,0,3)+RCM(o,0,1)*RCM(n,1,3)+
          RCM(o,0,2)*RCM(n,2,3)+RCM(o,0,3)*RCM(n,3,3);
  RC(1,0)=RCM(o,1,0)*RCM(n,0,0)+RCM(o,1,1)*RCM(n,1,0)+
          RCM(o,1,2)*RCM(n,2,0)+RCM(o,1,3)*RCM(n,3,0);
  RC(1,1)=RCM(o,1,0)*RCM(n,0,1)+RCM(o,1,1)*RCM(n,1,1)+
          RCM(o,1,2)*RCM(n,2,1)+RCM(o,1,3)*RCM(n,3,1);
  RC(1,2)=RCM(o,1,0)*RCM(n,0,2)+RCM(o,1,1)*RCM(n,1,2)+
          RCM(o,1,2)*RCM(n,2,2)+RCM(o,1,3)*RCM(n,3,2);
  RC(1,3)=RCM(o,1,0)*RCM(n,0,3)+RCM(o,1,1)*RCM(n,1,3)+
          RCM(o,1,2)*RCM(n,2,3)+RCM(o,1,3)*RCM(n,3,3);
  RC(2,0)=RCM(o,2,0)*RCM(n,0,0)+RCM(o,2,1)*RCM(n,1,0)+
          RCM(o,2,2)*RCM(n,2,0)+RCM(o,2,3)*RCM(n,3,0);
  RC(2,1)=RCM(o,2,0)*RCM(n,0,1)+RCM(o,2,1)*RCM(n,1,1)+
          RCM(o,2,2)*RCM(n,2,1)+RCM(o,2,3)*RCM(n,3,1);
  RC(2,2)=RCM(o,2,0)*RCM(n,0,2)+RCM(o,2,1)*RCM(n,1,2)+
          RCM(o,2,2)*RCM(n,2,2)+RCM(o,2,3)*RCM(n,3,2);
  RC(2,3)=RCM(o,2,0)*RCM(n,0,3)+RCM(o,2,1)*RCM(n,1,3)+
          RCM(o,2,2)*RCM(n,2,3)+RCM(o,2,3)*RCM(n,3,3);
  RC(3,0)=RCM(o,3,0)*RCM(n,0,0)+RCM(o,3,1)*RCM(n,1,0)+
          RCM(o,3,2)*RCM(n,2,0)+RCM(o,3,3)*RCM(n,3,0);
  RC(3,1)=RCM(o,3,0)*RCM(n,0,1)+RCM(o,3,1)*RCM(n,1,1)+
          RCM(o,3,2)*RCM(n,2,1)+RCM(o,3,3)*RCM(n,3,1);
  RC(3,2)=RCM(o,3,0)*RCM(n,0,2)+RCM(o,3,1)*RCM(n,1,2)+
          RCM(o,3,2)*RCM(n,2,2)+RCM(o,3,3)*RCM(n,3,2);
  RC(3,3)=RCM(o,3,0)*RCM(n,0,3)+RCM(o,3,1)*RCM(n,1,3)+
          RCM(o,3,2)*RCM(n,2,3)+RCM(o,3,3)*RCM(n,3,3);
}

static void vector_multiply (float * v, float * m)
{
  float o[4];
  guint i;
  
  for (i = 0; i < 4; i++) o[i] = v[i];
  
  v[0]=RC(0,0)*o[0]+RC(0,1)*o[1]+RC(0,2)*o[2]+RC(0,3)*o[3];
  v[1]=RC(1,0)*o[0]+RC(1,1)*o[1]+RC(1,2)*o[2]+RC(1,3)*o[3];
  v[2]=RC(2,0)*o[0]+RC(2,1)*o[1]+RC(2,2)*o[2]+RC(2,3)*o[3];
  v[3]=RC(3,0)*o[0]+RC(3,1)*o[1]+RC(3,2)*o[2]+RC(3,3)*o[3];
}

void gfs_gl_get_frustum (GfsGlViewParams * vp, GList * symmetries, GfsFrustum * f)
{
  GLint v[4];
  float p[16];
  int i;

  f->res = 2.*vp->res;
  f->symmetries = symmetries;
  guint n = 1;
  while (symmetries) {
    n *= 2;
    symmetries = symmetries->next;
  }
  f->s = g_malloc (n*sizeof (FttVector));

  glGetIntegerv (GL_VIEWPORT, v);
  f->width = v[2];
  glGetFloatv (GL_MODELVIEW_MATRIX, f->m);
  glGetFloatv (GL_PROJECTION_MATRIX, f->p);
  for (i = 0; i < 16; i++) p[i] = f->p[i];
  matrix_multiply (p, f->m);

  /* right */
  f->n[0][0] = p[3] - p[0];
  f->n[0][1] = p[7] - p[4];
  f->n[0][2] = p[11] - p[8];
  f->d[0]    = p[15] - p[12];
   
  /* left */
  f->n[1][0] = p[3] + p[0];
  f->n[1][1] = p[7] + p[4];
  f->n[1][2] = p[11] + p[8];
  f->d[1]    = p[15] + p[12];
  
  /* top */
  f->n[2][0] = p[3] - p[1];
  f->n[2][1] = p[7] - p[5];
  f->n[2][2] = p[11] - p[9];
  f->d[2]    = p[15] - p[13];

  /* bottom */
  f->n[3][0] = p[3] + p[1];
  f->n[3][1] = p[7] + p[5];
  f->n[3][2] = p[11] + p[9];
  f->d[3]    = p[15] + p[13];
  
  /* front */
  f->n[4][0] = p[3] + p[2];
  f->n[4][1] = p[7] + p[6];
  f->n[4][2] = p[11] + p[10];
  f->d[4]    = p[15] + p[14];
  
  /* back */
  f->n[5][0] = p[3] - p[2];
  f->n[5][1] = p[7] - p[6];
  f->n[5][2] = p[11] - p[10];
  f->d[5]    = p[15] - p[14];
  
  for (i = 0; i < 6; i++) {
    gdouble n = gts_vector_norm (f->n[i]);
    if (n > 0.) {
      f->n[i][0] /= n; f->n[i][1] /= n; f->n[i][2] /= n;
      f->d[i] /= n;
    }
  }
}

void gfs_gl_frustum_free (GfsFrustum * f)
{
  g_return_if_fail (f != NULL);
  g_free (f->s);
}

static guint create_symmetries (FttVector * f, GList * s, FttVector * p)
{
  guint j = 0, n;

  f[j++] = *p; n = j;
  while (s) {
    guint i;
    for (i = 0; i < n; i++)
      gfs_gl_symmetry_transform (s->data, &f[i], &f[j++]);
    n = j;
    s = s->next;
  }
  return n;
}

/**
 * gfs_sphere_in_frustum:
 * @p: the sphere center.
 * @r: the sphere radius.
 * @f: the view frustum.
 * 
 * Returns: GTS_OUT if the sphere is outside the view frustum, GTS_IN
 * if it is inside, GTS_ON if it is partly inside.
 */
GtsIntersect gfs_sphere_in_frustum (FttVector * p, gdouble r, GfsFrustum * f)
{
  GtsIntersect I1 = GTS_OUT;

  g_return_val_if_fail (p != NULL, GTS_OUT);
  g_return_val_if_fail (f != NULL, GTS_OUT);

  guint j, n = create_symmetries (f->s, f->symmetries, p);
  for (j = 0; j < n; j++) {
    p = &f->s[j];
    guint i;
    GtsIntersect I = GTS_IN;
    for (i = 0; i < 6; i++) {
      gdouble d = f->n[i][0]*p->x + f->n[i][1]*p->y + f->n[i][2]*p->z + f->d[i];
      if (d < -r) {
	I = GTS_OUT;
	break;
      }
      if (d < r)
	I = GTS_ON;
    }
    if (I == GTS_IN)
      return GTS_IN;
    if (I == GTS_ON)
      I1 = GTS_ON;
  }
  return I1;
}

/**
 * gfs_sphere_is_small:
 * @c: the sphere center.
 * @r: the sphere radius.
 * @f: the view frustum.
 * 
 * Returns: %TRUE if the screen size (in pixels) of the projected
 * sphere is smaller than the resolution of @f, %FALSE otherwise.
 */
gboolean gfs_sphere_is_small (FttVector * c, gdouble r, GfsFrustum * f)
{
  g_return_val_if_fail (c != NULL, FALSE);
  g_return_val_if_fail (f != NULL, FALSE);

  guint j, n = create_symmetries (f->s, f->symmetries, c);
  for (j = 0; j < n; j++) {
    c = &f->s[j];
    float v[4];
    v[0] = c->x; v[1] = c->y; v[2] = c->z; v[3] = 1.;
    vector_multiply (v, f->m);
    v[0] = r;
    vector_multiply (v, f->p);
    float rp = v[3] == 0. ? 0 : v[0]*f->width/v[3];
    if (rp >= f->res)
      return FALSE;
  }
  return TRUE;
}

static void cell_traverse_visible_no_check (FttCell * root,
					    GfsFrustum * f,
					    gint maxlevel,
					    FttCellTraverseFunc func,
					    gpointer data)
{
  if (FTT_CELL_IS_LEAF (root) || ftt_cell_level (root) == maxlevel)
    (* func) (root, data);
  else {
    gdouble r = ftt_cell_size (root)*GFS_DIAGONAL;
    FttVector p;
    
    ftt_cell_pos (root, &p);
    if (gfs_sphere_is_small (&p, r, f))
      (* func) (root, data);
    else {
      struct _FttOct * children = root->children;
      guint n;
      
      for (n = 0; n < FTT_CELLS; n++) {
	FttCell * c = &(children->cell[n]);
	if (!FTT_CELL_IS_DESTROYED (c))
	  cell_traverse_visible_no_check (c, f, maxlevel, func, data);
      }
    }
  }
}

static void cell_traverse_visible (FttCell * root,
				   GfsFrustum * f,
				   gint maxlevel,
				   FttCellTraverseFunc func,
				   gpointer data)
{
  gdouble r = ftt_cell_size (root)*GFS_DIAGONAL;
  FttVector p;
  GtsIntersect i;

  ftt_cell_pos (root, &p);
  i = gfs_sphere_in_frustum (&p, r, f);
  if (i == GTS_OUT)
    return;
  if (FTT_CELL_IS_LEAF (root) ||
      ftt_cell_level (root) == maxlevel || 
      gfs_sphere_is_small (&p, r, f))
    (* func) (root, data);
  else if (i == GTS_IN)
    cell_traverse_visible_no_check (root, f, maxlevel, func, data);
  else {
    struct _FttOct * children = root->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);
      if (!FTT_CELL_IS_DESTROYED (c))
	cell_traverse_visible (c, f, maxlevel, func, data);
    }
  }
}

static void box_traverse_visible (GfsBox * b, gpointer * data)
{
  cell_traverse_visible (b->root, data[0], *((gint *)data[3]), data[1], data[2]);
}

/**
 * gfs_gl_cell_traverse_visible:
 * @gl: a #GfsGl.
 * @f: a view frustum.
 * @func: a used-defined function.
 * @data: user data to pass to @func.
 *
 * Traverse the cells of @gl which are visible.
 */
void gfs_gl_cell_traverse_visible (GfsGl * gl,
				   GfsFrustum * f,
				   FttCellTraverseFunc func,
				   gpointer data)
{
  gpointer dat[4];

  g_return_if_fail (gl != NULL);
  g_return_if_fail (f != NULL);
  g_return_if_fail (func != NULL);

  dat[0] = f;
  dat[1] = func;
  dat[2] = data;
  dat[3] = &gl->maxlevel;
  gts_container_foreach (GTS_CONTAINER (gl->sim), (GtsFunc) box_traverse_visible, dat);
}

typedef struct {
  GfsFrustum * f;
  gint maxlevel;
  gboolean (* condition) (FttCell *, gpointer);
  gpointer datacon;
  FttCellTraverseFunc func;
  gpointer data;
} ConditionParams;

static void cell_traverse_visible_condition_no_check (FttCell * root, ConditionParams * par)
{
  if (!(* par->condition) (root, par->datacon))
    return;
  if (FTT_CELL_IS_LEAF (root) || ftt_cell_level (root) == par->maxlevel)
    (* par->func) (root, par->data);
  else {
    gdouble r = ftt_cell_size (root)*SLIGHTLY_LARGER*GFS_DIAGONAL;
    FttVector p;
    
    ftt_cell_pos (root, &p);
    if (gfs_sphere_is_small (&p, r, par->f))
      (* par->func) (root, par->data);
    else {
      struct _FttOct * children = root->children;
      guint n;
      
      for (n = 0; n < FTT_CELLS; n++) {
	FttCell * c = &(children->cell[n]);
	if (!FTT_CELL_IS_DESTROYED (c))
	  cell_traverse_visible_condition_no_check (c, par);
      }
    }
  }
}

static void cell_traverse_visible_condition (FttCell * root, ConditionParams * par)
{
  gdouble r = ftt_cell_size (root)*SLIGHTLY_LARGER*GFS_DIAGONAL;
  FttVector p;
  GtsIntersect i;

  if (!(* par->condition) (root, par->datacon))
    return;
  ftt_cell_pos (root, &p);
  i = gfs_sphere_in_frustum (&p, r, par->f);
  if (i == GTS_OUT)
    return;
  if (FTT_CELL_IS_LEAF (root) ||
      ftt_cell_level (root) == par->maxlevel || 
      gfs_sphere_is_small (&p, r, par->f))
    (* par->func) (root, par->data);
  else if (i == GTS_IN)
    cell_traverse_visible_condition_no_check (root, par);
  else {
    struct _FttOct * children = root->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);
      if (!FTT_CELL_IS_DESTROYED (c))
	cell_traverse_visible_condition (c, par);
    }
  }
}

static void box_traverse_visible_condition (GfsBox * b, ConditionParams * par)
{
  cell_traverse_visible_condition (b->root, par);
}

/**
 * gfs_gl_cell_traverse_visible_condition:
 * @gl: a #GfsGl.
 * @f: a view frustum.
 * @condition: the condition to verify.
 * @datacon: user data to pass to @condition.
 * @func: a used-defined function.
 * @data: user data to pass to @func.
 *
 * Traverse the cells of @gl which are visible and verify @condition.
 */
void gfs_gl_cell_traverse_visible_condition (GfsGl * gl,
					     GfsFrustum * f,
					     gboolean (* condition) (FttCell *, gpointer),
					     gpointer datacon,
					     FttCellTraverseFunc func,
					     gpointer data)
{
  ConditionParams par;

  g_return_if_fail (gl != NULL);
  g_return_if_fail (f != NULL);
  g_return_if_fail (condition != NULL);
  g_return_if_fail (func != NULL);

  par.f = f;
  par.condition = condition;
  par.datacon = datacon;
  par.func = func;
  par.data = data;
  par.maxlevel = gl->maxlevel;
  gts_container_foreach (GTS_CONTAINER (gl->sim), (GtsFunc) box_traverse_visible_condition, &par);
}

static gboolean is_mixed (FttCell * cell, gpointer data)
{
  return GFS_IS_MIXED (cell);
}

/**
 * gfs_gl_cell_traverse_visible_mixed:
 * @gl: a #GfsGl.
 * @f: a view frustum.
 * @func: a used-defined function.
 * @data: user data to pass to @func.
 *
 * Traverse the cells of @gl which are visible and mixed.
 */
void gfs_gl_cell_traverse_visible_mixed (GfsGl * gl,
					 GfsFrustum * f,
					 FttCellTraverseFunc func,
					 gpointer data)
{
  g_return_if_fail (gl != NULL);
  g_return_if_fail (f != NULL);
  g_return_if_fail (func != NULL);

  gfs_gl_cell_traverse_visible_condition (gl, f, is_mixed, NULL, func, data);
}

static void cell_traverse_visible_boundary_no_check (FttCell * root,
						     GfsFrustum * f,
						     FttDirection d,
						     gint maxlevel,
						     FttCellTraverseFunc func,
						     gpointer data)
{
  if (FTT_CELL_IS_LEAF (root) || ftt_cell_level (root) == maxlevel)
    (* func) (root, data);
  else {
    gdouble r = ftt_cell_size (root)*GFS_DIAGONAL;
    FttVector p;
    
    ftt_cell_pos (root, &p);
    if (gfs_sphere_is_small (&p, r, f))
      (* func) (root, data);
    else {
      FttCellChildren child;
      guint n;
      
      ftt_cell_children_direction (root, d, &child);
      for (n = 0; n < FTT_CELLS/2; n++)
	if (child.c[n])
	  cell_traverse_visible_boundary_no_check (child.c[n], f, d, maxlevel, func, data);
    }
  }
}

static void cell_traverse_visible_boundary (FttCell * root,
					    GfsFrustum * f,
					    FttDirection d,
					    gint maxlevel,
					    FttCellTraverseFunc func,
					    gpointer data)
{
  gdouble r = ftt_cell_size (root)*GFS_DIAGONAL;
  FttVector p;
  GtsIntersect i;

  ftt_cell_pos (root, &p);
  i = gfs_sphere_in_frustum (&p, r, f);
  if (i == GTS_OUT)
    return;
  if (FTT_CELL_IS_LEAF (root) ||
      ftt_cell_level (root) == maxlevel ||
      gfs_sphere_is_small (&p, r, f))
    (* func) (root, data);
  else if (i == GTS_IN)
    cell_traverse_visible_boundary_no_check (root, f, d, maxlevel, func, data);
  else {
    FttCellChildren child;
    guint n;

    ftt_cell_children_direction (root, d, &child);
    for (n = 0; n < FTT_CELLS/2; n++)
      if (child.c[n])
	cell_traverse_visible_boundary (child.c[n], f, d, maxlevel, func, data);
  }
}

static void box_traverse_visible_boundary (GfsBox * b, gpointer * data)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (!b->neighbor[d] || GFS_IS_BOUNDARY (b->neighbor[d]))
      cell_traverse_visible_boundary (b->root, data[0], d, *((gint *)data[3]), data[1], data[2]);
}

/**
 * gfs_gl_cell_traverse_visible_boundary:
 * @gl: a #GfsGl.
 * @f: a view frustum.
 * @func: a used-defined function.
 * @data: user data to pass to @func.
 *
 * Traverse the boundary cells of @gl which are visible.
 */
void gfs_gl_cell_traverse_visible_boundary (GfsGl * gl,
					    GfsFrustum * f,
					    FttCellTraverseFunc func,
					    gpointer data)
{
  gpointer dat[4];

  g_return_if_fail (gl != NULL);
  g_return_if_fail (f != NULL);
  g_return_if_fail (func != NULL);

  dat[0] = f;
  dat[1] = func;
  dat[2] = data;
  dat[3] = &gl->maxlevel;
  gts_container_foreach (GTS_CONTAINER (gl->sim), (GtsFunc) box_traverse_visible_boundary, dat);
}

void gfs_gl_init (void)
{
  gfs_gl_label_class ();
  gfs_gl_cells_class ();
  gfs_gl_fractions_class ();
  gfs_gl_boundaries_class ();
  gfs_gl_squares_class ();
  gfs_gl_linear_class ();
  gfs_gl_isoline_class ();
  gfs_gl_solid_class ();
  gfs_gl_vof_class ();
  gfs_gl_solid_class ();
  gfs_gl_levels_class ();
  gfs_gl_vectors_class ();
    gfs_gl_streamlines_class ();
  gfs_gl_ellipses_class ();
  gfs_gl_location_class ();
    gfs_gl_height_class ();
  gfs_gl_locate_class ();
  gfs_gl_pipes_class ();
  gfs_gl_clip_plane_class ();
  gfs_gl_clip_plane_class ();
  gfs_gl_cut_plane_class ();
#if (!FTT_2D)
  gfs_gl_isosurface_class ();
#endif /* 3D */
  gfs_gl_symmetry_class ();
    gfs_gl_periodic_class ();
}

void gfs_gl_init_gl (void)
{
  GLfloat light0_pos[4]   = { 0.0, 0.0, 50.0, 0.0 };
  GLfloat light0_color[4] = { 1., 1., 1., 1.0 }; /* white light */

  glDisable (GL_CULL_FACE);
  glEnable (GL_DEPTH_TEST);
  glEnable (GL_NORMALIZE);

  /* speedups */
  glEnable (GL_DITHER);
  glShadeModel (GL_SMOOTH);
  glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);
  glHint (GL_POLYGON_SMOOTH_HINT, GL_FASTEST);

  /* light */
  glLightfv (GL_LIGHT0, GL_POSITION, light0_pos);
  glLightfv (GL_LIGHT0, GL_DIFFUSE,  light0_color);
  glEnable (GL_LIGHT0);
  glEnable (GL_LIGHTING);

  glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable (GL_COLOR_MATERIAL);

  /* fonts */
#if HAVE_FTGL
  if (default_font)
    ftglDestroyFont (default_font);
  if (default_raster_font)
    ftglDestroyFont (default_raster_font);
  default_font = default_raster_font = NULL;
#endif /* HAVE_FTGL */
}

/* GfsGl: Object */

static void gl_read (GtsObject ** o, GtsFile * fp)
{
  GfsGl * gl = GFS_GL (*o);
  gchar * shading = NULL;
  GtsFileVariable var[] = {
    {GTS_FLOAT,  "r",           TRUE, &gl->lc.r},
    {GTS_FLOAT,  "g",           TRUE, &gl->lc.g},
    {GTS_FLOAT,  "b",           TRUE, &gl->lc.b},
    {GTS_STRING, "shading",     TRUE, &shading},
    {GTS_INT,    "maxlevel",    TRUE, &gl->maxlevel},
    {GTS_FLOAT,  "font_size",   TRUE, &gl->font_size},
    {GTS_INT,    "raster_font", TRUE, &gl->use_raster_font},
    {GTS_FLOAT,  "line_width",  TRUE, &gl->line_width},
    {GTS_NONE}
  };

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (class)");
    return;
  }
  gts_file_next_token (fp);

  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR) {
    g_free (shading);
    return;
  }

  if (var[3].set) {
    if (!strcmp (shading, "Constant"))
      gl->shading = GFS_GL_CONSTANT;
    else if (!strcmp (shading, "Flat"))
      gl->shading = GFS_GL_FLAT;
    else if (!strcmp (shading, "Smooth"))
      gl->shading = GFS_GL_SMOOTH;
    else if (!strcmp (shading, "CSmooth"))
      gl->shading = GFS_GL_CSMOOTH;
    else {
      gts_file_variable_error (fp, var, "shading", "unknown shading `%s'", shading);
      g_free (shading);
      return;
    }
    g_free (shading);
  }

  gfs_gl_set_raster_font (gl, gl->use_raster_font);
}

static void gl_write (GtsObject * o, FILE * fp)
{
  GfsGl * gl = GFS_GL (o);

  g_assert (strlen (o->klass->info.name) > 5);
  fprintf (fp, "%s {\n"
	   "  r = %g g = %g b = %g\n"
	   "  shading = %s\n"
	   "  maxlevel = %d\n"
	   "  font_size = %g\n"
	   "  raster_font = %d\n"
	   "  line_width = %g\n"
	   "}",
	   &o->klass->info.name[5],
	   gl->lc.r, gl->lc.g, gl->lc.b,
	   gl->shading == GFS_GL_CONSTANT ? "Constant" :
	   gl->shading == GFS_GL_FLAT ?     "Flat" :
	   gl->shading == GFS_GL_SMOOTH ?   "Smooth" :
	   gl->shading == GFS_GL_CSMOOTH ?  "CSmooth" : 
	                                    "Unknown",
	   gl->maxlevel,
	   gl->font_size,
	   gl->use_raster_font,
	   gl->line_width);
}

static void gl_set_simulation (GfsGl * gl, GfsSimulation * sim)
{
  gl->sim = sim;
}

static gboolean gl_relevant (GfsSimulation * sim)
{
  return TRUE;
}

static void gl_class_init (GfsGlClass * klass)
{
  klass->set_simulation = gl_set_simulation;
  klass->relevant = gl_relevant;
  GTS_OBJECT_CLASS (klass)->read = gl_read;
  GTS_OBJECT_CLASS (klass)->write = gl_write;
}

static void gl_init (GfsGl * gl)
{
  GtsColor c = { 0., 0., 0. };

  gl->shading = GFS_GL_CONSTANT;
  gl->lc = c;
  gl->maxlevel = -1;
  gl->format = GFSGL_SCREEN;
  gl->font_size = 1.;
  gl->use_raster_font = 1;
  gl->line_width = 1.;
}

GfsGlClass * gfs_gl_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_info = {
      "GfsGl",
      sizeof (GfsGl),
      sizeof (GfsGlClass),
      (GtsObjectClassInitFunc) gl_class_init,
      (GtsObjectInitFunc) gl_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &gfs_gl_info);
  }

  return klass;
}

GfsGl * gfs_gl_new (GfsGlClass * klass)
{
  GfsGl * object;

  g_return_val_if_fail (klass != NULL, NULL);

  object = GFS_GL (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

/**
 * gfs_gl_new_from_file:
 * @fp: a #GtsFile pointer.
 *
 * Returns: a new #GfsGl object read from @fp.
 */
GfsGl * gfs_gl_new_from_file (GtsFile * fp)
{
  GtsObjectClass * klass;
  GfsGl * gl;
  GtsObject * o;

  g_return_val_if_fail (fp != NULL, NULL);

  klass = gts_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gchar * ename = g_strconcat ("GfsGl", fp->token->str, NULL);
    klass = gts_object_class_from_name (ename);
    g_free (ename);
  }
  if (klass == NULL || !gts_object_class_is_from_class (klass, gfs_gl_class ()))
    return NULL;
  gl = gfs_gl_new (GFS_GL_CLASS (klass));
  o = GTS_OBJECT (gl);
  (* klass->read) (&o, fp);
  if (fp->type == GTS_ERROR) {
    gts_object_destroy (o);
    return NULL;
  }
  else
    return gl;
}

void gfs_gl_set_simulation (GfsGl * gl, GfsSimulation * sim)
{
  g_return_if_fail (gl != NULL);
  g_return_if_fail (sim != NULL);

  (* GFS_GL_CLASS (GTS_OBJECT (gl)->klass)->set_simulation) (gl, sim);
}

void gfs_gl_draw (GfsGl * gl, GfsFrustum * f)
{
  g_return_if_fail (gl != NULL);

  if (gl->sim && GFS_GL_CLASS (GTS_OBJECT (gl)->klass)->draw) {
    glColor3f (gl->lc.r, gl->lc.g, gl->lc.b);
    gl2psLineWidth (gl->line_width*gl->p->lw);
    glLineWidth (gl->line_width);
    (* GFS_GL_CLASS (GTS_OBJECT (gl)->klass)->draw) (gl, f);
  }
}

void gfs_gl_set_raster_font (GfsGl * gl, gboolean raster)
{
  g_return_if_fail (gl != NULL);

  gl->use_raster_font = raster;
}

void gfs_gl_set_font_size (GfsGl * gl, gfloat size)
{
  g_return_if_fail (gl != NULL);

  gl->font_size = size;
}

/* GfsGlLabel: Object */

static void gl_label_destroy (GtsObject * o)
{
  GfsGlLabel * l = GFS_GL_LABEL (o);

  if (l->formatted_label != l->label)
    g_free (l->formatted_label);
  g_free (l->label);

  (* GTS_OBJECT_CLASS (gfs_gl_label_class ())->parent_class->destroy) (o);
}

static void gl_label_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_gl_label_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsGlLabel * gl = GFS_GL_LABEL (*o);
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x", TRUE, &gl->p.x },
    {GTS_DOUBLE, "y", TRUE, &gl->p.y },
    {GTS_DOUBLE, "z", TRUE, &gl->p.z },
    {GTS_STRING, "label", TRUE, &gl->label },
    {GTS_INT,  "symbol", TRUE, &gl->symbol },
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
}

static void gl_label_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_gl_label_class ())->parent_class->write) (o, fp);
  GfsGlLabel * gl = GFS_GL_LABEL (o);
  fprintf (fp, 
	   " {\n"
	   "  x = %g y = %g z = %g\n"
	   "  label = \"%s\"\n"
	   "  symbol = %d\n"
	   "}",
	   gl->p.x, gl->p.y, gl->p.z,
	   gl->label, gl->symbol);
}

static void gl_text_setup (GfsGl * gl)
{
#if HAVE_FTGL
  if (!default_font) {
    default_font = ftglCreatePolygonFont (DEFAULT_FONT);
    default_raster_font = ftglCreateTextureFont (DEFAULT_FONT);
    if (!default_font || !default_raster_font)
      g_warning ("cannot create default FTGL font: %s", DEFAULT_FONT);
    else {
      ftglSetFontFaceSize (default_font,        DEFAULT_FONT_SIZE, DEFAULT_FONT_SIZE);
      ftglSetFontFaceSize (default_raster_font, DEFAULT_FONT_SIZE, DEFAULT_FONT_SIZE);
      /* fixme: for some weird reason using display lists doesn't work for polygon fonts */
      ftglSetFontDisplayList (default_font, 0);
    }
  }
#endif /* HAVE_FTGL */
}

static void gl_draw_text (GfsGl * gl, const gchar * text, 
			  gdouble x, gdouble y, gdouble z, 
			  gdouble size)
{
#if HAVE_FTGL
  if (text) {
    gl_text_setup (gl);
    glMatrixMode (GL_MODELVIEW);
    glPushMatrix ();
    glNormal3f (0., 0., 1.);
    glTranslatef (x, y, z);
    size /= DEFAULT_FONT_SIZE;
    glScalef (size, size, size);
    ftglRenderFont (gfs_gl_vector_format (gl) || !gl->use_raster_font ? 
		    default_font : default_raster_font,
		    text,
		    FTGL_RENDER_ALL);
    glPopMatrix ();
  }
#endif /* HAVE_FTGL */
}

static void gl_text_bbox (GfsGl * gl, const gchar * text, 
			  gdouble size,
			  float bbox[6])
{
#if HAVE_FTGL
  if (text) {
    gl_text_setup (gl);
    ftglGetFontBBox (gfs_gl_vector_format (gl) || !gl->use_raster_font ? 
		     default_font : default_raster_font,
		     text, strlen (text), bbox);
    int i;
    for (i = 0; i < 6; i++)
      bbox[i] *= size/DEFAULT_FONT_SIZE;
  }
#endif /* HAVE_FTGL */
}

static void gl_label_draw (GfsGl * gl, GfsFrustum * f)
{
  GfsGlLabel * l = GFS_GL_LABEL (gl);
  FttVector p = l->p;
  gfs_simulation_map (gl->sim, &p);
  gdouble size;
  if (gl->maxlevel >= 0)
    size = ftt_level_size (gl->maxlevel);
  else {
    FttCell * cell = gfs_domain_locate (GFS_DOMAIN (gl->sim), p, -1, NULL);
    size = cell ? ftt_cell_size (cell) : 1./64.;
  }

  gl->size = 0;
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glTranslatef (0., 0., gl->p->lc);
  if (l->symbol) {
    glMatrixMode (GL_MODELVIEW);
    glPushMatrix ();
    glTranslated (p.x, p.y, p.z);
    glScaled (size, size, size);
    draw_sphere ();
    glPopMatrix ();
    p.x += size/2.; p.y += size/2.;
  }
  glNormal3f (0., 0., 1.);
  gl_draw_text (gl, l->formatted_label, p.x, p.y, p.z, gl->font_size*size);
  glMatrixMode (GL_PROJECTION);
  glPopMatrix ();
}

void gfs_gl_label_set_label (GfsGlLabel * gl, const gchar * label, GfsSimulation * sim)
{
  g_return_if_fail (gl != NULL);
  g_return_if_fail (label != NULL);
  g_return_if_fail (sim != NULL);

  if (label != gl->label) {
    if (gl->formatted_label != gl->label)
      g_free (gl->formatted_label);
    gl->formatted_label = NULL;
    g_free (gl->label);
    gl->label = g_strdup (label);
  }

  gboolean dynamic = FALSE, parallel = FALSE;
  GSList * formats = gfs_format_new (gl->label, NULL, &dynamic, &parallel);
  if (dynamic || parallel) {
    if (gl->formatted_label != gl->label)
      g_free (gl->formatted_label);
    gl->formatted_label = 
      gfs_format_string (formats, GFS_DOMAIN (sim)->pid, sim->time.i, sim->time.t);
  }
  else
    gl->formatted_label = gl->label;
  gfs_format_destroy (formats);  
}

static void gl_label_set_simulation (GfsGl * gl, GfsSimulation * sim)
{  
  (*GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_label_class ())->parent_class)->set_simulation) 
    (gl, sim);

  gfs_gl_label_set_label (GFS_GL_LABEL (gl), GFS_GL_LABEL (gl)->label, sim);
}

static void gl_label_class_init (GfsGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gl_label_destroy;
  GTS_OBJECT_CLASS (klass)->read = gl_label_read;
  GTS_OBJECT_CLASS (klass)->write = gl_label_write;
  klass->set_simulation = gl_label_set_simulation;
  klass->draw = gl_label_draw;
}

static void label_init (GfsGl * gl)
{
  gl->maxlevel = 4;
  gl->font_size = 3;
  GFS_GL_LABEL (gl)->label = g_strdup ("Label");
}

GfsGlClass * gfs_gl_label_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_label_info = {
      "GfsGlLabel",
      sizeof (GfsGl2D),
      sizeof (GfsGl2DClass),
      (GtsObjectClassInitFunc) gl_label_class_init,
      (GtsObjectInitFunc) label_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_class ()), &gfs_gl_label_info);
  }

  return klass;
}

/* GfsGl2D: Object */

static void gl2D_read (GtsObject ** o, GtsFile * fp)
{
  GfsGl2D * gl = GFS_GL2D (*o);
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "n.x", TRUE},
    {GTS_DOUBLE, "n.y", TRUE},
    {GTS_DOUBLE, "n.z", TRUE},
    {GTS_DOUBLE, "pos", TRUE},
    {GTS_NONE}
  };

  (* GTS_OBJECT_CLASS (gfs_gl2D_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  var[0].data = &gl->n.x;
  var[1].data = &gl->n.y;
  var[2].data = &gl->n.z;
  var[3].data = &gl->pos;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  gfs_gl2D_update_plane (gl);
}

static void gl2D_write (GtsObject * o, FILE * fp)
{
  GfsGl2D * gl = GFS_GL2D (o);

  (* GTS_OBJECT_CLASS (gfs_gl2D_class ())->parent_class->write) (o, fp);

  fprintf (fp, " {\n"
	   "  n.x = %g n.y = %g n.z = %g\n"
	   "  pos = %g\n"
	   "}",
	   gl->n.x, gl->n.y, gl->n.z,
	   gl->pos);
}

static void gl2D_update_plane (GfsGl2D * gl)
{
  gdouble n;
  GtsVector Q0 = {0., 0., 0.};
  GtsVector Q1 = {0., 0., 0.};
  gdouble max;
  guint d = 0;
  
  n = gts_vector_norm (&gl->n.x);
  g_assert (n > 0.);
  gts_vector_normalize (&gl->n.x);

  /* build a vector orthogonal to the constraint */
  max = gl->n.x*gl->n.x;
  if (gl->n.y*gl->n.y > max) { max = gl->n.y*gl->n.y; d = 1; }
  if (gl->n.z*gl->n.z > max) { max = gl->n.z*gl->n.z; d = 2; }
  switch (d) {
  case 0: Q0[0] = - gl->n.z/gl->n.x; Q0[2] = 1.0; break;
  case 1: Q0[1] = - gl->n.z/gl->n.y; Q0[2] = 1.0; break;
  case 2: Q0[2] = - gl->n.x/gl->n.z; Q0[0] = 1.0; break;
  }
  
  /* build a second vector orthogonal to the first and to the constraint */
  gts_vector_cross (Q1, &gl->n.x, Q0);
  
  gl->p[0].x = gl->pos*gl->n.x;
  gl->p[0].y = gl->pos*gl->n.y;
  gl->p[0].z = gl->pos*gl->n.z;
  
  gl->p[1].x = gl->p[0].x + Q0[0];
  gl->p[1].y = gl->p[0].y + Q0[1];
  gl->p[1].z = gl->p[0].z + Q0[2];
  
  gl->p[2].x = gl->p[0].x + Q1[0];
  gl->p[2].y = gl->p[0].y + Q1[1];
  gl->p[2].z = gl->p[0].z + Q1[2];
}

static gdouble gl2D_pick (GfsGl * gl, GfsGlRay * r)
{
  GfsGl2D * gl2 = GFS_GL2D (gl);
  GtsVector AB;
  gdouble ABn;

  gts_vector_init (AB, &r->a, &r->b);
  ABn = gts_vector_scalar (AB, &gl2->n.x);
  if (fabs (ABn) < 1e-6) {
    gl2->picked = NULL;
    return GFS_NODATA;
  }
  else {
    GtsVector AP;
    gdouble a;

    AP[0] = gl2->n.x*gl2->pos - r->a.x;
    AP[1] = gl2->n.y*gl2->pos - r->a.y;
    AP[2] = gl2->n.z*gl2->pos - r->a.z;
    a = gts_vector_scalar (AP, &gl2->n.x)/ABn;
    gl2->pickedpos.x = r->a.x*(1. - a) + a*r->b.x;
    gl2->pickedpos.y = r->a.y*(1. - a) + a*r->b.y;
    gl2->pickedpos.z = r->a.z*(1. - a) + a*r->b.z;
    gl2->picked = gfs_domain_locate (GFS_DOMAIN (gl->sim), gl2->pickedpos, gl->maxlevel, NULL);
    return gl2->picked ? a : GFS_NODATA;
  }
}

static void gl2D_class_init (GfsGlClass * klass)
{
  klass->pick = gl2D_pick;
  GTS_OBJECT_CLASS (klass)->read = gl2D_read;
  GTS_OBJECT_CLASS (klass)->write = gl2D_write;
  GFS_GL2D_CLASS (klass)->update_plane = gl2D_update_plane;
}

static void gl2D_init (GfsGl2D * object)
{
  object->n.x = 0.; object->n.y = 0.; object->n.z = 1.;
  object->pos = 0.;

  object->p[0].x = object->p[0].y = object->p[0].z = 0.;
  object->p[1].x = 1.; object->p[1].y = object->p[1].z = 0.;
  object->p[2].y = 1.; object->p[2].x = object->p[2].z = 0.;

  gfs_gl2D_update_plane (object);
}

GfsGl2DClass * gfs_gl2D_class (void)
{
  static GfsGl2DClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl2D_info = {
      "GfsGl2D",
      sizeof (GfsGl2D),
      sizeof (GfsGl2DClass),
      (GtsObjectClassInitFunc) gl2D_class_init,
      (GtsObjectInitFunc) gl2D_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_class ()), &gfs_gl2D_info);
  }

  return klass;
}

void gfs_gl2D_update_plane (GfsGl2D * gl)
{
  g_return_if_fail (gl != NULL);
  
  (* GFS_GL2D_CLASS (GTS_OBJECT (gl)->klass)->update_plane) (gl);
}

/* GfsGlSymmetry: Object */

void gfs_gl_symmetry_transform (GfsGl * gl, FttVector * p, FttVector * t)
{
  g_return_if_fail (gl != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (t != NULL);

  GLfloat in[4];
  in[0] = p->x; in[1] = p->y; in[2] = p->z; in[3] = 1.;
  vector_multiply (in, GFS_GL_SYMMETRY (gl)->m);
  t->x = in[0]; t->y = in[1]; t->z = in[2];
}

void gfs_gl_symmetry_apply (GList * symmetry, GLuint display_list)
{
  if (symmetry) {
    guint length = g_list_length (symmetry);
    GLuint symmetry_list = glGenLists (length);
    GLuint list = symmetry_list;
    GList * i = symmetry;
    while (i) {
      glNewList (list, GL_COMPILE);
      glCallList (display_list);
      glMatrixMode (GL_MODELVIEW);
      glPushMatrix ();
      glMultMatrixf (GFS_GL_SYMMETRY (i->data)->m);
      glCallList (display_list);
      glPopMatrix ();
      glEndList ();
      display_list = list++;
      i = i->next;
    }
    glCallList (display_list);
    glDeleteLists (symmetry_list, length);
  }
  else
    glCallList (display_list);
}

static void gl_symmetry_update (GfsGl2D * gl)
{
  GLfloat * s = GFS_GL_SYMMETRY (gl)->m, t[16];

  gts_vector_normalize (&gl->n.x);
  if (GFS_GL_SYMMETRY (gl)->periodic) {
    s[0] = 1.;                    s[1] = 0.;                      s[2] = 0.;
    s[3] = 0.;
    s[4] = s[1];                  s[5] = 1.; s[6] = 0.;       
    s[7] = 0.;
    s[8] = s[2];                  s[9] = s[6];                    s[10] = 1.; 
    s[11] = 0.;
    s[12] = 0.;                   s[13] = 0.;                     s[14] = 0.;                    
    s[15] = 1.;

    t[0] = 1.;  t[1] = 0.;   t[2] = 0.;  t[3] = 0.;
    t[4] = 0.;  t[5] = 1.;   t[6] = 0.;  t[7] = 0.;
    t[8] = 0.;  t[9] = 0.;   t[10] = 1.; t[11] = 0.;
    t[12] = gl->n.x*gl->pos; 
    t[13] = gl->n.y*gl->pos;  
    t[14] = gl->n.z*gl->pos; 
    t[15] = 1.;
  }
  else {
    s[0] = 1. - 2.*gl->n.x*gl->n.x; s[1] = - 2.*gl->n.x*gl->n.y;  s[2] = - 2.*gl->n.x*gl->n.z;
    s[3] = 0.;
    s[4] = s[1];                  s[5] = 1. - 2.*gl->n.y*gl->n.y; s[6] = - 2.*gl->n.y*gl->n.z;
    s[7] = 0.;
    s[8] = s[2];                  s[9] = s[6];                    s[10] = 1. - 2.*gl->n.z*gl->n.z; 
    s[11] = 0.;
    s[12] = 0.;                   s[13] = 0.;                     s[14] = 0.;                    
    s[15] = 1.;

    t[0] = 1.;  t[1] = 0.;   t[2] = 0.;  t[3] = 0.;
    t[4] = 0.;  t[5] = 1.;   t[6] = 0.;  t[7] = 0.;
    t[8] = 0.;  t[9] = 0.;   t[10] = 1.; t[11] = 0.;
    t[12] = - 2.*gl->n.x*gl->pos; 
    t[13] = - 2.*gl->n.y*gl->pos;  
    t[14] = - 2.*gl->n.z*gl->pos; 
    t[15] = 1.;
  }
  matrix_multiply (s, t);
}

static void gl_symmetry_class_init (GfsGl2DClass * klass)
{
  GFS_GL_CLASS (klass)->pick = NULL;
  klass->update_plane = gl_symmetry_update;
}

GfsGlClass * gfs_gl_symmetry_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_symmetry_info = {
      "GfsGlSymmetry",
      sizeof (GfsGlSymmetry),
      sizeof (GfsGl2DClass),
      (GtsObjectClassInitFunc) gl_symmetry_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl2D_class ()), &gfs_gl_symmetry_info);
  }

  return klass;
}

/* GfsGlPeriodic: Object */

static void gl_periodic_init (GfsGlSymmetry * gl)
{
  gl->periodic = TRUE;
}

GfsGlClass * gfs_gl_periodic_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_periodic_info = {
      "GfsGlPeriodic",
      sizeof (GfsGlSymmetry),
      sizeof (GfsGl2DClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gl_periodic_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_symmetry_class ()), 
				  &gfs_gl_periodic_info);
  }

  return klass;
}

/* GfsGlCells: Object */

static void gl_cells_draw (GfsGl * gl, GfsFrustum * f)
{
  gl->size = 0;
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glTranslatef (0., 0., gl->p->lc);
  gfs_gl_normal (gl);
  gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc) gl_cell, gl);
  glPopMatrix ();
}

static void gl_cells_class_init (GfsGlClass * klass)
{
  klass->draw = gl_cells_draw;
}

GfsGlClass * gfs_gl_cells_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_cells_info = {
      "GfsGlCells",
      sizeof (GfsGl2D),
      sizeof (GfsGl2DClass),
      (GtsObjectClassInitFunc) gl_cells_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl2D_class ()), &gfs_gl_cells_info);
  }

  return klass;
}

/* GfsGlFractions: Object */

static gboolean gl_fractions_relevant (GfsSimulation * sim)
{
  return (sim->solids->items != NULL);
}

static void gl_fractions_class_init (GfsGlClass * klass)
{
  klass->draw = gl_fractions_draw;
  klass->relevant = gl_fractions_relevant;
}

static void gl_fractions_init (GfsGl * gl)
{
  GtsColor c = { 0., 0., 1.};

  gl->lc = c;
}

GfsGlClass * gfs_gl_fractions_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_fractions_info = {
      "GfsGlFractions",
      sizeof (GfsGl),
      sizeof (GfsGlClass),
      (GtsObjectClassInitFunc) gl_fractions_class_init,
      (GtsObjectInitFunc) gl_fractions_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_class ()), &gfs_gl_fractions_info);
  }

  return klass;
}

/* GfsGlBoundaries: Object */

static void gl_boundaries_draw (GfsGl * gl, GfsFrustum * f)
{
  gl->size = 0;
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glTranslatef (0., 0., gl->p->lc);
  glBegin (GL_LINES);
  gfs_gl_cell_traverse_visible_boundary (gl, f, (FttCellTraverseFunc) gl_boundaries, gl);
  glEnd ();
  glPopMatrix ();
}

static void gl_boundaries_class_init (GfsGlClass * klass)
{
  klass->draw = gl_boundaries_draw;
}

static void gl_boundaries_init (GfsGl * gl)
{
  GtsColor c = { 0., 0., 0.};

  gl->lc = c;
}

GfsGlClass * gfs_gl_boundaries_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_boundaries_info = {
      "GfsGlBoundaries",
      sizeof (GfsGl),
      sizeof (GfsGlClass),
      (GtsObjectClassInitFunc) gl_boundaries_class_init,
      (GtsObjectInitFunc) gl_boundaries_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_class ()), &gfs_gl_boundaries_info);
  }

  return klass;
}

/* GfsGlLevels: Object */

static void gl_levels_destroy (GtsObject * object)
{
  GfsGlLevels * gl = GFS_GL_LEVELS (object);

  if (gl->v)
    gts_object_destroy (GTS_OBJECT (gl->v));

  (* GTS_OBJECT_CLASS (gfs_gl_levels_class ())->parent_class->destroy) (object);
}

static void set_level (FttCell * cell, GfsVariable * v)
{
  GFS_VALUE (cell, v) = ftt_cell_level (cell);
}

static void gl_levels_set_simulation (GfsGl * object, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsGlLevels * gl = GFS_GL_LEVELS (object);

  (*GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_levels_class ())->parent_class)->set_simulation)
    (object, sim);

  if (gl->v)
    gts_object_destroy (GTS_OBJECT (gl->v));
  gl->v = gfs_temporary_variable (domain);

  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) set_level, gl->v);
  gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_get_from_below_intensive, gl->v);
  gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, gl->v);
}

static void gl_levels_draw (GfsGl * gl, GfsFrustum * f)
{
  gl->size = 0;
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glTranslatef (0., 0., gl->p->lc);
  glBegin (GL_LINES);
  gfs_gl_normal (gl);
  gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc) gl_face, gl);
  glEnd ();
  glPopMatrix ();
}

static void min_max_level (FttCell * cell, gint * l)
{
  guint level = ftt_cell_level (cell);
  if (level < l[0]) l[0] = level;
  if (level > l[1]) l[1] = level;
}

static gboolean gl_levels_relevant (GfsSimulation * sim)
{
  gint l[2] = { G_MAXINT, 0 };

  gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) min_max_level, l);
  return (l[1] > l[0]);
}

static void gl_levels_class_init (GfsGlClass * klass)
{
  klass->set_simulation = gl_levels_set_simulation;
  klass->draw = gl_levels_draw;
  klass->pick = NULL;
  klass->relevant = gl_levels_relevant;
  GTS_OBJECT_CLASS (klass)->destroy = gl_levels_destroy;
}

GfsGlClass * gfs_gl_levels_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_levels_info = {
      "GfsGlLevels",
      sizeof (GfsGlLevels),
      sizeof (GfsGl2DClass),
      (GtsObjectClassInitFunc) gl_levels_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl2D_class ()),
				  &gfs_gl_levels_info);
  }

  return klass;
}

/* GfsGlVarFunc: object */

GfsGlVarFunc * gfs_gl_var_func_new (void)
{
  GfsGlVarFunc * vf = g_malloc0 (sizeof (GfsGlVarFunc));
  return vf;
}

void gfs_gl_var_func_destroy (GfsGlVarFunc * vf)
{
  g_return_if_fail (vf != NULL);

  if (vf->v && vf->v != gfs_function_get_variable (vf->f))
    gts_object_destroy (GTS_OBJECT (vf->v));
  if (vf->f)
    gts_object_destroy (GTS_OBJECT (vf->f));
  g_free (vf);
}

static void update_v (FttCell * cell, gpointer * data)
{
  GFS_VALUE (cell, GFS_VARIABLE (data[1])) = gfs_function_value (GFS_FUNCTION (data[0]), cell);
}

GtsFile * gfs_gl_var_func_set (GfsGlVarFunc * vf, 
			       GfsSimulation * sim, 
			       const gchar * func,
			       GString * expr,
			       GfsVariableClass * klass)
{
  GfsDomain * domain;
  GfsFunction * f;
  GtsFile * fp;

  g_return_val_if_fail (vf != NULL, NULL);
  g_return_val_if_fail (sim != NULL, NULL);
  g_return_val_if_fail (func != NULL, NULL);
  if (!klass)
    klass = gfs_variable_class ();
  domain = GFS_DOMAIN (sim);

  fp = gts_file_new_from_string (func);
  f = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_read (f, domain, fp);
  gfs_pending_functions_compilation (fp);
  if (fp->type == GTS_ERROR) {
    gts_object_destroy (GTS_OBJECT (f));
    return fp;
  }
  gts_file_destroy (fp);
  
  GfsVariable * v;
  if (!(v = gfs_function_get_variable (f)) ||
      !gts_object_class_is_from_class (GTS_OBJECT (v)->klass, klass) ||
      (gfs_variable_is_dimensional (v) && sim->physical_params.L != 1.)) {
    v = gfs_variable_new (klass, domain, NULL, NULL);
    gfs_catch_floating_point_exceptions ();
    gpointer data[2] = { f, v };
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) update_v, data);
    if (gfs_restore_floating_point_exceptions ()) {
      fp = gts_file_new_from_string (func);
      gts_file_error (fp, "Floating-point exception");
      gts_object_destroy (GTS_OBJECT (v));
      gts_object_destroy (GTS_OBJECT (f));
      return fp;
    }
    gfs_event_init (GFS_EVENT (v), sim);
    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) v->fine_coarse, v);
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, v);
  }

  if (vf->v && vf->v != gfs_function_get_variable (vf->f))
    gts_object_destroy (GTS_OBJECT (vf->v));
  if (vf->f)
    gts_object_destroy (GTS_OBJECT (vf->f));
  vf->f = f;
  vf->v = v;

  if (expr && expr->str != func) {
    g_free (expr->str);
    expr->str = g_strdup (func);
    expr->len = strlen (expr->str);
  }
  return NULL;
}

/* GfsGlScalar: Object */

static void gl_scalar_read (GtsObject ** o, GtsFile * fp)
{
  GfsGlScalar * gl = GFS_GL_SCALAR (*o);
  gchar * cmap;
  GtsFileVariable var[] = {
    {GTS_INT,    "amin", TRUE, &gl->amin},
    {GTS_DOUBLE, "min",  TRUE, &gl->min},
    {GTS_INT,    "amax", TRUE, &gl->amax},
    {GTS_DOUBLE, "max",  TRUE, &gl->max},
    {GTS_STRING, "cmap", TRUE, &cmap},
    {GTS_INT,    "show", TRUE, &gl->show},
    {GTS_NONE}
  };

  (* GTS_OBJECT_CLASS (gfs_gl_scalar_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  g_string_free (gl->expr, TRUE);
  if (!(gl->expr = gfs_function_expression (fp, NULL)))
    return;
  gts_file_next_token (fp);

  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (var[4].set) {
    gfs_colormap_destroy (gl->cmap);
    if (!strcmp (cmap, "Cool"))
      gl->cmap = gfs_colormap_cool_warm ();
    else if (!strcmp (cmap, "Gray"))
      gl->cmap = gfs_colormap_gray ();
    else if (!strcmp (cmap, "Jet"))
      gl->cmap = gfs_colormap_jet ();
    else
      gts_file_error (fp, "unknown colormap '%s'", cmap);
    g_free (cmap);
  }
}

static void gl_scalar_write (GtsObject * o, FILE * fp)
{
  GfsGlScalar * gl = GFS_GL_SCALAR (o);
  
  (* GTS_OBJECT_CLASS (gfs_gl_scalar_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s {\n  amin = %d", gl->expr->str, gl->amin);
  if (gl->amin)
    fputc ('\n', fp);
  else
    fprintf (fp, " min = %g\n", gl->min);
  fprintf (fp, "  amax = %d", gl->amax);
  if (gl->amax)
    fputc ('\n', fp);
  else
    fprintf (fp, " max = %g\n", gl->max);
  if (gl->show)
    fputs ("  show = 1\n", fp);
  fprintf (fp, "  cmap = %s\n}", gl->cmap->name);
}

static void gl_scalar_destroy (GtsObject * object)
{
  GfsGlScalar * gl = GFS_GL_SCALAR (object);

  gfs_gl_var_func_destroy (gl->vf);
  g_string_free (gl->expr, TRUE);
  if (gl->cmap)
    gfs_colormap_destroy (gl->cmap);

  (* GTS_OBJECT_CLASS (gfs_gl_scalar_class ())->parent_class->destroy) (object);
}

static void min_max (FttCell * cell, GfsGlScalar * gl)
{
  gdouble v = GFS_VALUE (cell, gl->v);
  if (v < G_MAXDOUBLE && v > gl->amaxv)
    gl->amaxv = v;
  if (v < gl->aminv)
    gl->aminv = v;
}

static GtsFile * gl_scalar_set (GfsGlScalar * gl, const gchar * func)
{
  GtsFile * fp;

  if ((fp = gfs_gl_var_func_set (gl->vf, GFS_GL (gl)->sim, func, gl->expr, NULL)))
    return fp;

  gl->v = gl->vf->v;
  gl->amaxv = -G_MAXDOUBLE;
  gl->aminv =  G_MAXDOUBLE;
  gfs_domain_cell_traverse (GFS_DOMAIN (GFS_GL (gl)->sim),
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) min_max, gl);
  gfs_all_reduce (GFS_DOMAIN (GFS_GL (gl)->sim), gl->amaxv, MPI_DOUBLE, MPI_MAX);
  gfs_all_reduce (GFS_DOMAIN (GFS_GL (gl)->sim), gl->aminv, MPI_DOUBLE, MPI_MIN);
  if (gl->amax) gl->max = gl->amaxv;
  if (gl->amin) gl->min = gl->aminv;

  return NULL;
}

static void gl_scalar_set_simulation (GfsGl * object, GfsSimulation * sim)
{
  GfsGlScalar * gls = GFS_GL_SCALAR (object);
  GtsFile * fp = NULL;

  (*GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_scalar_class ())->parent_class)->set_simulation)
    (object, sim);

  if (gls->expr->str[0] == '\0' || (fp = gfs_gl_scalar_set (gls, gls->expr->str))) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    
    if (domain->variables) {
      GfsVariable * v = domain->variables->data;
      gfs_gl_scalar_set (gls, v->name);
    }
    else
      gfs_gl_scalar_set (gls, "0");
  }
  if (fp)
    gts_file_destroy (fp);
}

static gchar * gl_scalar_properties (GfsGlScalar * gl, int len)
{
  gchar * s = g_strndup (gl->expr->str, len + 3);
  if (strlen (s) > len) {
    guint i;
    for (i = 0; i < 3; i++)
      s[len + i] = '.';
  }
  return s;
}

static void gl_scalar_draw (GfsGl * gl, GfsFrustum * f)
{
  GfsGlScalar * gls = GFS_GL_SCALAR (gl);

  /* draw colorscale */
  if (gls->show && gls->max > gls->min) {
    glMatrixMode (GL_PROJECTION);
    glPushMatrix ();
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    int width = viewport[2], height = viewport[3];
    glLoadIdentity ();
    gluOrtho2D (0, width, 0, height);
    glMatrixMode (GL_MODELVIEW);
    glPushMatrix ();
    glLoadIdentity ();

    int i, np = 10;
    float border = gl->font_size*height/20.;
    float cwidth = border;
    float cheight = (gl->font_size*height - 3.*border)/np;

    gchar s[80];
    snprintf (s, 80, "% 6.3e", M_PI);
    float bbox[6];
    gl_text_bbox (gl, s, cheight/2., bbox);
    float twidth  = bbox[3] - bbox[0];
    float theight = bbox[4] - bbox[1];

    glTranslatef (width - border/2. - cwidth - border/4. - twidth, border, 0.);

    glColor3f (gl->lc.r, gl->lc.g, gl->lc.b);
    gchar * title = gl_scalar_properties (gls, strlen (s) - 3.);
    gl_text_bbox (gl, title, cheight/2., bbox);
    float Twidth  = bbox[3] - bbox[0];
    gl_draw_text (gl, title, 
		  cwidth + border/4. + (twidth - Twidth)/2., 
		  np*cheight + cheight/10., 0., cheight/2.);
    g_free (title);
    
    for (i = 0; i < np; i++) {
      GtsColor c = gfs_colormap_color (gls->cmap, (0.5 + i)/(float)np);
      glColor3f (c.r, c.g, c.b);
      glRectf (0, i*cheight, cwidth, (i + 1)*cheight);
      snprintf (s, 80, "% 6.3e", gls->min + (0.5 + i)*(gls->max - gls->min)/(float)np);
      glColor3f (gl->lc.r, gl->lc.g, gl->lc.b);
      gl_draw_text (gl, s, cwidth + border/4., i*cheight + (cheight - theight)/2., 0., cheight/2.);
    }

    glPopMatrix ();
    glMatrixMode (GL_PROJECTION);
    glPopMatrix ();
  }
}

static void gl_scalar_class_init (GfsGlClass * klass)
{
  klass->set_simulation = gl_scalar_set_simulation;
  klass->draw = gl_scalar_draw;
  GTS_OBJECT_CLASS (klass)->read = gl_scalar_read;
  GTS_OBJECT_CLASS (klass)->write = gl_scalar_write;
  GTS_OBJECT_CLASS (klass)->destroy = gl_scalar_destroy;
  GFS_GL_SCALAR_CLASS (klass)->set_scalar = gl_scalar_set;
}

static void gl_scalar_init (GfsGlScalar * object)
{
  object->expr = g_string_new ("");
  object->vf = gfs_gl_var_func_new ();
  object->amax = object->amin = TRUE;
  object->min = 0.;
  object->max = 1.;
  object->cmap = gfs_colormap_jet ();
}

GfsGlScalarClass * gfs_gl_scalar_class (void)
{
  static GfsGlScalarClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_scalar_info = {
      "GfsGlScalar",
      sizeof (GfsGlScalar),
      sizeof (GfsGlScalarClass),
      (GtsObjectClassInitFunc) gl_scalar_class_init,
      (GtsObjectInitFunc) gl_scalar_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl2D_class ()),
				  &gfs_gl_scalar_info);
  }

  return klass;
}

GtsFile * gfs_gl_scalar_set (GfsGlScalar * gl, const gchar * func)
{
  g_return_val_if_fail (gl != NULL, NULL);
  g_return_val_if_fail (func != NULL, NULL);

  return (* GFS_GL_SCALAR_CLASS (GTS_OBJECT (gl)->klass)->set_scalar) (gl, func);
}

/* GfsGlSquares: Object */

static void gl_squares_class_init (GfsGlClass * klass)
{
  klass->draw = gl_squares_draw;
}

GfsGlClass * gfs_gl_squares_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_squares_info = {
      "GfsGlSquares",
      sizeof (GfsGlSquares),
      sizeof (GfsGlScalarClass),
      (GtsObjectClassInitFunc) gl_squares_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_scalar_class ()),
				  &gfs_gl_squares_info);
  }

  return klass;
}

/* GfsGlLinear: Object */

static void gl_linear_destroy (GtsObject * o)
{
  GfsGlLinear * gl = GFS_GL_LINEAR (o);

  gfs_gl_var_func_destroy (gl->vf);
  g_string_free (gl->expr, TRUE);
  if (gl->nx) {
    gts_object_destroy (GTS_OBJECT (gl->nx));
    gts_object_destroy (GTS_OBJECT (gl->ny));
  }

  (* GTS_OBJECT_CLASS (gfs_gl_linear_class ())->parent_class->destroy) (o);
}

static void gl_linear_read (GtsObject ** o, GtsFile * fp)
{
  GfsGlLinear * gl = GFS_GL_LINEAR (*o);

  (* GTS_OBJECT_CLASS (gfs_gl_linear_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '{') { /* for backward compatibility */
    GtsFileVariable var[] = {
      {GTS_INT,    "reversed",   TRUE},
      {GTS_INT,    "use_scalar", TRUE},
      {GTS_NONE}
    };

    g_string_free (gl->expr, TRUE);
    if (!(gl->expr = gfs_function_expression (fp, NULL)))
      return;
    gts_file_next_token (fp);
    
    var[0].data = &gl->reversed;
    var[1].data = &gl->use_scalar;
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
  }
  else { /* this is the old format for linear */
    g_warning ("obsolete GfsGlLinear/GfsGlIsoline syntax");
    if (!GFS_IS_GL_ISOLINE (gl)) {
      gdouble scale;
      GtsFileVariable var[] = {
	{GTS_DOUBLE, "scale", TRUE},
	{GTS_NONE}
      };
      var[0].data = &scale;
      gts_file_assign_variables (fp, var);
      if (fp->type == GTS_ERROR)
	return;
    }
    gl->use_scalar = (GfsVariable *) TRUE;
  }
}

static void gl_linear_write (GtsObject * o, FILE * fp)
{
  GfsGlLinear * gl = GFS_GL_LINEAR (o);
  
  (* GTS_OBJECT_CLASS (gfs_gl_linear_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s {\n"
	   "  reversed = %d\n"
	   "  use_scalar = %d\n"
	   "}",
	   gl->expr->str,
	   gl->reversed, 
	   (gl->use_scalar != NULL));
}

static void gl_linear_set_simulation (GfsGl * object, GfsSimulation * sim)
{
  GfsGlLinear * gls = GFS_GL_LINEAR (object);
  GtsFile * fp = NULL;
  
  if (GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_linear_class ())->parent_class)->set_simulation)
    (*GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_linear_class ())->parent_class)->set_simulation)
      (object, sim);

  if (gls->nx) {
    gts_object_destroy (GTS_OBJECT (gls->nx));
    gts_object_destroy (GTS_OBJECT (gls->ny));
    gls->nx = gls->ny = NULL;
  }

  if (gls->expr->str[0] == '\0' || (fp = gfs_gl_linear_set (gls, gls->expr->str)))
    gfs_gl_linear_set (gls, "0");
  if (fp)
    gts_file_destroy (fp);
}

static void gl_linear_class_init (GfsGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gl_linear_destroy;
  GTS_OBJECT_CLASS (klass)->read  = gl_linear_read;
  GTS_OBJECT_CLASS (klass)->write = gl_linear_write;
  klass->set_simulation = gl_linear_set_simulation;
  klass->draw = gl_linear_draw;
}

static void gl_linear_init (GfsGlLinear * gl)
{
  GtsColor c = { 1., 1., 1. };

  gl->expr = g_string_new ("");
  gl->vf = gfs_gl_var_func_new ();
  GFS_GL (gl)->lc = c;
  gl->use_scalar = (GfsVariable *) TRUE;
}

GfsGlClass * gfs_gl_linear_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_linear_info = {
      "GfsGlLinear",
      sizeof (GfsGlLinear),
      sizeof (GfsGlScalarClass),
      (GtsObjectClassInitFunc) gl_linear_class_init,
      (GtsObjectInitFunc) gl_linear_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_scalar_class ()),
				  &gfs_gl_linear_info);
  }

  return klass;
}

static void normals (FttCell * cell, GfsGlLinear * gl)
{
  gdouble size = ftt_cell_size (cell);
  GFS_VALUE (cell, gl->nx) = gfs_center_gradient (cell, FTT_X, gl->vf->v->i)/size;
  GFS_VALUE (cell, gl->ny) = gfs_center_gradient (cell, FTT_Y, gl->vf->v->i)/size;
}

GtsFile * gfs_gl_linear_set (GfsGlLinear * gl, const gchar * func)
{
  GtsFile * fp;

  g_return_val_if_fail (gl != NULL, NULL);
  g_return_val_if_fail (func != NULL, NULL);

  if ((fp = gfs_gl_var_func_set (gl->vf, GFS_GL (gl)->sim, func, gl->expr, NULL)))
    return fp;

  if (gfs_function_get_constant_value (gl->vf->f) != 0.) {
    GfsDomain * domain = GFS_DOMAIN (GFS_GL (gl)->sim);
    if (!gl->nx) {
      gl->nx = gfs_temporary_variable (domain);
      gl->nx->component = FTT_X;
      gl->ny = gfs_temporary_variable (domain);
      gl->nx->component = FTT_Y;
    }
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) normals, gl);
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, gl->nx);
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, gl->ny);
  }

  return NULL;
}

/* GfsGlIsoline: Object */

static gboolean iso_cuts_cell (FttCell * cell, gpointer data)
{
  IsoParams * par = data;
  return (par->level >= GFS_VALUE (cell, par->min) &&
	  par->level <= GFS_VALUE (cell, par->max));
}

/**
 * gfs_gl_cell_traverse_visible_iso:
 * @gl: a #GfsGl2D.
 * @f: a view frustum.
 * @min: a variable storing the minimum level.
 * @max: a variable storing the maximum level.
 * @level: the level of the isoline/surface.
 * @func: a used-defined function.
 * @data: user data to pass to @func.
 *
 * Traverse the cells of @gl which are visible and whose bounding
 * sphere is intersected by the isoline/surface defined by @min, @max
 * and @level.
 */
void gfs_gl_cell_traverse_visible_iso (GfsGl * gl,
				       GfsFrustum * f,
				       GfsVariable * min,
				       GfsVariable * max,
				       gdouble level,
				       FttCellTraverseFunc func,
				       gpointer data)
{
  IsoParams p = { min, max, level };

  g_return_if_fail (gl != NULL);
  g_return_if_fail (f != NULL);
  g_return_if_fail (min != NULL);
  g_return_if_fail (max != NULL);
  g_return_if_fail (func != NULL);

  gfs_gl_cell_traverse_visible_condition (gl, f, iso_cuts_cell, &p, func, data);
}

static void gl_isoline_destroy (GtsObject * o)
{
  GfsGlIsoline * gl = GFS_GL_ISOLINE (o);

  if (gl->min)
    gts_object_destroy (GTS_OBJECT (gl->min));
  if (gl->max)
    gts_object_destroy (GTS_OBJECT (gl->max));
  g_array_free (gl->levels, TRUE);
  g_free (gl->ls);

  (* GTS_OBJECT_CLASS (gfs_gl_isoline_class ())->parent_class->destroy) (o);
}

static void gl_isoline_update_levels (GfsGl * gl)
{
  GfsGlScalar * gls = GFS_GL_SCALAR (gl);
  GfsGlIsoline * gli = GFS_GL_ISOLINE (gl);

  g_array_set_size (gli->levels, 0);
  if (gli->n > 0.) {
    guint i;
    for (i = 0; i < gli->n; i++) {
      gdouble l = gls->min + (i + 1)*(gls->max - gls->min)/(gli->n + 1.);
      g_array_append_val (gli->levels, l);
    }
  }

  if (gli->ls) {
    gchar ** list, ** s;

    s = list = g_strsplit (gli->ls, ",", -1);
    while (*s) {
      char * end;
      gdouble l = strtod (*s, &end);
      if (*end == '\0')
	g_array_append_val (gli->levels, l);
      s++;
    }
    g_strfreev (list);
  }
}

static void gl_isoline_read (GtsObject ** o, GtsFile * fp)
{
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "n",      TRUE},
    {GTS_STRING, "levels", TRUE},
    {GTS_NONE}
  };

  (* GTS_OBJECT_CLASS (gfs_gl_isoline_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  var[0].data = &GFS_GL_ISOLINE (*o)->n;
  var[1].data = &GFS_GL_ISOLINE (*o)->ls;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  gl_isoline_update_levels (GFS_GL (*o));
}

static void gl_isoline_write (GtsObject * o, FILE * fp)
{
  const gchar * levels = GFS_GL_ISOLINE (o)->ls;

  (* GTS_OBJECT_CLASS (gfs_gl_isoline_class ())->parent_class->write) (o, fp);

  fprintf (fp, " {\n  n = %g", GFS_GL_ISOLINE (o)->n);
  if (levels && levels[0] != '\0')
    fprintf (fp, " levels = %s", levels);
  fprintf (fp, "\n}");
}

static void min_max_iso (FttCell * cell, GfsGlIsoline * gl)
{
  gdouble min = G_MAXDOUBLE, max = -G_MAXDOUBLE;
  guint i;

  if (FTT_CELL_IS_LEAF (cell)) {
    for (i = 0; i < NCORNERS; i++) {
      gdouble v = gfs_cell_corner_value (cell, corner[i], GFS_GL_SCALAR (gl)->v, -1);
      if (v < min) min = v;
      if (v > max) max = v;
    }
  }
  else {
    FttCellChildren child;

    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i]) {
	gdouble vmin = GFS_VALUE (child.c[i], gl->min);
	gdouble vmax = GFS_VALUE (child.c[i], gl->max);
	if (vmin < min) min = vmin;
	if (vmax > max) max = vmax;
      }
  }
  GFS_VALUE (cell, gl->min) = min;
  GFS_VALUE (cell, gl->max) = max;
}

static GtsFile * gl_isoline_set_scalar (GfsGlScalar * gl, const gchar * func)
{
  GtsFile * fp = (*GFS_GL_SCALAR_CLASS 
		  (GTS_OBJECT_CLASS (gfs_gl_isoline_class ())->parent_class)->set_scalar)
    (gl, func);
  if (fp != NULL)
    return fp;
  
  gfs_domain_cell_traverse (GFS_DOMAIN (GFS_GL (gl)->sim),
			    FTT_POST_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) min_max_iso, gl);  
  return fp;
}

static void gl_isoline_set_simulation (GfsGl * object, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsGlIsoline * gls = GFS_GL_ISOLINE (object);

  if (gls->min)
    gts_object_destroy (GTS_OBJECT (gls->min));
  gls->min = gfs_temporary_variable (domain);
  if (gls->max)
    gts_object_destroy (GTS_OBJECT (gls->max));
  gls->max = gfs_temporary_variable (domain);

  (*GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_isoline_class ())->parent_class)->set_simulation)
    (object, sim);
}

static void gl_isoline_class_init (GfsGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gl_isoline_destroy;
  GTS_OBJECT_CLASS (klass)->read = gl_isoline_read;
  GTS_OBJECT_CLASS (klass)->write = gl_isoline_write;
  klass->set_simulation = gl_isoline_set_simulation;
  GFS_GL_SCALAR_CLASS (klass)->set_scalar = gl_isoline_set_scalar;
  klass->draw = gl_isoline_draw;
  klass->pick = NULL;
}

static void gl_isoline_init (GfsGlIsoline * gl)
{
  GtsColor c = { 0., 0., 0. };
  GFS_GL (gl)->lc = c;
  gl->levels = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gl->n = 10.;
}

GfsGlClass * gfs_gl_isoline_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_isoline_info = {
      "GfsGlIsoline",
      sizeof (GfsGlIsoline),
      sizeof (GfsGlScalarClass),
      (GtsObjectClassInitFunc) gl_isoline_class_init,
      (GtsObjectInitFunc) gl_isoline_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_linear_class ()),
				  &gfs_gl_isoline_info);
  }

  return klass;
}

void gfs_gl_isoline_set_levels (GfsGlIsoline * gl, const gchar * levels)
{
  g_return_if_fail (gl != NULL);
  g_return_if_fail (levels != NULL);

  g_free (gl->ls);
  gl->ls = g_strdup (levels);
}

/* GfsGlVOF: Object */

static void gl_vof_destroy (GtsObject * o)
{
  GfsGlVOF * gl = GFS_GL_VOF (o);

  gfs_gl_var_func_destroy (gl->vf);
  g_string_free (gl->expr, TRUE);

  (* GTS_OBJECT_CLASS (gfs_gl_vof_class ())->parent_class->destroy) (o);
}

static void gl_vof_read (GtsObject ** o, GtsFile * fp)
{
  GfsGlVOF * gl = GFS_GL_VOF (*o);
  GtsFileVariable var[] = {
    {GTS_INT, "reversed",    TRUE},
    {GTS_INT, "use_scalar",  TRUE},
    {GTS_INT, "draw_edges",  TRUE},
    {GTS_INT, "interpolate", TRUE},
    {GTS_NONE}
  };

  (* GTS_OBJECT_CLASS (gfs_gl_vof_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  g_string_free (gl->expr, TRUE);
  if (!(gl->expr = gfs_function_expression (fp, NULL)))
    return;
  gts_file_next_token (fp);

  var[0].data = &gl->reversed;
  var[1].data = &gl->use_scalar;
  var[2].data = &gl->draw_edges;
  var[3].data = &gl->interpolate;
  gts_file_assign_variables (fp, var);
}

static void gl_vof_write (GtsObject * o, FILE * fp)
{
  GfsGlVOF * gl = GFS_GL_VOF (o);

  (* GTS_OBJECT_CLASS (gfs_gl_vof_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s {\n"
	   "  reversed = %d\n"
	   "  use_scalar = %d\n"
	   "  draw_edges = %d\n"
	   "  interpolate = %d\n"
	   "}",
	   gl->expr->str,
	   gl->reversed, 
	   (gl->use_scalar != NULL),
	   gl->draw_edges,
	   gl->interpolate);
}

GtsFile * gfs_gl_vof_set (GfsGlVOF * gl, const gchar * func)
{
  GtsFile * fp;

  g_return_val_if_fail (gl != NULL, NULL);
  g_return_val_if_fail (func != NULL, NULL);

  if ((fp = gfs_gl_var_func_set (gl->vf, GFS_GL (gl)->sim, func, gl->expr, 
				 GFS_VARIABLE_CLASS (gfs_variable_tracer_vof_class ()))))
    return fp;
  return NULL;
}

static void gl_vof_set_simulation (GfsGl * object, GfsSimulation * sim)
{
  GfsGlVOF * gls = GFS_GL_VOF (object);
  GtsFile * fp = NULL;

  (*GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_vof_class ())->parent_class)->set_simulation)
    (object, sim);

  if (gls->expr->str[0] == '\0' || (fp = gfs_gl_vof_set (gls, gls->expr->str))) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    
    if (domain->variables) {
      GfsVariable * v = NULL;
      GSList * i = domain->variables;
      
      while (i && !v) {
	if (GFS_IS_VARIABLE_TRACER_VOF (i->data))
	  v = i->data;
	i = i->next;
      }
      if (!v)
	v = domain->variables->data;
      gfs_gl_vof_set (gls, v->name);
    }
    else
      gfs_gl_vof_set (gls, "0");
  }
  if (fp)
    gts_file_destroy (fp);
}

static void gl_vof_class_init (GfsGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gl_vof_destroy;
  GTS_OBJECT_CLASS (klass)->read = gl_vof_read;
  GTS_OBJECT_CLASS (klass)->write = gl_vof_write;
  klass->set_simulation = gl_vof_set_simulation;
  klass->draw = gl_vof_draw;
  klass->cut = gl_vof_cut;
  klass->pick = NULL;
}

static void gl_vof_init (GfsGlVOF * gl)
{
  gl->expr = g_string_new ("");
  gl->vf = gfs_gl_var_func_new ();
#if !FTT_2D
  GFS_GL (gl)->lc.r = 1.; GFS_GL (gl)->lc.g = 1.; GFS_GL (gl)->lc.b = 1.;
#endif
}

GfsGlClass * gfs_gl_vof_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_vof_info = {
      "GfsGlVOF",
      sizeof (GfsGlVOF),
      sizeof (GfsGlScalarClass),
      (GtsObjectClassInitFunc) gl_vof_class_init,
      (GtsObjectInitFunc) gl_vof_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_scalar_class ()),
				  &gfs_gl_vof_info);
  }

  return klass;
}

/* GfsGlVectors: Object */

static void gl_vectors_read (GtsObject ** o, GtsFile * fp)
{
  GfsGlVectors * gl = GFS_GL_VECTORS (*o);
  FttComponent c;
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "scale",      TRUE},
    {GTS_INT,    "use_scalar", TRUE},
    {GTS_NONE}
  };

  (* GTS_OBJECT_CLASS (gfs_gl_vectors_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  for (c = 0; c < FTT_DIMENSION; c++) {
    g_string_free (gl->expr[c], TRUE);
    if (!(gl->expr[c] = gfs_function_expression (fp, NULL)))
      return;
    gts_file_next_token (fp);
  }

  var[0].data = &gl->scale;
  var[1].data = &gl->use_scalar;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;
  if (var[0].set)
    gl->already_set = TRUE;
}

static void gl_vectors_write (GtsObject * o, FILE * fp)
{
  GfsGlVectors * gl = GFS_GL_VECTORS (o);
  FttComponent c;

  (* GTS_OBJECT_CLASS (gfs_gl_vectors_class ())->parent_class->write) (o, fp);

  for (c = 0; c < FTT_DIMENSION; c++)
    fprintf (fp, " %s", gl->expr[c]->str);
  fprintf (fp, 
	   " {\n"
	   "  scale = %g\n"
	   "  use_scalar = %d\n"
	   "}",
	   gl->scale, gl->use_scalar);
}

static void gl_vectors_destroy (GtsObject * object)
{
  GfsGlVectors * gl = GFS_GL_VECTORS (object);
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++) {
    gfs_gl_var_func_destroy (gl->vf[c]);
    g_string_free (gl->expr[c], TRUE);
  }

  (* GTS_OBJECT_CLASS (gfs_gl_vectors_class ())->parent_class->destroy) (object);
}

static void maxv (FttCell * cell, GfsGlVectors * gl)
{
  FttComponent c;
  gdouble n2 = 0.;
  gdouble size = ftt_cell_size (cell);

  for (c = 0; c < FTT_DIMENSION; c++)
    n2 += GFS_VALUE (cell, gl->v[c])*GFS_VALUE (cell, gl->v[c]);
  if (n2 > gl->max)
    gl->max = n2;
  if (size < gl->h)
    gl->h = size;
}

static void update_norm (GfsGlVectors * gl)
{
  gl->max = -G_MAXDOUBLE;
  gl->h = G_MAXDOUBLE;
  gfs_domain_cell_traverse (GFS_DOMAIN (GFS_GL (gl)->sim),
			    FTT_POST_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) maxv, gl);
  gl->max = gl->max >= 0. ? sqrt (gl->max) : 1.;
  if (!gl->already_set) {
    gl->scale = gl->max > 0. ? gl->h/gl->max : 1.;
    gl->already_set = TRUE;
  }
}

static void gl_vectors_set_simulation (GfsGl * object, GfsSimulation * sim)
{
  GfsGlVectors * gls = GFS_GL_VECTORS (object);
  GfsDomain * domain;
  FttComponent c;

  (*GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_vectors_class ())->parent_class)->set_simulation)
    (object, sim);

  domain = GFS_DOMAIN (sim);
  for (c = 0; c < FTT_DIMENSION; c++) {
    GtsFile * fp = NULL;

    if (gls->expr[c]->str[0] == '\0' ||
	(fp = gfs_gl_var_func_set (gls->vf[c], sim, gls->expr[c]->str, NULL, NULL))) {
      GfsVariable ** u = gfs_domain_velocity (domain);

      if (u)
	gfs_gl_var_func_set (gls->vf[c], sim, u[c]->name, gls->expr[c], NULL);
      else if (domain->variables) {
	GfsVariable * v = domain->variables->data;
	gfs_gl_var_func_set (gls->vf[c], sim, v->name, gls->expr[c], NULL);
      }
      else
	gfs_gl_var_func_set (gls->vf[c], sim, "0", gls->expr[c], NULL);
    }
    if (fp)
      gts_file_destroy (fp);
    gls->v[c] = gls->vf[c]->v;
  }
  update_norm (gls);
}

static void gl_vectors_draw (GfsGl * gl, GfsFrustum * f)
{
  GfsGlVectors * gls = GFS_GL_VECTORS (gl);
  
  gl->size = 0;
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glTranslatef (0., 0., gl->p->lc);
  glBegin (GL_LINES);
  gfs_gl_normal (gl);
  gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc) gl_vector, gl);
  glEnd ();
  glPopMatrix ();

  if (gls->use_scalar)
    (* GFS_GL_CLASS (GTS_OBJECT (gl)->klass->parent_class)->draw) (gl, f);  
}

static void gl_vectors_class_init (GfsGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_vectors_read;
  GTS_OBJECT_CLASS (klass)->write = gl_vectors_write;
  GTS_OBJECT_CLASS (klass)->destroy = gl_vectors_destroy;
  klass->set_simulation = gl_vectors_set_simulation;
  klass->draw = gl_vectors_draw;
  klass->pick = NULL;
}

static void gl_vectors_init (GfsGlVectors * object)
{
  FttComponent c;
  
  for (c = 0; c < FTT_DIMENSION; c++) {
    object->vf[c] = gfs_gl_var_func_new ();
    object->expr[c] = g_string_new ("");
  }
  object->scale = 1.;
}

GfsGlClass * gfs_gl_vectors_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_vectors_info = {
      "GfsGlVectors",
      sizeof (GfsGlVectors),
      sizeof (GfsGlScalarClass),
      (GtsObjectClassInitFunc) gl_vectors_class_init,
      (GtsObjectInitFunc) gl_vectors_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_scalar_class ()),
				  &gfs_gl_vectors_info);
  }

  return klass;
}

GtsFile * gfs_gl_vectors_set (GfsGlVectors * gl, FttComponent c, const gchar * func)
{
  GtsFile * fp;

  g_return_val_if_fail (gl != NULL, NULL);
  g_return_val_if_fail (c < FTT_DIMENSION, NULL);
  g_return_val_if_fail (func != NULL, NULL);

  if ((fp = gfs_gl_var_func_set (gl->vf[c], GFS_GL (gl)->sim, func, gl->expr[c], NULL)))
    return fp;

  gl->v[c] = gl->vf[c]->v;
  update_norm (gl);

  return NULL;
}

/* GfsGlStreamline: Object */

static void gl_streamline_destroy (GtsObject * o)
{
  gfs_streamline_destroy (GFS_GL_STREAMLINE (o)->l);
  glDeleteLists (GFS_GL_STREAMLINE (o)->list, 1);

  (* GTS_OBJECT_CLASS (gfs_gl_streamline_class ())->parent_class->destroy) (o);
}

static void gl_streamline_read (GtsObject ** o, GtsFile * fp)
{
  GfsGlStreamline * gl = GFS_GL_STREAMLINE (*o);
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x", TRUE},
    {GTS_DOUBLE, "y", TRUE},
    {GTS_DOUBLE, "z", TRUE},
    {GTS_NONE}
  };

  var[0].data = &gl->c.x;
  var[1].data = &gl->c.y;
  var[2].data = &gl->c.z;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;
}

static void gl_streamline_write (GtsObject * o, FILE * fp)
{
  GfsGlStreamline * gl = GFS_GL_STREAMLINE (o);

  fprintf (fp,
	   " {\n"
	   "    x = %.16f y = %.16f z = %.16f\n"
	   "  }",
	   gl->c.x, gl->c.y, gl->c.z);
}

static void gl_streamline_class_init (GtsObjectClass * klass)
{
  klass->destroy = gl_streamline_destroy;
  klass->read = gl_streamline_read;
  klass->write = gl_streamline_write;
}

GtsObjectClass * gfs_gl_streamline_class (void)
{
  static GtsObjectClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_streamline_info = {
      "GfsGlStreamline",
      sizeof (GfsGlStreamline),
      sizeof (GtsObjectClass),
      (GtsObjectClassInitFunc) gl_streamline_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (gts_object_class (), &gfs_gl_streamline_info);
  }

  return klass;
}

static void gfs_gl_streamline_write (GfsGlStreamline * s, FILE * fp)
{
  (* GTS_OBJECT (s)->klass->write) (GTS_OBJECT (s), fp);
}

/* GfsGlStreamlines: Object */

static void gl_streamlines_destroy (GtsObject * o)
{
  GfsGlStreamlines * gl = GFS_GL_STREAMLINES (o);

  gfs_gl_streamlines_reset (gl);
  g_list_foreach (gl->stream, (GFunc) gts_object_destroy, NULL);
  g_list_free (gl->stream);

  if (gl->s)
    gts_object_destroy (GTS_OBJECT (gl->s));

  (* GTS_OBJECT_CLASS (gfs_gl_streamlines_class ())->parent_class->destroy) (o);
}

static void gl_streamlines_read (GtsObject ** o, GtsFile * fp)
{
  GfsGlStreamlines * gl = GFS_GL_STREAMLINES (*o);
  GtsFileVariable var[] = {
    {GTS_INT,    "show_cells", TRUE},
    {GTS_DOUBLE, "dmin",       TRUE},
    {GTS_DOUBLE, "radius",     TRUE},
    {GTS_NONE}
  };

  (* GTS_OBJECT_CLASS (gfs_gl_streamlines_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  
  var[0].data = &gl->show_cells;
  var[1].data = &gl->dmin;
  var[2].data = &gl->radius;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return;
  }
  fp->scope_max++;
  gts_file_next_token (fp);
  while (fp->type == '\n')
    gts_file_next_token (fp);
  while (fp->type == '{') {
    GtsObject * o = gts_object_new (gfs_gl_streamline_class ());

    (* o->klass->read) (&o, fp);
    if (fp->type == GTS_ERROR) {
      gts_object_destroy (o);
      return;
    }
    gl->stream = g_list_append (gl->stream, o);
  }
  while (fp->type == '\n')
    gts_file_next_token (fp);
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);
}

static void gl_streamlines_write (GtsObject * o, FILE * fp)
{
  GfsGlStreamlines * gl = GFS_GL_STREAMLINES (o);

  (* GTS_OBJECT_CLASS (gfs_gl_streamlines_class ())->parent_class->write) (o, fp);

  fprintf (fp, " { show_cells = %d dmin = %g radius = %g } {\n ", 
	   gl->show_cells, gl->dmin, gl->radius);
  g_list_foreach (gl->stream, (GFunc) gfs_gl_streamline_write, fp);
  fputs ("\n}", fp);
}

static void gl_streamlines_set_simulation (GfsGl * object, GfsSimulation * sim)
{
  GfsGlStreamlines * gls = GFS_GL_STREAMLINES (object);
  
  gfs_gl_streamlines_reset (gls);

  (*GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_streamlines_class ())->parent_class)->set_simulation)
    (object, sim);

  if (gls->s)
    gts_object_destroy (GTS_OBJECT (gls->s));
  gls->s = gfs_temporary_variable (GFS_DOMAIN (sim));
  gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, gls->s);
}

static gboolean cell_overlaps_segment (FttCell * cell, gpointer data)
{
  GList * l = data;
  GtsBBox bb;
  GtsSegment s;

  ftt_cell_bbox (cell, &bb);
  s.v1 = l->data; s.v2 = l->next->data;
  return gts_bbox_overlaps_segment (&bb, &s);
}

static void add_segment (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GList * l = data[1];

  if (FTT_CELL_IS_LEAF (cell))
    GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, v)) = 
      g_slist_prepend (GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, v)), l);
  else
    GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, v)) = l;
}

static gdouble distance2 (GList * i, GtsPoint * p, GtsSegment * closest)
{
  gdouble d1 = G_MAXDOUBLE, d2 = G_MAXDOUBLE;
  GtsSegment s;

  s.v1 = i->data;
  if (i->next) {
    s.v2 = i->next->data;
    d1 = gts_point_segment_distance2 (p, &s);
  }
  if (i->prev) {
    s.v2 = i->prev->data;
    d2 = gts_point_segment_distance2 (p, &s);
  }
  if (closest) {
    closest->v1 = s.v1;
    closest->v2 = d1 < d2 ? i->next->data : i->prev->data;
  }
  return MIN (d1, d2);
}

static gboolean distance_is_larger (GList * a, GList * b, guint d)
{
  while (a && a != b && --d)
    a = a->next;
  return !d;
}

static gdouble cell_distance2 (FttCell * cell, GtsPoint * p, gpointer data)
{
  GfsGlStreamlines * gls = ((gpointer *) data)[0];
  GSList * i = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, gls->s));

  if (i == NULL)
    return G_MAXDOUBLE;
  else if (!FTT_CELL_IS_LEAF (cell))
    return ftt_cell_point_distance2_min (cell, p);
  else {
    GfsGlStreamline * s = ((gpointer *) data)[1];
    GList * l = ((gpointer *) data)[2];
    gdouble dmin = G_MAXDOUBLE, h = ftt_cell_size (cell);

    while (i) {
      GList * j = i->data;
      gboolean same_stream = (GTS_OBJECT (j->data)->reserved == s);
      
      if (!same_stream || distance_is_larger (l, j, 16)) {
	gdouble d = distance2 (j, p, NULL);

	if (d < dmin)
	  dmin = d;
	if (same_stream && d < h*h/100.)
	  return 0.;
      }
      i = i->next;
    }
    return dmin;
  }
}

static gboolean stop (FttCell * cell, GList * l, gpointer data)
{
  GfsGlStreamlines * gls = ((gpointer *) data)[0];
  GfsGlStreamline * s = ((gpointer *) data)[1];
  gboolean stop = FALSE;

  ((gpointer *) data)[2] = l;
  stop = (gfs_domain_cell_point_distance2 (GFS_DOMAIN (GFS_GL (gls)->sim), l->data,
					   cell_distance2, data, NULL) 
	  <= gls->dmin*gls->dmin/4.);
  if (l->next) {
    gpointer datum[2];
    
    datum[0] = gls->s;
    datum[1] = l;
    gfs_domain_cell_traverse_condition (GFS_DOMAIN (GFS_GL (gls)->sim), 
					FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
					(FttCellTraverseFunc) add_segment, datum,
					cell_overlaps_segment, l);
  }
  GTS_OBJECT (l->data)->reserved = s;
  return stop;
}

static void gl_streamline_update_display_list (GfsGlStreamline * s,
					       GfsGlStreamlines * gls)
{
  if (s->list == 0) {
    s->list = glGenLists (1);
    if (s->list == 0)
      g_warning ("No available OpenGL display list!");
  }
  glNewList (s->list, GL_COMPILE);
  gl_streamline_draw (s, gls);
  glEndList ();
}

static void gl_streamline_update (GfsGlStreamline * s,
				  GfsGlStreamlines * gls)
{
  GtsPoint p;
  gpointer data[3];
  
  if (s->l)
    return;

  data[0] = gls;
  data[1] = s;
  p.x = s->c.x; p.y = s->c.y; p.z = s->c.z;
  if (gfs_domain_cell_point_distance2 (GFS_DOMAIN (GFS_GL (gls)->sim), &p,
				       cell_distance2, data, NULL) > gls->dmin*gls->dmin) {
    s->l = gfs_streamline_new (GFS_DOMAIN (GFS_GL (gls)->sim),
			       GFS_GL_VECTORS (gls)->v, s->c,
			       GFS_GL_VECTORS (gls)->use_scalar ? GFS_GL_SCALAR (gls)->v : NULL,
			       GFS_GL_SCALAR (gls)->min, GFS_GL_SCALAR (gls)->max, FALSE,
			       stop, data);
    if (s->l && !s->l->next) {
      gfs_streamline_destroy (s->l);
      s->l = NULL;
    }
    else
      gl_streamline_update_display_list (s, gls);
  }
}

static void gl_streamlines_draw (GfsGl * gl, GfsFrustum * f)
{
  GfsGlStreamlines * gls = GFS_GL_STREAMLINES (gl);
  GList * i = gls->stream;

  gl->size = 0;
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glTranslatef (0., 0., gl->p->lc/2.);
  gfs_gl_normal (gl);
  if (gls->show_cells) {
    glColor3f (0.3, 0.3, 0.3);
    gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc) gl_cell, gl);
  }
  glTranslatef (0., 0., gl->p->lc/2.);
  while (i) {
    GfsGlStreamline * s = i->data;

    if (i == gls->selected)
      glColor3f (1., 0., 0.);
    else if (!GFS_GL_VECTORS (gl)->use_scalar)
      glColor3f (gl->lc.r, gl->lc.g, gl->lc.b);
    gl_streamline_update (s, gls);
    if (s->l) {
      if (gl->format == GFSGL_PPM_OFFSCREEN)
	gl_streamline_draw (s, gls);
      else
	glCallList (s->list);
    }
    gl->size++;
    i = i->next;
  }
  glPopMatrix ();
}

static gdouble gl_streamlines_pick (GfsGl * gl, GfsGlRay * r)
{
  gint maxlevel = gl->maxlevel;
  gdouble a;

  gl->maxlevel = -1;
  a = gl2D_pick (gl, r);
  gl->maxlevel = maxlevel;
  if (a != GFS_NODATA) {
    GfsGl2D * gl2 = GFS_GL2D (gl);
    GfsGlStreamlines * gls = GFS_GL_STREAMLINES (gl);
    GSList * i = GFS_DOUBLE_TO_POINTER (GFS_VALUE (gl2->picked, gls->s));
    gdouble dmin = G_MAXDOUBLE;
    GfsGlStreamline * closest = NULL;
    GtsPoint p;

    p.x = gl2->pickedpos.x; p.y = gl2->pickedpos.y; p.z = gl2->pickedpos.z;
    while (i) {
      GList * j = i->data;
      gdouble d = distance2 (j, &p, NULL);

      if (d < dmin) {
	dmin = d;
	closest = GTS_OBJECT (j->data)->reserved;
      }
      i = i->next;
    }
    gls->candidate = closest ? g_list_find (gls->stream, closest) : NULL;
  }
  return a;
}

static void gl_streamlines_class_init (GfsGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gl_streamlines_destroy;
  GTS_OBJECT_CLASS (klass)->read = gl_streamlines_read;
  GTS_OBJECT_CLASS (klass)->write = gl_streamlines_write;
  klass->set_simulation = gl_streamlines_set_simulation;
  klass->draw = gl_streamlines_draw;
  klass->pick = gl_streamlines_pick;
}

static void gl_streamlines_init (GfsGlStreamlines * gl)
{
  GtsColor c = { 1., 1., 1. };

  GFS_GL (gl)->lc = c;
  gl->show_cells = TRUE;
}

GfsGlClass * gfs_gl_streamlines_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_streamlines_info = {
      "GfsGlStreamlines",
      sizeof (GfsGlStreamlines),
      sizeof (GfsGlScalarClass),
      (GtsObjectClassInitFunc) gl_streamlines_class_init,
      (GtsObjectInitFunc) gl_streamlines_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_vectors_class ()),
				  &gfs_gl_streamlines_info);
  }

  return klass;
}

static void reset_segments (FttCell * cell, GfsGlStreamlines * gl)
{
  if (FTT_CELL_IS_LEAF (cell)) {
    GSList * l = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, gl->s));
    g_slist_free (l);
  }
  GFS_VALUE (cell, gl->s) = 0.;
}

void gfs_gl_streamlines_reset (GfsGlStreamlines * gl)
{
  GList * i;

  g_return_if_fail (gl != NULL);

  if (GFS_GL (gl)->sim)
    gfs_domain_cell_traverse (GFS_DOMAIN (GFS_GL (gl)->sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) reset_segments, gl);
  i = gl->stream;
  while (i) {
    GfsGlStreamline * s = i->data;
    gfs_streamline_destroy (s->l);
    s->l = NULL;    
    i = i->next;
  }
}

static void remove_segment (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];

  if (FTT_CELL_IS_LEAF (cell)) {
    GList * i = data[1];
    GSList * l = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, v));
    l = g_slist_remove (l, i);
    l = g_slist_remove (l, i->next);
    GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, v)) = l;
  }
  else {
    FttCellChildren child;
    gboolean empty = TRUE;
    guint i;
    
    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS && empty; i++)
      if (child.c[i] && GFS_VALUE (child.c[i], v) != 0.)
	empty = FALSE;
    if (empty)
      GFS_VALUE (cell, v) = 0.;
  }
}

void gfs_gl_streamlines_reset_selected (GfsGlStreamlines * gl)
{
  GfsGlStreamline * s;
  gpointer datum[2];
  GList * i;

  g_return_if_fail (gl != NULL);
  g_return_if_fail (gl->selected != NULL);

  s = gl->selected->data;
  datum[0] = gl->s;
  i = s->l;
  while (i && i->next) {
    datum[1] = i;
    gfs_domain_cell_traverse_condition (GFS_DOMAIN (GFS_GL (gl)->sim), 
					FTT_POST_ORDER, FTT_TRAVERSE_ALL, -1,
					(FttCellTraverseFunc) remove_segment, datum,
					cell_overlaps_segment, i);
    i = i->next;
  }
  gfs_streamline_destroy (s->l);
  s->l = NULL;
}

GfsGlStreamline * gfs_gl_streamlines_add (GfsGlStreamlines * gl, FttVector p)
{
  GfsGlStreamline * s;

  g_return_val_if_fail (gl != NULL, NULL);

  s = GFS_GL_STREAMLINE (gts_object_new (gfs_gl_streamline_class ()));
  s->c = p;
  gl_streamline_update (s, gl);
  if (s->l) {
    gl->stream = g_list_append (gl->stream, s);
    gl->selected = g_list_last (gl->stream);
    return s;
  }
  else {
    gts_object_destroy (GTS_OBJECT (s));
    return NULL;
  }
}

gboolean gfs_gl_streamlines_remove_selected (GfsGlStreamlines * gl)
{
  g_return_val_if_fail (gl != NULL, FALSE);

  if (gl->selected) {
    gfs_gl_streamlines_reset_selected (gl);
    gts_object_destroy (gl->selected->data);
    gl->stream = g_list_remove_link (gl->stream, gl->selected);
    g_list_free (gl->selected);
    gl->selected = NULL;
    return TRUE;
  }
  return FALSE;
}

void gfs_gl_streamlines_update_display_lists (GfsGlStreamlines * gl)
{
  g_return_if_fail (gl != NULL);

  g_list_foreach (gl->stream, (GFunc) gl_streamline_update_display_list, gl);
}

gdouble gfs_gl_streamlines_closest (GfsGlStreamlines * gl, 
				    FttVector * p,
				    GtsPoint * closest)
{
  FttCell * cell = NULL;
  GtsPoint p1;
  gdouble dmin;
  gpointer data[3];

  g_return_val_if_fail (gl != NULL, G_MAXDOUBLE);
  g_return_val_if_fail (p != NULL, G_MAXDOUBLE);
  g_return_val_if_fail (closest != NULL, G_MAXDOUBLE);

  data[0] = gl; data[1] = data[2] = NULL;
  p1.x = p->x; p1.y = p->y; p1.z = p->z;
  dmin = gfs_domain_cell_point_distance2 (GFS_DOMAIN (GFS_GL (gl)->sim), &p1,
					  cell_distance2, data, &cell);
  if (cell) {
    GSList * i = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, gl->s));
    GtsSegment segment;
    gdouble dmin = G_MAXDOUBLE;
    
    while (i) {
      GList * j = i->data;
      GtsSegment s1;
      gdouble d = distance2 (j, &p1, &s1);
      
      if (d < dmin) {
	dmin = d;
	segment = s1;
      }
      i = i->next;
    }
    gts_point_segment_closest (&p1, &segment, closest);
  }
  return dmin;
}

static GList * streamline_step (GList * i, GtsPoint * p, gdouble ds)
{
  if (!i || !i->next)
    return NULL;
  if (p->x == G_MAXDOUBLE) {
    GtsPoint * p1 = i->data;
    p->x = p1->x; p->y = p1->y; p->z = p1->z;
    return i;
  }
  else {
    gdouble l = gts_point_distance (i->data, i->next->data);
    gdouble s;
  
    g_assert (l > 0.);
    s = (gts_point_distance (i->data, p) + ds)/l;
    while (s > 1. && i->next->next) {
      s = l*(s - 1.);
      i = i->next;
      l = gts_point_distance (i->data, i->next->data);
      g_assert (l > 0.);
      s /= l;
    }
    if (s > 1.)
      return NULL;
    else {
      GtsPoint * p1 = i->data;
      GtsPoint * p2 = i->next->data;
      p->x = p1->x + s*(p2->x - p1->x);
      p->y = p1->y + s*(p2->y - p1->y);
      p->z = p1->z + s*(p2->z - p1->z);
      return i;
    }
  }
}

static void push_stream (GfsGlStreamline * s, GtsFifo * fifo)
{
  gts_fifo_push (fifo, s);
}

void gfs_gl_streamlines_evenly_spaced (GfsGlStreamlines * gl,
				       gboolean (* callback) (GfsGlStreamlines *, gpointer),
				       gpointer data)
{
  GtsFifo * fifo;
  GfsGlStreamline * s;
  gboolean stop = FALSE;

  g_return_if_fail (gl != NULL);
  g_return_if_fail (gl->dmin > 0.);

  fifo = gts_fifo_new ();
  g_list_foreach (gl->stream, (GFunc) push_stream, fifo);
  while ((s = gts_fifo_pop (fifo)) && !stop) {
    GtsPoint p;
    GList * i = s->l;

    g_assert (i);
    p.x = G_MAXDOUBLE;
    while ((i = streamline_step (i, &p, gl->dmin/10.))) {
      GtsVector n;
      FttVector p1;
      GfsGlStreamline * s1;

      gts_vector_init (n, GTS_POINT (i->data), GTS_POINT (i->next->data));
      gts_vector_normalize (n);

      p1.x = p.x + n[1]*gl->dmin;
      p1.y = p.y - n[0]*gl->dmin;
      p1.z = 0.;
      if ((s1 = gfs_gl_streamlines_add (gl, p1))) {
	if (callback)
	  stop |= (* callback) (gl, data);
	gts_fifo_push (fifo, s1);
      }

      p1.x = p.x - n[1]*gl->dmin;
      p1.y = p.y + n[0]*gl->dmin;
      p1.z = 0.;
      if ((s1 = gfs_gl_streamlines_add (gl, p1))) {
	if (callback)
	  stop |= (* callback) (gl, data);	
	gts_fifo_push (fifo, s1);
      }
    }
  }
  gts_fifo_destroy (fifo);
  gl->selected = NULL;
}

/* GfsGlEllipses: Object */

static void gl_ellipses_read (GtsObject ** o, GtsFile * fp)
{
  GfsGlEllipses * gl = GFS_GL_ELLIPSES (*o);
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "scale",      TRUE},
    {GTS_INT,    "use_scalar", TRUE},
    {GTS_NONE}
  };
  guint i;

  (* GTS_OBJECT_CLASS (gfs_gl_ellipses_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  for (i = 0; i < 4; i++) {
    g_string_free (gl->expr[i], TRUE);
    if (!(gl->expr[i] = gfs_function_expression (fp, NULL)))
      return;
    gts_file_next_token (fp);
  }

  var[0].data = &gl->scale;
  var[1].data = &gl->use_scalar;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;
  if (var[0].set)
    gl->already_set = TRUE;
}

static void gl_ellipses_write (GtsObject * o, FILE * fp)
{
  GfsGlEllipses * gl = GFS_GL_ELLIPSES (o);
  guint i;

  (* GTS_OBJECT_CLASS (gfs_gl_ellipses_class ())->parent_class->write) (o, fp);

  for (i = 0; i < 4; i++)
    fprintf (fp, " %s", gl->expr[i]->str);
  fprintf (fp, " {\n"
	   "  scale = %g\n"
	   "  use_scalar = %d\n"
	   "}",
	   gl->scale, gl->use_scalar);
}

static void gl_ellipses_destroy (GtsObject * o)
{
  GfsGlEllipses * gl = GFS_GL_ELLIPSES (o);
  guint i;

  for (i = 0; i < 4; i++) {
    gfs_gl_var_func_destroy (gl->vf[i]);
    g_string_free (gl->expr[i], TRUE);
  }

  (* GTS_OBJECT_CLASS (gfs_gl_ellipses_class ())->parent_class->destroy) (o);
}

static void maxe (FttCell * cell, GfsGlEllipses * gl)
{
  guint i;
  gdouble n2 = 0.;
  gdouble size = ftt_cell_size (cell);

  for (i = 0; i < 4; i++)
    n2 += GFS_VALUE (cell, gl->v[i])*GFS_VALUE (cell, gl->v[i]);
  if (n2 > gl->max)
    gl->max = n2;
  if (size < gl->h)
    gl->h = size;
}

static void update_ellipses_norm (GfsGlEllipses * gl)
{
  gl->max = -G_MAXDOUBLE;
  gl->h = G_MAXDOUBLE;
  gfs_domain_cell_traverse (GFS_DOMAIN (GFS_GL (gl)->sim),
			    FTT_POST_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) maxe, gl);
  gl->max = gl->max >= 0. ? sqrt (gl->max) : 1.;
  if (!gl->already_set) {
    gl->scale = gl->max > 0. ? gl->h/gl->max : 1.;
    gl->already_set = TRUE;
  }
}

static void gl_ellipses_set_simulation (GfsGl * object, GfsSimulation * sim)
{
  GfsGlEllipses * gle = GFS_GL_ELLIPSES (object);
  GfsDomain * domain;
  FttComponent c;

  (*GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_ellipses_class ())->parent_class)->set_simulation)
    (object, sim);

  domain = GFS_DOMAIN (sim);
  for (c = 0; c < 4; c++) {
    GtsFile * fp = NULL;

    if (gle->expr[c]->str[0] == '\0' ||
	(fp = gfs_gl_var_func_set (gle->vf[c], sim, gle->expr[c]->str, NULL, NULL))) {
      if (domain->variables) {
	GfsVariable * v = domain->variables->data;
	gfs_gl_var_func_set (gle->vf[c], sim, v->name, gle->expr[c], NULL);
      }
      else
	gfs_gl_var_func_set (gle->vf[c], sim, "0", gle->expr[c], NULL);
    }
    if (fp)
      gts_file_destroy (fp);
    gle->v[c] = gle->vf[c]->v;
  }
  update_ellipses_norm (gle);
}

static void gl_ellipses_draw (GfsGl * gl, GfsFrustum * f)
{
  GfsGlEllipses * gls = GFS_GL_ELLIPSES (gl);
  
  gl->size = 0;
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glTranslatef (0., 0., gl->p->lc);
  gfs_gl_normal (gl);
  gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc) gl_ellipse, gl);
  glPopMatrix ();

  if (gls->use_scalar)
    (* GFS_GL_CLASS (GTS_OBJECT (gl)->klass->parent_class)->draw) (gl, f);  
}

static void gl_ellipses_class_init (GfsGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_ellipses_read;
  GTS_OBJECT_CLASS (klass)->write = gl_ellipses_write;
  GTS_OBJECT_CLASS (klass)->destroy = gl_ellipses_destroy;
  klass->set_simulation = gl_ellipses_set_simulation;
  klass->draw = gl_ellipses_draw;
  klass->pick = NULL;
}

static void gl_ellipses_init (GfsGlEllipses * object)
{
  FttComponent c;
  
  for (c = 0; c < 4; c++) {
    object->vf[c] = gfs_gl_var_func_new ();
    object->expr[c] = g_string_new ("");
  }
  object->scale = 1.;
}

GfsGlClass * gfs_gl_ellipses_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_ellipses_info = {
      "GfsGlEllipses",
      sizeof (GfsGlEllipses),
      sizeof (GfsGlScalarClass),
      (GtsObjectClassInitFunc) gl_ellipses_class_init,
      (GtsObjectInitFunc) gl_ellipses_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_scalar_class ()),
				  &gfs_gl_ellipses_info);
  }

  return klass;
}

GtsFile * gfs_gl_ellipses_set (GfsGlEllipses * gl, guint i, const gchar * func)
{
  GtsFile * fp;
  
  g_return_val_if_fail (gl != NULL, NULL);
  g_return_val_if_fail (i < 4, NULL);
  g_return_val_if_fail (func != NULL, NULL);

  if ((fp = gfs_gl_var_func_set (gl->vf[i], GFS_GL (gl)->sim, func, gl->expr[i], NULL)))
    return fp;

  gl->v[i] = gl->vf[i]->v;
  update_ellipses_norm (gl);

  return NULL;
}

/* GfsGlLocation: Object */

static void gl_location_read (GtsObject ** o, GtsFile * fp)
{
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "size",  TRUE},
    {GTS_INT,    "label", TRUE},
    {GTS_NONE}
  };

  (* GTS_OBJECT_CLASS (gfs_gl_location_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  
  var[0].data = &GFS_GL_LOCATION (*o)->size;
  var[1].data = &GFS_GL_LOCATION (*o)->label;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;
}

static void gl_location_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_gl_location_class ())->parent_class->write) (o, fp);

  fprintf (fp, " {\n  size = %g\n  label = %d\n}", 
	   GFS_GL_LOCATION (o)->size, GFS_GL_LOCATION (o)->label);
}

static void gl_location_draw (GfsGl * gl, GfsFrustum * f)
{
  GSList * i;
  gdouble size;

  i = gl->sim->events->items;
  while (i) {
    if (GFS_IS_OUTPUT_LOCATION (i->data)) {
      GfsOutputLocation * l = i->data;
      guint j;

      for (j = 0; j < l->p->len; j++) {
	FttVector p = g_array_index (l->p, FttVector, j);
	gfs_simulation_map (gl->sim, &p);
	FttCell * cell = gfs_domain_locate (GFS_DOMAIN (gl->sim), p, -1, NULL);

	if (gl->maxlevel >= 0)
	  size = GFS_GL_LOCATION (gl)->size*ftt_level_size (gl->maxlevel);
	else if (cell)
	  size = GFS_GL_LOCATION (gl)->size*ftt_cell_size (cell);
	else
	  size = GFS_GL_LOCATION (gl)->size/64.;

	if (cell)
	  glColor3f (gl->lc.r, gl->lc.g, gl->lc.b);
	else
	  glColor3f (1., 0., 0.);

	glMatrixMode (GL_MODELVIEW);
	glPushMatrix ();
	glTranslated (p.x, p.y, p.z);
	glScaled (size, size, size);
	draw_sphere ();
	glPopMatrix ();
	
	if (j == 0 && GFS_GL_LOCATION (gl)->label)
	  gl_draw_text (gl, l->label ? l->label : GFS_OUTPUT (l)->format,
			p.x + size/2., p.y + size/2., p.z + size/2.,
			size*gl->font_size);
      }
    }
    i = i->next;
  }
}

static void gl_location_init (GfsGl * gl)
{
  GtsColor c = { 0., 0., 1. };

  gl->lc = c;
  gl->font_size = 2.;

  GFS_GL_LOCATION (gl)->size = 1.;
}

static gboolean gl_location_relevant (GfsSimulation * sim)
{
  GSList * i = sim->events->items;

  while (i) {
    if (GFS_IS_OUTPUT_LOCATION (i->data))
      return TRUE;
    i = i->next;
  }
  return FALSE;
}

static void gl_location_class_init (GfsGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_location_read;
  GTS_OBJECT_CLASS (klass)->write = gl_location_write;
  klass->draw = gl_location_draw;
  klass->relevant = gl_location_relevant;
}

GfsGlClass * gfs_gl_location_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_location_info = {
      "GfsGlLocation",
      sizeof (GfsGlLocation),
      sizeof (GfsGlClass),
      (GtsObjectClassInitFunc) gl_location_class_init,
      (GtsObjectInitFunc) gl_location_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_class ()), &gfs_gl_location_info);
  }

  return klass;
}

/* GfsGlHeight: Object */

typedef struct {
  GfsGl * gl;
  GfsVariableTracerVOFHeight * t;
} HFData;

static void draw_height (GfsGl * gl, FttComponent c, 
			 FttVector p, gdouble size,
			 gdouble h, const char * label)
{
  glColor3f (gl->lc.r, gl->lc.g, gl->lc.b);
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix ();
  FttVector o = p;
  (&o.x)[c] += h;
  glTranslated (o.x, o.y, o.z);
  glScaled (size, size, size);
  draw_sphere ();
  glPopMatrix ();

  if (GFS_GL_LOCATION (gl)->label)
    gl_draw_text (gl, label,
		  o.x + size/2., o.y + size/2., o.z + size/2.,
		  size*gl->font_size);
  
  gl->size++;
}

static void gl_height (FttCell * cell, HFData * hf)
{
  GfsGl * gl = hf->gl;
  FttVector p;
  gdouble h = ftt_cell_size (cell);
  ftt_cell_pos (cell, &p);
  
  gdouble size = GFS_GL_LOCATION (gl)->size;
  if (gl->maxlevel >= 0)
    size *= ftt_level_size (gl->maxlevel);
  else
    size *= h;

  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++) {
    static gchar hb[3][4] = {"Hbx", "Hby", "Hbz"};
    static gchar ht[3][4] = {"Htx", "Hty", "Htz"};
    if (GFS_HAS_DATA (cell, hf->t->hb[c]))
      draw_height (gl, c, p, size, GFS_VALUE (cell, hf->t->hb[c])*h, hb[c]);
    if (GFS_HAS_DATA (cell, hf->t->ht[c]))
      draw_height (gl, c, p, size, - GFS_VALUE (cell, hf->t->ht[c])*h, ht[c]);
  }
}

static void gl_height_draw (GfsGl * gl, GfsFrustum * f)
{
  HFData hf = { gl, NULL };
  GSList * i = GFS_DOMAIN (gl->sim)->variables;

  while (i && !hf.t) {
    if (GFS_IS_VARIABLE_TRACER_VOF_HEIGHT (i->data))
      hf.t = i->data;
    i = i->next;
  }
  if (hf.t == NULL)
    return;

  gl->size = 0;
  gfs_gl_cell_traverse_visible (gl, f, (FttCellTraverseFunc) gl_height, &hf);
}

static void gl_height_init (GfsGl * gl)
{
  GtsColor c = { 0., 1., 0. };
  gl->lc = c;
  GFS_GL_LOCATION (gl)->size = 0.2;
}

static gboolean gl_height_relevant (GfsSimulation * sim)
{
  GSList * i = GFS_DOMAIN (sim)->variables;

  while (i) {
    if (GFS_IS_VARIABLE_TRACER_VOF_HEIGHT (i->data))
      return TRUE;
    i = i->next;
  }
  return FALSE;
}

static void gl_height_class_init (GfsGlClass * klass)
{
  klass->draw = gl_height_draw;
  klass->relevant = gl_height_relevant;
}

GfsGlClass * gfs_gl_height_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsGlHeight",
      sizeof (GfsGlLocation),
      sizeof (GfsGlClass),
      (GtsObjectClassInitFunc) gl_height_class_init,
      (GtsObjectInitFunc) gl_height_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_location_class ()), &info);
  }

  return klass;
}

/* GfsGlLocate: Object */

static void gl_locate_read (GtsObject ** o, GtsFile * fp)
{
  GfsGlLocate * gl = GFS_GL_LOCATE (*o);
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x", TRUE},
    {GTS_DOUBLE, "y", TRUE},
    {GTS_DOUBLE, "z", TRUE},
    {GTS_NONE}
  };

  (* GTS_OBJECT_CLASS (gfs_gl_locate_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  var[0].data = &gl->p.x;
  var[1].data = &gl->p.y;
  var[2].data = &gl->p.z;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;
}

static void gl_locate_write (GtsObject * o, FILE * fp)
{
  GfsGlLocate * gl = GFS_GL_LOCATE (o);

  (* GTS_OBJECT_CLASS (gfs_gl_locate_class ())->parent_class->write) (o, fp);

  fprintf (fp, " { x = %g y = %g z = %g }", gl->p.x, gl->p.y, gl->p.z);
}

static void gl_locate_draw (GfsGl * gl, GfsFrustum * f)
{
  FttVector p = GFS_GL_LOCATE (gl)->p;
  gfs_simulation_map (gl->sim, &p);
  FttCell * cell = gfs_domain_locate (GFS_DOMAIN (gl->sim), p, gl->maxlevel, NULL);
  if (cell) {
    gl->size = 0;
    glMatrixMode (GL_PROJECTION);
    glPushMatrix ();
    glTranslatef (0., 0., gl->p->lc);
    glNormal3d (0., 0., 1.);
    gl_locate (cell, gl);
    glPopMatrix ();
  }
}

static void gl_locate_class_init (GfsGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gl_locate_read;
  GTS_OBJECT_CLASS (klass)->write = gl_locate_write;  
  klass->draw = gl_locate_draw;
}

GfsGlClass * gfs_gl_locate_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_locate_info = {
      "GfsGlLocate",
      sizeof (GfsGlLocate),
      sizeof (GfsGlClass),
      (GtsObjectClassInitFunc) gl_locate_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_class ()), &gfs_gl_locate_info);
  }

  return klass;
}

/* GfsGlPipes: Object */

static void gl_pipes_draw (GfsGl * gl, GfsFrustum * f)
{
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glTranslatef (0., 0., gl->p->lc);
  glNormal3d (0., 0., 1.);
  glColor3f (gl->lc.r, gl->lc.g, gl->lc.b);
  
  GSList * i = gl->sim->events->items;
  while (i) {
    if (GFS_IS_SOURCE_PIPE (i->data)) {
      GfsSourcePipe * p = i->data;
      FttVector start = p->start, end = p->end;
      gfs_simulation_map (gl->sim, &start);
      gfs_simulation_map (gl->sim, &end);
      FttVector t = { end.x - start.x, end.y - start.y };
      gdouble dt = ftt_vector_norm (&t);
      if (dt > 0.) {
	t.x /= dt; t.y /= dt;
	glMatrixMode (GL_MODELVIEW);
	glPushMatrix ();
	glTranslatef (start.x, start.y, 0.);
	glRotatef (atan2 (t.y, t.x)*180./M_PI, 0., 0., 1.);
	/* fixme: metric is not taken into account properly */
	gdouble h = p->diameter/gl->sim->physical_params.L/2.;
	glBegin (GL_QUADS);
	glVertex2d (dt, h);
	glVertex2d (0., h);
	glVertex2d (0., -h);
	glVertex2d (dt, -h);
	glEnd ();
	if (GFS_EVENT (p)->name) {
	  gchar * label = GFS_EVENT (p)->name;
	  float bbox[6];
	  float size = dt/strlen(label)*gl->font_size;
	  gl_text_bbox (gl, label, size, bbox);
	  float twidth  = bbox[3] - bbox[0];
	  float theight = bbox[4] - bbox[1];	
	  gl_draw_text (gl, label, dt/2. - twidth/2., theight/2., 0., size);
	}
	glPopMatrix ();
	glMatrixMode (GL_PROJECTION);
      }
    }
    i = i->next;
  }

  glPopMatrix ();
}

static void gl_pipes_init (GfsGl * gl)
{
  GtsColor c = { 1., 0., 0. };
  gl->lc = c;
  gl->font_size = 1.;
}

static gboolean gl_pipes_relevant (GfsSimulation * sim)
{
  GSList * i = sim->events->items;

  while (i) {
    if (GFS_IS_SOURCE_PIPE (i->data))
      return TRUE;
    i = i->next;
  }
  return FALSE;
}

static void gl_pipes_class_init (GfsGlClass * klass)
{
  klass->draw = gl_pipes_draw;
  klass->relevant = gl_pipes_relevant;
}

GfsGlClass * gfs_gl_pipes_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsGlPipes",
      sizeof (GfsGl),
      sizeof (GfsGlClass),
      (GtsObjectClassInitFunc) gl_pipes_class_init,
      (GtsObjectInitFunc) gl_pipes_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_class ()), &info);
  }

  return klass;
}

/* GfsGlClipPlane: Object */

static void gl_clip_plane_destroy (GtsObject * o)
{
  gint i = GFS_GL_CLIP_PLANE (o)->i;

  if (i >= 0) {
    glDisable (GL_CLIP_PLANE0 + i);
    GFS_GL (o)->p->cp[i] = FALSE;
  }

  (* GTS_OBJECT_CLASS (gfs_gl_clip_plane_class ())->parent_class->destroy) (o);
}

static void gl_clip_plane_draw (GfsGl * gl, GfsFrustum * f)
{
  if (GFS_GL_CLIP_PLANE (gl)->i >= 0) {
    guint cp = GL_CLIP_PLANE0 + GFS_GL_CLIP_PLANE (gl)->i;
    
    if (GFS_GL_CLIP_PLANE (gl)->disabled)
      glDisable (cp);
    else {
      GfsGl2D * gl2 = GFS_GL2D (gl);
      GLdouble eq[] = {-gl2->n.x, -gl2->n.y, -gl2->n.z, gl2->pos};
      
      glClipPlane (cp, eq);
      glEnable (cp);
    }
  }
}

static void gl_clip_plane_set_simulation (GfsGl * object, GfsSimulation * sim)
{
  (*GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_clip_plane_class ())->parent_class)->set_simulation)
    (object, sim);

  if (GFS_GL_CLIP_PLANE (object)->i < 0) {
    guint i;
  
    g_assert (object->p);
    for (i = 0; i < 6; i++)
      if (!object->p->cp[i]) {
	object->p->cp[i] = TRUE;
	GFS_GL_CLIP_PLANE (object)->i = i;
	return;
      }
    g_warning ("too many clipping planes!");
  }
}

static void gl_clip_plane_class_init (GfsGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gl_clip_plane_destroy;
  klass->set_simulation = gl_clip_plane_set_simulation;
  klass->draw = gl_clip_plane_draw;
  klass->pick = NULL;
}

static void gl_clip_plane_init (GfsGlClipPlane * cp)
{
  cp->i = -1;
}

GfsGlClass * gfs_gl_clip_plane_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_clip_plane_info = {
      "GfsGlClipPlane",
      sizeof (GfsGlClipPlane),
      sizeof (GfsGl2DClass),
      (GtsObjectClassInitFunc) gl_clip_plane_class_init,
      (GtsObjectInitFunc) gl_clip_plane_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl2D_class ()), &gfs_gl_clip_plane_info);
  }

  return klass;
}

void gfs_gl_clip_plane_disable (GfsGlClipPlane * gl)
{
  gboolean disabled = gl->disabled;
  gl->disabled = TRUE;
  gfs_gl_draw (GFS_GL (gl), NULL);
  gl->disabled = disabled;
}

/* GfsGlCutPlane: Object */

static void gl_cut_plane (FttCell * cell, GfsGl * gl)
{
  GList * i = GFS_GL_CUT_PLANE (gl)->list;
  while (i) {
    if (GFS_GL_CLASS (GTS_OBJECT (i->data)->klass)->cut)
      (* GFS_GL_CLASS (GTS_OBJECT (i->data)->klass)->cut) (i->data, cell, GFS_GL2D (gl));
    i = i->next;
  }
}

static void gfs_gl_cut_plane_draw (GfsGl * gl, GfsFrustum * f)
{
  gl->size = 0;
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glTranslatef (0., 0., gl->p->lc);
  glBegin (GL_LINES);
  gfs_gl_normal (gl);
  gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc) gl_cut_plane, gl);
  glEnd ();
  glPopMatrix ();
}

static void gl_cut_plane_class_init (GfsGlClass * klass)
{
  klass->draw = gfs_gl_cut_plane_draw;
  klass->pick = NULL;
}

GfsGlClass * gfs_gl_cut_plane_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_cut_plane_info = {
      "GfsGlCutPlane",
      sizeof (GfsGlCutPlane),
      sizeof (GfsGl2DClass),
      (GtsObjectClassInitFunc) gl_cut_plane_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl2D_class ()), &gfs_gl_cut_plane_info);
  }

  return klass;
}

/* GfsGlViewParams */

void gfs_gl_view_params_init (GfsGlViewParams * p)
{
  g_return_if_fail (p != NULL);

  p->do_init = TRUE;
  p->beginx = p->beginy = 0;
  p->dx = p->dy = 0;
  p->tx = p->ty = 0;
  p->sx = p->sy = p->sz = 1.;
  p->quat[0] = p->quat[1] = p->quat[2] = 0; p->quat[3] = 1;
  p->dquat[0] = p->dquat[1] = p->dquat[2] = 0; p->dquat[3] = 1;
  p->fov = 30.;
  gfs_gl_trackball (p->quat , 0.0, 0.0, 0.0, 0.0);

  p->bg.r = 0.3; p->bg.g = 0.4; p->bg.b = 0.6;
  p->res = p->base_res = 1.;
  p->lc = 0.001;
  p->reactivity = 0.1;
  p->motion = FALSE;

  guint i;
  for (i = 0; i < 6; i++)
    p->cp[i] = FALSE;
}

GfsGlViewParams * gfs_gl_view_params_new (void)
{
  GfsGlViewParams * p = g_malloc (sizeof (GfsGlViewParams));
  gfs_gl_view_params_init (p);
  return p;
}

void gfs_gl_view_params_write (GfsGlViewParams * p, FILE * fp)
{
  g_return_if_fail (p != NULL);
  g_return_if_fail (fp != NULL);
  
  fprintf (fp, 
	   "View {\n"
	   "  tx = %g ty = %g\n"
	   "  sx = %g sy = %g sz = %g\n"
	   "  q0 = %g q1 = %g q2 = %g q3 = %g\n"
	   "  fov = %g\n"
	   "  r = %g g = %g b = %g\n"
	   "  res = %g\n"
	   "  lc = %g\n"
	   "  reactivity = %g\n"
	   "}",
	   p->tx, p->ty,
	   p->sx, p->sy, p->sz,
	   p->quat[0], p->quat[1], p->quat[2], p->quat[3], 
	   p->fov,
	   p->bg.r, p->bg.g, p->bg.b,
	   p->base_res,
	   p->lc,
	   p->reactivity);
}

void gfs_gl_view_params_read (GfsGlViewParams * p, GtsFile * fp)
{
  GtsFileVariable var[] = {
    {GTS_FLOAT, "tx",         TRUE},
    {GTS_FLOAT, "ty",         TRUE},
    {GTS_FLOAT, "q0",         TRUE},
    {GTS_FLOAT, "q1",         TRUE},
    {GTS_FLOAT, "q2",         TRUE},
    {GTS_FLOAT, "q3",         TRUE},
    {GTS_FLOAT, "fov",        TRUE},
    {GTS_FLOAT, "r",          TRUE},
    {GTS_FLOAT, "g",          TRUE},
    {GTS_FLOAT, "b",          TRUE},
    {GTS_FLOAT, "res",        TRUE},
    {GTS_FLOAT, "lc",         TRUE},
    {GTS_FLOAT, "reactivity", TRUE},
    {GTS_FLOAT, "sx",         TRUE},
    {GTS_FLOAT, "sy",         TRUE},
    {GTS_FLOAT, "sz",         TRUE},
    {GTS_NONE}
  };

  g_return_if_fail (p != NULL);
  g_return_if_fail (fp != NULL);
  
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (\"View\")");
    return;
  }
  if (strcmp (fp->token->str, "View")) {
    gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  var[0].data = &p->tx;
  var[1].data = &p->ty;
  var[2].data = &p->quat[0];
  var[3].data = &p->quat[1];
  var[4].data = &p->quat[2];
  var[5].data = &p->quat[3];
  var[6].data = &p->fov;
  var[7].data = &p->bg.r;
  var[8].data = &p->bg.g;
  var[9].data = &p->bg.b;
  var[10].data = &p->base_res;
  var[11].data = &p->lc;
  var[12].data = &p->reactivity;
  var[13].data = &p->sx;
  var[14].data = &p->sy;
  var[15].data = &p->sz;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;
}

void gfs_gl2ps_params_init (GfsGl2PSParams * p)
{
  g_return_if_fail (p != NULL);

  p->format = GFSGL_PPM;
  p->options = (GL2PS_SIMPLE_LINE_OFFSET |
		GL2PS_SILENT |
		GL2PS_BEST_ROOT |
		GL2PS_OCCLUSION_CULL |
		GL2PS_USE_CURRENT_VIEWPORT |
		GL2PS_TIGHT_BOUNDING_BOX);
  p->width = p->height = 0;
  p->sort = GL2PS_SIMPLE_SORT;
  p->lw = 1.;
}

void gfs_gl2ps_params_read (GfsGl2PSParams * p, GtsFile * fp)
{
  g_return_if_fail (p != NULL);
  g_return_if_fail (fp != NULL);

  gchar * format, * orientation, * sort;
  GtsFileVariable var[] = {
    {GTS_STRING, "format",      TRUE, &format},
    {GTS_STRING, "orientation", TRUE, &orientation},
    {GTS_FLOAT,  "line_width",  TRUE, &p->lw},
    {GTS_UINT,   "width",       TRUE, &p->width},
    {GTS_UINT,   "height",      TRUE, &p->height},
    {GTS_STRING, "sort",        TRUE, &sort},
    {GTS_NONE}
  };

  gfs_gl2ps_params_init (p);
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (var[0].set) {
    if (!strcmp (format, "PPM"))
      p->format = GFSGL_PPM;
    else if (!strcmp (format, "PS"))
      p->format = GL2PS_PS;
    else if (!strcmp (format, "EPS"))
      p->format = GL2PS_EPS;
    else if (!strcmp (format, "PDF"))
      p->format = GL2PS_PDF;
    else if (!strcmp (format, "SVG"))
      p->format = GL2PS_SVG;
    else if (!strcmp (format, "Gnuplot"))
      p->format = GFSGL_GNUPLOT;
    else if (!strcmp (format, "OBJ"))
      p->format = GFSGL_OBJ;
    else if (!strcmp (format, "KML"))
      p->format = GFSGL_KML;
    else if (!strcmp (format, "Latex"))
      p->format = GL2PS_TEX;
    else {
      gts_file_variable_error (fp, var, "format", "unknown format `%s'", format);
      g_free (format);
      return;
    }
    g_free (format);
  }

  if (var[1].set) {
    if (!strcmp (orientation, "Portrait"))
      p->options &= ~GL2PS_LANDSCAPE;
    else if (!strcmp (orientation, "Landscape"))
      p->options |= GL2PS_LANDSCAPE;
    else {
      gts_file_variable_error (fp, var, "orientation", "unknown orientation `%s'", fp->token->str);
      g_free (orientation);
      return;
    }
    g_free (orientation);
  }
  else
    p->options &= ~GL2PS_LANDSCAPE;  

  if (var[5].set) {
    if (!strcmp (sort, "BSP"))
      p->sort = GL2PS_BSP_SORT;
    else if (!strcmp (sort, "Simple"))
      p->sort = GL2PS_SIMPLE_SORT;
    else if (!strcmp (sort, "None"))
      p->sort = GL2PS_NO_SORT;
    else {
      gts_file_variable_error (fp, var, "sort", "unknown option `%s'", sort);
      g_free (sort);
      return;
    }
    g_free (sort);
  }
}

void gfs_gl2ps_params_write (GfsGl2PSParams * p, FILE * fp)
{
  g_return_if_fail (p != NULL);
  g_return_if_fail (fp != NULL);

  fprintf (fp, " { format = %s orientation = %s line_width = %g width = %d height = %d sort = %s }",
	   p->format == GFSGL_PPM ?     "PPM" :
	   p->format == GL2PS_PS ?      "PS" :
	   p->format == GL2PS_EPS ?     "EPS" :
	   p->format == GL2PS_PDF ?     "PDF" :
	   p->format == GL2PS_SVG ?     "SVG" :
	   p->format == GFSGL_GNUPLOT ? "Gnuplot" :
	   p->format == GFSGL_OBJ ?     "OBJ" :
	   p->format == GFSGL_KML ?     "KML" :
	   p->format == GL2PS_TEX ?     "Latex" : 
	   "?",
	   p->options & GL2PS_LANDSCAPE ? "Landscape" : "Portrait",
	   p->lw, p->width, p->height,
	   p->sort == GL2PS_SIMPLE_SORT ? "Simple" :
	   p->sort == GL2PS_BSP_SORT ?    "BSP" :
	   p->sort == GL2PS_NO_SORT ?     "None" :
	   "?");
}

void gfs_gl_write_image (FILE * fp, const GLubyte * buffer, guint width, guint height)
{
  gint i, x, y;
  const GLubyte *ptr = buffer;

  g_return_if_fail (fp != NULL);
  g_return_if_fail (buffer != NULL);

  fprintf (fp, "P6 %d %d 255\n", width, height);
  for (y = height - 1; y >= 0; y--) {
    for (x = 0; x < width; x++) {
      i = (y*width + x)*4;
      fputc (ptr[i], fp);   /* write red */
      fputc (ptr[i+1], fp); /* write green */
      fputc (ptr[i+2], fp); /* write blue */
    }
  }
}

typedef struct {
  GList * symmetries;
  FttVector * s;
  gdouble max;
} ExtentData;

static void update_max (FttVector * p, ExtentData * d)
{
  guint i, n = create_symmetries (d->s, d->symmetries, p);
  for (i = 0; i < n; i++) {
    FttComponent c;
    gdouble m = 0.;
    for (c = 0; c < FTT_DIMENSION; c++)
      m += (&d->s[i].x)[c]*(&d->s[i].x)[c];
    if (m > d->max) d->max = m;
  }
}

static void extent (GfsBox * b, ExtentData * d)
{
  FttCell * cell = b->root;
  gdouble h = ftt_cell_size (cell)/2.;
  FttVector p, o;

  ftt_cell_pos (cell, &p);
#if FTT_2D
  o.x = p.x + h; o.y = p.y + h; update_max (&o, d);
  o.x = p.x + h; o.y = p.y - h; update_max (&o, d);
  o.x = p.x - h; o.y = p.y + h; update_max (&o, d);
  o.x = p.x - h; o.y = p.y - h; update_max (&o, d);
#else /* 3D */
  o.x = p.x + h; o.y = p.y + h; o.z = p.z + h; update_max (&o, d);
  o.x = p.x + h; o.y = p.y - h; o.z = p.z + h; update_max (&o, d);
  o.x = p.x - h; o.y = p.y + h; o.z = p.z + h; update_max (&o, d);
  o.x = p.x - h; o.y = p.y - h; o.z = p.z + h; update_max (&o, d);
  o.x = p.x + h; o.y = p.y + h; o.z = p.z - h; update_max (&o, d);
  o.x = p.x + h; o.y = p.y - h; o.z = p.z - h; update_max (&o, d);
  o.x = p.x - h; o.y = p.y + h; o.z = p.z - h; update_max (&o, d);
  o.x = p.x - h; o.y = p.y - h; o.z = p.z - h; update_max (&o, d);
#endif /* 3D */
}

/**
 * gfs_gl_domain_extent:
 * @domain: a #GfsDomain.
 * @symmetries: a list of #GfsGlSymmetry.
 *
 * Returns: the maximum extent along any dimension of @domain.
 */
gdouble gfs_gl_domain_extent (GfsDomain * domain, GList * symmetries)
{
  if (domain == NULL)
    return 1.;

  ExtentData d;
  d.max = 0.;
  d.symmetries = symmetries;
  guint n = 1;
  while (symmetries) {
    n *= 2;
    symmetries = symmetries->next;
  }
  d.s = g_malloc (n*sizeof (FttVector));
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) extent, &d);
  g_free (d.s);
  gfs_all_reduce (domain, d.max, MPI_DOUBLE, MPI_MAX);
  return 2.*sqrt (d.max);
}

/**
 * gfs_gl_feedback_begin:
 * @buffersize: the size of the feedback buffer.
 *
 * Returns: a newly setup openGL feedback buffer.
 */
GfsGlFeedback * gfs_gl_feedback_begin (guint buffersize)
{
  GfsGlFeedback * f;

  g_return_val_if_fail (buffersize > 0, NULL);

  f = g_malloc (sizeof (GfsGlFeedback));
  f->feedback = g_malloc (sizeof (GLfloat)*buffersize);
  glFeedbackBuffer (buffersize, GL_3D_COLOR, f->feedback);
  glRenderMode (GL_FEEDBACK);
  return f;
}

/**
 * gfs_gl_feedback_end:
 * @f: the #GfsGlFeedback.
 * @sim: a #GfsSimulation.
 * @fp: a file pointer.
 * @format: the file format.
 *
 * Writes to @fp the contents of @f.
 *
 * Returns: %FALSE if @f has overflowed, %TRUE otherwise.
 */
gboolean gfs_gl_feedback_end (GfsGlFeedback * f, GfsSimulation * sim, FILE * fp, GfsGlFormat format)
{
  GLint size;

  g_return_val_if_fail (f != NULL, FALSE);
  g_return_val_if_fail (sim != NULL, FALSE);
  g_return_val_if_fail (fp != NULL, FALSE);

  size = glRenderMode (GL_RENDER);
  if (size >= 0) {
    GLint used = size;
    GLfloat * current = f->feedback;
    GLdouble model[16], proj[16];
    GLint viewport[4];
    FttVector p, lastp = { G_MAXDOUBLE, G_MAXDOUBLE, G_MAXDOUBLE };
    guint nps = 0, np = 0, linetoken = FALSE;

    /* Header */
    switch (format) {
    case GFSGL_KML:
      fputs ("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
	     "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n"
	     "<Placemark>\n"
	     "<name>GfsView</name>\n"
	     "<MultiGeometry>\n", fp);
      break;
    default:
      break;
    }

    /* Body */
    glGetDoublev (GL_MODELVIEW_MATRIX, model);
    glGetDoublev (GL_PROJECTION_MATRIX, proj);
    glGetIntegerv (GL_VIEWPORT, viewport);
    while (used > 0)
      switch ((GLint) *current) {

      case GL_POINT_TOKEN :
	current++; used--;
	gluUnProject (current[0], current[1], current[2], model, proj, viewport, &p.x, &p.y, &p.z);
	gfs_simulation_map_inverse (sim, &p);
	switch (format) {
	case GFSGL_GNUPLOT:
	  fprintf (fp, "%g %g %g\n\n", p.x, p.y, p.z); break;
	case GFSGL_OBJ:
	  fprintf (fp, "v %g %g %g\np -1\n", p.x, p.y, p.z); 
	  np++;
	  nps = 0;
	  break;
	case GFSGL_KML:
	  if (linetoken) {
	    fputs ("</coordinates>\n</LineString>\n", fp);
	    linetoken = FALSE;
	  }
	  fprintf (fp, "<Point>\n<coordinates>%f,%f,%f</coordinates>\n</Point>\n", p.x, p.y, p.z);
	  break;
	default:
	  g_assert_not_reached ();
	}
	current += 7; used -= 7;
	break;

      case GL_LINE_TOKEN :
      case GL_LINE_RESET_TOKEN :
	current++; used--;
	gluUnProject (current[0], current[1], current[2], model, proj, viewport, &p.x, &p.y, &p.z);
	gfs_simulation_map_inverse (sim, &p);
	if (p.x != lastp.x || p.y != lastp.y || p.z != lastp.z) {
	  switch (format) {
	  case GFSGL_GNUPLOT:
	    fprintf (fp, "\n%g %g %g\n", p.x, p.y, p.z);
	    break;
	  case GFSGL_OBJ:
	    if (nps > 0) {
	      fputc ('l', fp);
	      int i;
	      for (i = nps; i <= np; i++)
		fprintf (fp, " %d", i);
	      fputc ('\n', fp);
	    }
	    nps = ++np;
	    fprintf (fp, "v %g %g %g\n", p.x, p.y, p.z);
	    break;
	  case GFSGL_KML:
	    if (linetoken)
	      fputs ("</coordinates>\n</LineString>\n", fp);
	    fprintf (fp, "<LineString>\n<coordinates>\n%f,%f,%f\n", p.x, p.y, p.z);
	    linetoken = TRUE;
	    break;
	  default:
	    g_assert_not_reached ();
	  }
	}
	current += 7; used -= 7;
	gluUnProject (current[0], current[1], current[2], model, proj, viewport, &p.x, &p.y, &p.z);
	gfs_simulation_map_inverse (sim, &p);
	lastp = p;
	switch (format) {
	case GFSGL_GNUPLOT:
	  fprintf (fp, "%g %g %g\n", p.x, p.y, p.z); 
	  break;
	case GFSGL_OBJ:
	  fprintf (fp, "v %g %g %g\n", p.x, p.y, p.z); 
	  np++;
	  break;
	case GFSGL_KML:
	  fprintf (fp, "%f,%f,%f\n", p.x, p.y, p.z);
	  break;
	default:
	  g_assert_not_reached ();
	}
	current += 7; used -= 7;
	break;

      case GL_POLYGON_TOKEN : {
	GLint count = (GLint) current[1], vcount = 0;
	current += 2;
	used -= 2;
	if (format == GFSGL_KML) {
	  if (linetoken) {
	    fputs ("</coordinates>\n</LineString>\n", fp);
	    linetoken = FALSE;
	  }
	  fputs ("<Polygon>\n<coordinates>\n", fp);
	}
	while (count > 0 && used > 0) {
	  gluUnProject (current[0], current[1], current[2], model, proj, viewport, 
			&p.x, &p.y, &p.z);
	  gfs_simulation_map_inverse (sim, &p);
	  switch (format) {
	  case GFSGL_GNUPLOT:
	    fprintf (fp, "%g %g %g\n", p.x, p.y, p.z); break;
	  case GFSGL_OBJ:
	    fprintf (fp, "v %g %g %g\n", p.x, p.y, p.z); 
	    np++;
	    nps = 0;
	    break;
	  case GFSGL_KML:
	    fprintf (fp, "%f,%f,%f\n", p.x, p.y, p.z); break;
	  default:
	    g_assert_not_reached ();
	  }
	  current += 7; used -= 7;
	  count--;
	  vcount++;
	}
	switch (format) {
	case GFSGL_KML:
	  fputs ("</coordinates>\n</Polygon>\n", fp); break;
	case GFSGL_OBJ:
	  fprintf (fp, "f -%d", vcount--);
	  while (vcount)
	    fprintf (fp, " -%d", vcount--);
	  /* fall through */
	case GFSGL_GNUPLOT: fputc ('\n', fp); break;
	default:
	  g_assert_not_reached ();
	}
	break;
      }
      }

    /* Footer */
    switch (format) {
    case GFSGL_OBJ:
      if (np > nps && nps > 0) {
	fputc ('l', fp);
	int i;
	for (i = nps; i <= np; i++)
	  fprintf (fp, " %d", i);
	fputc ('\n', fp);
      }
      break;
     
    case GFSGL_KML:
      if (linetoken) {
	fputs ("</coordinates>\n</LineString>\n", fp);
	linetoken = FALSE;
      }
      fputs ("</MultiGeometry>\n"
	     "</Placemark>\n"
	     "</kml>\n", fp);
      break;
    
    default:
      break;
    }
  }
  g_free (f->feedback);
  g_free (f);
  return (size >= 0);
}
