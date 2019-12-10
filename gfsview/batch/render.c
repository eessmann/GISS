/* Gerris - The GNU Flow Solver
 * Copyright (C) 2010 National Institute of Water and
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

#include <string.h>
#include <GL/osmesa.h>
#if defined(__APPLE__)
#  include <OpenGL/glu.h>
#else
#  include <GL/glu.h>
#endif

#include "config.h"
#ifdef HAVE_MPI
#  include <mpi.h>
#endif /* HAVE_MPI */

#include "render.h"
#include "gl/trackball.h"

static GList * get_symmetries (GList * i)
{
  GList * symmetry = NULL;

  while (i) {
    if (GFS_IS_GL_SYMMETRY (i->data))
      symmetry = g_list_append (symmetry, i->data);
    i = i->next;
  }
  return symmetry;
}

static void view_draw (GfsGlViewParams * view,
		       GfsGl2PSParams * p,
		       GfsDomain * domain,
		       GList * list,
		       guint width, guint height)
{
  GLfloat m[4][4];
  gdouble max;
  GList * i;

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  GList * symmetries = get_symmetries (list);
  max = gfs_gl_domain_extent (domain, symmetries);
  gluPerspective (view->fov, width/(float)height, 1., 1. + 2.*max);
  glMatrixMode (GL_MODELVIEW);
	    
  glLoadIdentity ();
  glTranslatef (view->tx, view->ty, - (1. + max));
  gfs_gl_build_rotmatrix (m, view->quat);
  glMultMatrixf (&m[0][0]);
  glScalef (view->sx, view->sy, view->sz);
	    
  glClearColor (view->bg.r, view->bg.g, view->bg.b, 0.);
  glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  GfsFrustum frustum;
  gfs_gl_get_frustum (view, symmetries, &frustum);
  GLuint display_list = glGenLists (1);

  glNewList (display_list, GL_COMPILE);
  i = list;
  while (i) {
    GfsGl * gl = i->data;

    if (GFS_IS_GL_CLIP_PLANE (gl)) {
      gl->format = p->format;
      gfs_gl_clip_plane_disable (GFS_GL_CLIP_PLANE (gl));
    }
    i = i->next;
  }
  i = list;
  while (i) {
    if (GFS_IS_GL_CUT_PLANE (i->data))
      GFS_GL_CUT_PLANE (i->data)->list = list;
    i = i->next;
  }

  GSList * clip = NULL;
  gboolean firstclip = TRUE;
  i = list;
  while (i) {
    GfsGl * gl = i->data;
    gl->format = p->format;
    if (GFS_IS_GL_CLIP_PLANE (gl)) {
      if (firstclip) {
	g_slist_foreach (clip, (GFunc) gfs_gl_clip_plane_disable, NULL);
	g_slist_free (clip); clip = NULL;
	firstclip = FALSE;	  
      }
      gfs_gl_draw (gl, &frustum);
      clip = g_slist_prepend (clip, gl);
    }
    else {
      gfs_gl_draw (gl, &frustum);
      firstclip = TRUE;
    }
    i = i->next;
  }  
  g_slist_free (clip);
  glEndList();

  gfs_gl_symmetry_apply (symmetries, display_list);
  gfs_gl_frustum_free (&frustum);
  g_list_free (symmetries);
  glDeleteLists (display_list, 1);
  glFinish ();
}

#ifdef HAVE_MPI
#if FTT_2D
static void compose_image (GLubyte * bg, GLubyte * fg,
			   guint width, guint height)
{
  int i;
  for (i = 0; i < 4*width*height; i += 4)
    if (bg[i+3] == 0) {
      bg[i] = fg[i]; bg[i+1] = fg[i+1]; bg[i+2] = fg[i+2]; bg[i+3] = fg[i+3];
    }
}
#else /* 3D */
static void compose_image (GLubyte * bg, guint32 * bgdepth,
			   GLubyte * fg, guint32 * fgdepth,
			   guint width, guint height)
{
  int i;
  for (i = 0; i < 4*width*height; i += 4)
    if (bgdepth[i/4] > fgdepth[i/4]) {
      bg[i] = fg[i]; bg[i+1] = fg[i+1]; bg[i+2] = fg[i+2]; bg[i+3] = fg[i+3];
      bgdepth[i/4] = fgdepth[i/4];
    }
}
#endif /* 3D */
#endif /* HAVE_MPI */

void gfs_gl_osmesa_render (GfsGl2PSParams * p, GfsSimulation * sim,
			   GfsGlViewParams * view, GList * list,
			   FILE * fptr,
			   gboolean parallel)
{
  OSMesaContext ctx;
  guint width = p->width > 0 ? p->width : 640;
  guint height = p->height > 0 ? p->height : 480;
  void * image = g_malloc (width*height*4*sizeof (GLubyte));

  /* OSMesa somehow generates floating-point exceptions... turn them
     off so that people don't blame gfsview! */
  gfs_disable_floating_point_exceptions ();
	    
  /* Create an RGBA-mode context for OSMesa */
#if OSMESA_MAJOR_VERSION * 100 + OSMESA_MINOR_VERSION >= 305
  /* specify Z, stencil, accum sizes */
  ctx = OSMesaCreateContextExt (OSMESA_RGBA, 32, 0, 0, NULL);
#else
  ctx = OSMesaCreateContext (OSMESA_RGBA, NULL);
#endif
  if (!ctx) {
    fprintf (stderr, "gfsview-batch: OSMesaCreateContext failed!\n");
    exit (1);
  }

  if (!OSMesaMakeCurrent (ctx, image, GL_UNSIGNED_BYTE, width, height)) {
    fprintf (stderr, "gfsview-batch: OSMesaMakeCurrent failed!\n");
    exit (1);
  }
	
  gfs_gl_init_gl ();

  if (sim) {
    switch (p->format) {

    case GFSGL_PPM_OFFSCREEN: case GFSGL_PPM_SCREEN: {
      view_draw (view, p, GFS_DOMAIN (sim), list, width, height);
#ifdef HAVE_MPI
      if (parallel && GFS_DOMAIN (sim)->pid >= 0) {
	int size = width*height*4;
#if FTT_2D
	if (GFS_DOMAIN (sim)->pid == 0) {
	  void * image1 = g_malloc (size);
	  int pe, npe;
	  MPI_Comm_size (MPI_COMM_WORLD, &npe);
	  for (pe = 1; pe < npe; pe++) {
	    MPI_Status status;
	    MPI_Recv (image1, size, MPI_BYTE, pe, 0, MPI_COMM_WORLD, &status);
	    compose_image (image, image1, width, height);
	  }
	  g_free (image1);
	} else
	  MPI_Send (image, size, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
#else /* 3D */
	void * depth1;
	GLint width1, height1, bytesPerValue;
	OSMesaGetDepthBuffer (ctx, &width1, &height1, &bytesPerValue, &depth1);
	g_assert (width == width1 && height == height1 && bytesPerValue == 4);
	if (GFS_DOMAIN (sim)->pid == 0) {
	  void * depth = g_malloc (size);
	  memcpy (depth, depth1, size);
	  depth1 = g_malloc (size);
	  void * image1 = g_malloc (size);
	  int pe, npe;
	  MPI_Comm_size (MPI_COMM_WORLD, &npe);
	  for (pe = 1; pe < npe; pe++) {
	    MPI_Status status;
	    MPI_Recv (image1, size, MPI_BYTE, pe, 0, MPI_COMM_WORLD, &status);
	    MPI_Recv (depth1, size, MPI_BYTE, pe, 0, MPI_COMM_WORLD, &status);
	    compose_image (image, depth, image1, depth1, width, height);
	  }
	  g_free (image1);
	  g_free (depth1);
	  g_free (depth);
	} else {
	  MPI_Send (image,  size, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
	  MPI_Send (depth1, size, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
	}
#endif /* 3D */
      }
#endif /* HAVE_MPI */
      gfs_gl_write_image (fptr, image, width, height);
      break;
    }

    case GFSGL_GNUPLOT: case GFSGL_OBJ: case GFSGL_KML: {
      guint buffsize = 0;
      gboolean done = FALSE;
      gfloat res = view->res;
      view->res = 0.;
      while (!done) {
	GfsGlFeedback * f;
	buffsize += 2048*2048;
	f = gfs_gl_feedback_begin (buffsize);
	view_draw (view, p, GFS_DOMAIN (sim), list, width, height);
	done = gfs_gl_feedback_end (f, sim, fptr, p->format);
      }
      view->res = res;
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
			buffsize, fptr, "");
	view->lw = p->lw;
	view_draw (view, p, GFS_DOMAIN (sim), list, width, height);
	state = gl2psEndPage();
      }
    }
    }
  }
	    
  g_free (image);
  fflush (fptr);

  OSMesaDestroyContext (ctx);

  gfs_enable_floating_point_exceptions ();
}
