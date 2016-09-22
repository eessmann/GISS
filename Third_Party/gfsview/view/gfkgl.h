/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2004 National * Institute of Water and
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

#ifndef __GFKGL_H__
#define __GFKGL_H__

#include <gfs.h>
#include <gtk/gtk.h>

#include "gl/gfsgl.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define GFK_IS_EDITED(o) (g_object_get_data (G_OBJECT (o), "edited"))

G_LOCK_EXTERN (gfk_gl_scripting);
extern gboolean gfk_gl_scripting;

G_LOCK_EXTERN (scripting_pending);

typedef enum {
  GFS_SAVE_EVENT, GFS_APPEND_EVENT, GFS_ECHO_EVENT
} GfkScriptingEvent;

typedef struct {
  GfkScriptingEvent event;
  GtkWidget * view;
  gpointer data;
} GfkScriptingMessage;

gboolean gfk_receive_scripting_message (gpointer data);

enum {
  VISIBLE_COLUMN,
  ICON_COLUMN,
  PROPERTIES_COLUMN,
  GL_COLUMN,
  SELECTED_COLUMN,
  N_COLUMNS
};

guint         gfk_decimal_digits           (double x, guint significant);

/* GfkGl: Header */

typedef struct _GfkGl           GfkGl;

struct _GfkGl {
  /*< private >*/
  GtsObject parent;
  gchar * props;

  /*< public >*/
  GfsGl * gl;
  GtkWidget * glarea, * list;
  GtkWidget * params, * properties, * color_selector, * font;
};

typedef struct _GfkGlClass    GfkGlClass;

struct _GfkGlClass {
  /*< private >*/
  GtsObjectClass parent_class;
  GfsGlClass * gl_class;

  /*< public >*/
  void          (* post_init)        (GfkGl *);
  void          (* set_simulation)   (GfkGl *, GfsSimulation *);
  void          (* update_interface) (GfkGl *);
  GtkWidget   * (* icon)             (GfkGlClass *);
  gchar       * (* name)             (GfkGlClass *);
  gchar       * (* properties)       (GfkGl *);
  gchar       * (* pickinfo)         (GfkGl *, gboolean);
};

#define GFK_GL(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGl,\
					         gfk_gl_class ())
#define GFK_GL_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfkGlClass,\
						 gfk_gl_class())
#define GFK_IS_GL(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_class ()))

GfkGlClass *  gfk_gl_class                 (void);
GfkGl *       gfk_gl_new                   (GfkGlClass * klass,
					    GtkWidget * glarea,
					    GtkWidget * list);
void          gfk_gl_expose                (GfkGl * gl);
void          gfk_gl_set_sensitive         (GfkGl * gl, 
					    GtkWidget * page, 
					    gboolean sensitive);
void          gfk_gl_prepend_params        (GfkGl * gl, 
					    GtkWidget * widget, 
					    GtkWidget * label);
void          gfk_gl_set_color             (GfkGl * gl, 
					    GtsColor c);
void          gfk_gl_update_properties     (GfkGl * gl);
void          gfk_gl_update_interface      (GfkGl * gl);
void          gfk_gl_set_simulation        (GfkGl * gl, 
					    GfsSimulation * sim);

/* GfkGlLabel: Header */

typedef struct _GfkGlLabel         GfkGlLabel;

struct _GfkGlLabel {
  /*< private >*/
  GfkGl parent;
  GtkWidget * label;
};

#define GFK_GL_LABEL(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlLabel,\
					         gfk_gl_label_class ())
#define GFK_IS_GL_LABEL(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_label_class ()))

GfkGlClass * gfk_gl_label_class  (void);

/* GfkGl2D: Header */

typedef struct _GfkGl2D         GfkGl2D;

struct _GfkGl2D {
  /*< private >*/
  GfkGl parent;
  FttVector n;
  GtkWidget * params;
  gchar * pickinfo;
};

#define GFK_GL2D(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGl2D,\
					         gfk_gl2D_class ())
#define GFK_IS_GL2D(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl2D_class ()))

GfkGlClass *  gfk_gl2D_class               (void);
void          gfk_gl2D_update_pos_bounds   (GfkGl2D * gl);

/* GfkGlSymmetry: Header */

#define GFK_IS_GL_SYMMETRY(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_symmetry_class ()))

GfkGlClass *  gfk_gl_symmetry_class             (void);

/* GfkGlPeriodic: Header */

#define GFK_IS_GL_PERIODIC(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_periodic_class ()))

GfkGlClass *  gfk_gl_periodic_class             (void);

/* GfkGlCells: Header */

typedef struct _GfkGlCells         GfkGlCells;

struct _GfkGlCells {
  /*< private >*/
  GfkGl2D parent;
  GtkWidget * cells;
  gboolean edit;
  guint level;
};

#define GFK_GL_CELLS(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlCells,\
					         gfk_gl_cells_class ())
#define GFK_IS_GL_CELLS(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_cells_class ()))

GfkGlClass * gfk_gl_cells_class  (void);

/* GfkGlFractions: Header */

#define GFK_IS_GL_FRACTIONS(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_fractions_class ()))

GfkGlClass * gfk_gl_fractions_class  (void);

/* GfkGlBoundaries: Header */

#define GFK_IS_GL_BOUNDARIES(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_boundaries_class ()))

GfkGlClass * gfk_gl_boundaries_class  (void);

/* GfkGlScalar: Header */

typedef struct _GfkGlScalar         GfkGlScalar;

struct _GfkGlScalar {
  /*< private >*/
  GfkGl2D parent;
  GtkWidget * scalar;
};

#define GFK_GL_SCALAR(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlScalar,\
					         gfk_gl_scalar_class ())
#define GFK_IS_GL_SCALAR(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_scalar_class ()))

GfkGlClass *  gfk_gl_scalar_class            (void);

/* GfkGlSquares: Header */

typedef struct _GfkGlSquares         GfkGlSquares;

struct _GfkGlSquares {
  /*< private >*/
  GfkGlScalar parent;

  /*< public >*/
};

#define GFK_GL_SQUARES(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlSquares,\
					         gfk_gl_squares_class ())
#define GFK_IS_GL_SQUARES(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_squares_class ()))

GfkGlClass * gfk_gl_squares_class  (void);

/* GfkGlLinear: Header */

typedef struct _GfkGlLinear         GfkGlLinear;

struct _GfkGlLinear {
  /*< private >*/
  GfkGlScalar parent;
  GtkWidget * scalar;

  /*< public >*/
};

#define GFK_GL_LINEAR(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlLinear,\
					         gfk_gl_linear_class ())
#define GFK_IS_GL_LINEAR(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_linear_class ()))

GfkGlClass * gfk_gl_linear_class  (void);

/* GfkGlIsoline: Header */

typedef struct _GfkGlIsoline         GfkGlIsoline;

struct _GfkGlIsoline {
  /*< private >*/
  GfkGlLinear parent;
};

#define GFK_GL_ISOLINE(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlIsoline,\
					         gfk_gl_isoline_class ())
#define GFK_IS_GL_ISOLINE(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_isoline_class ()))

GfkGlClass * gfk_gl_isoline_class  (void);

/* GfkGlSolid: Header */

#if (!FTT_2D)
typedef struct _GfkGlSolid         GfkGlSolid;

struct _GfkGlSolid {
  /*< private >*/
  GfkGlScalar parent;
  GtkWidget * solid;
};

#define GFK_GL_SOLID(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlSolid,\
					         gfk_gl_solid_class ())
#endif /* 3D */

#define GFK_IS_GL_SOLID(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_solid_class ()))

GfkGlClass * gfk_gl_solid_class  (void);

/* GfkGlLevels: Header */

typedef struct _GfkGlLevels         GfkGlLevels;

struct _GfkGlLevels {
  /*< private >*/
  GfkGl2D parent;
};

#define GFK_GL_LEVELS(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlLevels,\
					         gfk_gl_levels_class ())
#define GFK_IS_GL_LEVELS(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_levels_class ()))

GfkGlClass * gfk_gl_levels_class  (void);

/* GfkGlVectors: Header */

typedef struct _GfkGlVectors         GfkGlVectors;

struct _GfkGlVectors {
  /*< private >*/
  GfkGlScalar parent;
  GtkWidget * vector;
};

#define GFK_GL_VECTORS(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlVectors,\
					         gfk_gl_vectors_class ())
#define GFK_IS_GL_VECTORS(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_vectors_class ()))

GfkGlClass * gfk_gl_vectors_class        (void);

/* GfkGlStreamlines: Header */

typedef struct _GfkGlStreamlines         GfkGlStreamlines;

struct _GfkGlStreamlines {
  /*< private >*/
  GfkGlVectors parent;
  GtkWidget * stream;
  gboolean edit;
};

#define GFK_GL_STREAMLINES(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlStreamlines,\
					         gfk_gl_streamlines_class ())
#define GFK_IS_GL_STREAMLINES(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_streamlines_class ()))

GfkGlClass * gfk_gl_streamlines_class        (void);

/* GfkGlEllipses: Header */

typedef struct _GfkGlEllipses         GfkGlEllipses;

struct _GfkGlEllipses {
  /*< private >*/
  GfkGlScalar parent;
  GtkWidget * ellipse;
};

#define GFK_GL_ELLIPSES(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlEllipses,\
					         gfk_gl_ellipses_class ())
#define GFK_IS_GL_ELLIPSES(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_ellipses_class ()))

GfkGlClass * gfk_gl_ellipses_class            (void);

/* GfkGlLocation: Header */

typedef struct _GfkGlLocation         GfkGlLocation;

struct _GfkGlLocation {
  /*< private >*/
  GfkGl parent;
  GtkWidget * location;
};

#define GFK_GL_LOCATION(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlLocation,\
					         gfk_gl_location_class ())

#define GFK_IS_GL_LOCATION(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_location_class ()))

GfkGlClass * gfk_gl_location_class  (void);

/* GfkGlHeight: Header */

#define GFK_IS_GL_HEIGHT(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_height_class ()))

GfkGlClass * gfk_gl_height_class  (void);

/* GfkGlLocate: Header */

typedef struct _GfkGlLocate         GfkGlLocate;

struct _GfkGlLocate {
  /*< private >*/
  GfkGl parent;
  GtkWidget * locate;
};

#define GFK_GL_LOCATE(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlLocate,\
					         gfk_gl_locate_class ())

#define GFK_IS_GL_LOCATE(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_locate_class ()))

GfkGlClass * gfk_gl_locate_class  (void);

/* GfkGlPipes: Header */

#define GFK_IS_GL_PIPES(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_pipes_class ()))

GfkGlClass * gfk_gl_pipes_class  (void);

/* GfkGlInfo: Header */

typedef struct _GfkGlInfo         GfkGlInfo;

struct _GfkGlInfo {
  /*< private >*/
  GfkGl parent;
  GtkWidget * info;
};

#define GFK_GL_INFO(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlInfo,\
					         gfk_gl_info_class ())
#define GFK_IS_GL_INFO(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_info_class ()))

GfkGlClass *  gfk_gl_info_class            (void);

/* GfkGlView: Header */

#define GFK_GL_PACK_MSG(msg, p)   (memcpy (&(msg).data.b[0], &(p), sizeof (GfsGl2PSParams *)))
#define GFK_GL_UNPACK_MSG(msg, p) (memcpy (&(p), &(msg)->data.b[0], sizeof (GfsGl2PSParams *)))

GtkWidget *     gfk_gl_view                (GtkWidget * glarea);
void            gfk_gl_view_set_scripting  (GtkWidget * view, 
					    gboolean active);
void            gfk_gl_view_set_simulation (GtkWidget * view, 
					    GfsSimulation * sim,
					    const gchar * fname);
gboolean        gfk_gl_view_read_parameters (GtkWidget * view, 
					     GtsFile * fp,
					     gboolean discard);
GfsSimulation * gfk_gl_simulation_read     (const gchar * fname,
					    GtkWidget * view,
					    gboolean set);
void            gfk_gl_view_draw           (GtkWidget * view,
					    guint format);
void            gfk_gl_view_pick           (GtkWidget * view, 
					    GfsGlRay * ray,
					    gboolean motion);
void            gfk_gl_view_clear          (GtkWidget * view);
void            gfs_gl2ps                  (GfsGl2PSParams * p, 
					    FILE * fp,
					    const gchar * fname, 
					    GtkWidget * view);

/* GfkGlVOF: Header */

typedef struct _GfkGlVOF         GfkGlVOF;

struct _GfkGlVOF {
  /*< private >*/
  GfkGlScalar parent;
  GtkWidget * scalar;

  /*< public >*/
};

#define GFK_GL_VOF(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlVOF,\
					         gfk_gl_vof_class ())
#define GFK_IS_GL_VOF(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_vof_class ()))

GfkGlClass * gfk_gl_vof_class            (void);

/* GfkGlClipPlane: Header */

#define GFK_IS_GL_CLIP_PLANE(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_clip_plane_class ()))

GfkGlClass * gfk_gl_clip_plane_class            (void);

#if (!FTT_2D)

/* GfkGlIsosurface: Header */

typedef struct _GfkGlIsosurface         GfkGlIsosurface;

struct _GfkGlIsosurface {
  /*< private >*/
  GfkGlScalar parent;
  GtkWidget * scalar;

  /*< public >*/
};

#define GFK_GL_ISOSURFACE(obj)            GTS_OBJECT_CAST (obj,\
					         GfkGlIsosurface,\
					         gfk_gl_isosurface_class ())
#define GFK_IS_GL_ISOSURFACE(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_isosurface_class ()))

GfkGlClass * gfk_gl_isosurface_class            (void);

/* GfkGlCutPlane: Header */

#define GFK_IS_GL_CUT_PLANE(obj)         (gts_object_is_from_class (obj,\
						 gfk_gl_cut_plane_class ()))

GfkGlClass * gfk_gl_cut_plane_class            (void);

#endif /* 3D */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __GFKGL_H__ */
