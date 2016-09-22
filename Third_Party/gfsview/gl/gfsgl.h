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

#ifndef __GFSGL_H__
#define __GFSGL_H__

#include <gfs.h>

#include "gl2ps/gl2ps.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

FILE * gfs_gl_popen (const gchar * fname);

#define GFS_COLORMAP_TEXTURE_SAMPLES 256

typedef struct _GfsColormap GfsColormap;

struct _GfsColormap {
  GPtrArray * colors;
  gboolean reversed;
  GLfloat texture[3*GFS_COLORMAP_TEXTURE_SAMPLES];
  gchar * name;
};

GfsColormap *     gfs_colormap_jet                 (void);
GfsColormap *     gfs_colormap_cool_warm           (void);
GfsColormap *     gfs_colormap_gray                (void);
GtsColor          gfs_colormap_color               (GfsColormap * cmap, 
						    gdouble val);
void              gfs_colormap_texture             (GfsColormap * cmap);
void              gfs_colormap_destroy             (GfsColormap * colormap);

typedef struct _GfsGlViewParams GfsGlViewParams;

typedef struct {
  float m[16], p[16], res;
  GtsVector n[6];
  gdouble d[6];
  guint width;

  GList * symmetries;
  FttVector * s;
} GfsFrustum;

void          gfs_gl_get_frustum           (GfsGlViewParams * vp,
					    GList * symmetries,
					    GfsFrustum * f);
void          gfs_gl_frustum_free          (GfsFrustum * f);
GtsIntersect  gfs_sphere_in_frustum        (FttVector * p, 
					    gdouble r, 
					    GfsFrustum * f);
gboolean      gfs_sphere_is_small          (FttVector * c, 
					    gdouble r, 
					    GfsFrustum * f);

void          gfs_gl_init                  (void);
void          gfs_gl_init_gl               (void);

/* GfsGl: Header */

typedef struct _GfsGl           GfsGl;
typedef struct _GfsGl2D         GfsGl2D;

typedef enum {
  GFS_GL_CONSTANT,
  GFS_GL_FLAT,
  GFS_GL_SMOOTH,
  GFS_GL_CSMOOTH
} GfsGlShading;

struct _GfsGl {
  /*< private >*/
  GtsObject parent;

  /*< public >*/
  GfsSimulation * sim;
  GfsGlViewParams * p;
  guint size, format;
  
  GtsColor lc;
  GfsGlShading shading;
  gint maxlevel;

  gfloat font_size;
  gboolean use_raster_font;
  gfloat line_width;
};

typedef struct _GfsGlClass    GfsGlClass;
typedef struct {
  FttVector a, b;
} GfsGlRay;

struct _GfsGlClass {
  /*< private >*/
  GtsObjectClass parent_class;

  /*< public >*/
  void          (* set_simulation)  (GfsGl *, GfsSimulation *);
  void          (* draw)            (GfsGl *, GfsFrustum *);
  void          (* cut)             (GfsGl *, FttCell *, GfsGl2D *);
  gdouble       (* pick)            (GfsGl *, GfsGlRay *);
  gboolean      (* relevant)        (GfsSimulation *);
};

#define GFS_GL(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGl,\
					         gfs_gl_class ())
#define GFS_GL_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsGlClass,\
						 gfs_gl_class())
#define GFS_IS_GL(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_class ()))

GfsGlClass *  gfs_gl_class                 (void);
GfsGl *       gfs_gl_new                   (GfsGlClass * klass);
GfsGl *       gfs_gl_new_from_file         (GtsFile * fp);
void          gfs_gl_cell_traverse_visible (GfsGl * gl,
					    GfsFrustum * f,
					    FttCellTraverseFunc func,
					    gpointer data);
void          gfs_gl_cell_traverse_visible_condition (GfsGl * gl,
						      GfsFrustum * f,
						      gboolean (* condition) (FttCell *, gpointer),
						      gpointer datacon,
						      FttCellTraverseFunc func,
						      gpointer data);
void          gfs_gl_cell_traverse_visible_mixed (GfsGl * gl,
						  GfsFrustum * f,
						  FttCellTraverseFunc func,
						  gpointer data);
void          gfs_gl_cell_traverse_visible_boundary (GfsGl * gl,
						     GfsFrustum * f,
						     FttCellTraverseFunc func,
						     gpointer data);
void          gfs_gl_set_simulation        (GfsGl * gl, 
					    GfsSimulation * sim);
void          gfs_gl_draw                  (GfsGl * gl, 
					    GfsFrustum * f);
void          gfs_gl_set_raster_font       (GfsGl * gl, 
					    gboolean raster);
void          gfs_gl_set_font_size         (GfsGl * gl, 
					    gfloat size);

#if FTT_2D
# define      gfs_gl_cell_traverse_visible_plane gfs_gl_cell_traverse_visible
#endif /* 2D */

#define       gfs_gl_normal(gl)            (glNormal3d (GFS_GL2D (gl)->n.x,\
                                            GFS_GL2D (gl)->n.y,\
                                            GFS_GL2D (gl)->n.z))
#define       gfs_gl_vector_format(gl)     (gl->format != GFSGL_SCREEN && \
                                            gl->format != GFSGL_PPM_OFFSCREEN)

/* GfsGlLabel: Header */

typedef struct _GfsGlLabel GfsGlLabel;

struct _GfsGlLabel {
  /*< private >*/
  GfsGl parent;
  gchar * formatted_label;

  /*< public >*/
  FttVector p;
  gboolean symbol;
  gchar * label;
};

#define GFS_GL_LABEL(obj)            GTS_OBJECT_CAST (obj,\
							 GfsGlLabel,	\
							 gfs_gl_label_class ())
#define GFS_IS_GL_LABEL(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_label_class ()))

GfsGlClass *   gfs_gl_label_class     (void);
void           gfs_gl_label_set_label (GfsGlLabel * gl, 
				       const gchar * label, 
				       GfsSimulation * sim);

/* GfsGl2D: Header */

struct _GfsGl2D {
  /*< private >*/
  GfsGl parent;
  FttVector n;
  gdouble pos;
  FttCell * picked;
  FttVector pickedpos;

  /*< public >*/
  FttVector p[3];
};

typedef struct _GfsGl2DClass    GfsGl2DClass;

struct _GfsGl2DClass {
  /*< private >*/
  GfsGlClass parent_class;

  /*< public >*/
  void (* update_plane)  (GfsGl2D *);
};

#define GFS_GL2D(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGl2D,\
					         gfs_gl2D_class ())
#define GFS_GL2D_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsGl2DClass,\
						 gfs_gl2D_class())
#define GFS_IS_GL2D(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl2D_class ()))

GfsGl2DClass * gfs_gl2D_class             (void);
void           gfs_gl2D_update_plane      (GfsGl2D * gl);
void           gfs_gl2D_update_pos_bounds (GfsGl2D * gl);

/* GfsGlSymmetry: Header */

typedef struct _GfsGlSymmetry GfsGlSymmetry;

struct _GfsGlSymmetry {
  /*< private >*/
  GfsGl2D parent;

  /*< public >*/
  GLfloat m[16];
  gboolean periodic;
};

#define GFS_GL_SYMMETRY(obj)            GTS_OBJECT_CAST (obj,\
							 GfsGlSymmetry,	\
							 gfs_gl_symmetry_class ())
#define GFS_IS_GL_SYMMETRY(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_symmetry_class ()))

GfsGlClass *   gfs_gl_symmetry_class  (void);
void           gfs_gl_symmetry_apply  (GList * symmetry,
				       GLuint display_list);
void           gfs_gl_symmetry_transform (GfsGl * gl,
					  FttVector * p,
					  FttVector * t);

/* GfsGlPeriodic: Header */

GfsGlClass *   gfs_gl_periodic_class  (void);

/* GfsGlCells: Header */

#define GFS_IS_GL_CELLS(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_cells_class ()))

GfsGlClass * gfs_gl_cells_class  (void);

/* GfsGlFractions: Header */

#define GFS_IS_GL_FRACTIONS(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_fractions_class ()))

GfsGlClass * gfs_gl_fractions_class  (void);

/* GfsGlBoundaries: Header */

#define GFS_IS_GL_BOUNDARIES(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_boundaries_class ()))

GfsGlClass * gfs_gl_boundaries_class  (void);

/* GfsGlVarFunc: Header */

typedef struct {
  GfsVariable * v;
  GfsFunction * f;
} GfsGlVarFunc;

GfsGlVarFunc *   gfs_gl_var_func_new     (void);
void             gfs_gl_var_func_destroy (GfsGlVarFunc * vf);
GtsFile *        gfs_gl_var_func_set     (GfsGlVarFunc * vf, 
					  GfsSimulation * sim, 
					  const gchar * func,
					  GString * expr,
					  GfsVariableClass * klass);

/* GfsGlScalar: Header */

typedef struct _GfsGlScalar         GfsGlScalar;

struct _GfsGlScalar {
  /*< private >*/
  GfsGl2D parent;
  GfsGlVarFunc * vf;
  gboolean amin, amax;
  gdouble aminv, amaxv;

  /*< public >*/
  GString * expr;
  GfsVariable * v;
  gdouble min, max;
  GfsColormap * cmap;
  gboolean show;
};

typedef struct _GfsGlScalarClass    GfsGlScalarClass;

struct _GfsGlScalarClass {
  /*< private >*/
  GfsGl2DClass parent_class;

  /*< public >*/
  GtsFile * (* set_scalar)  (GfsGlScalar *, const gchar * func);
};

#define GFS_GL_SCALAR(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlScalar,\
					         gfs_gl_scalar_class ())
#define GFS_GL_SCALAR_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsGlScalarClass,\
						 gfs_gl_scalar_class())
#define GFS_IS_GL_SCALAR(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_scalar_class ()))

GfsGlScalarClass *  gfs_gl_scalar_class        (void);
GtsFile *           gfs_gl_scalar_set          (GfsGlScalar * gl,
						const gchar * func);

/* GfsGlSquares: Header */

typedef struct _GfsGlSquares         GfsGlSquares;

struct _GfsGlSquares {
  /*< private >*/
  GfsGlScalar parent;

  /*< public >*/
};

#define GFS_GL_SQUARES(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlSquares,\
					         gfs_gl_squares_class ())
#define GFS_IS_GL_SQUARES(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_squares_class ()))

GfsGlClass * gfs_gl_squares_class  (void);

/* GfsGlLinear: Header */

typedef struct _GfsGlLinear         GfsGlLinear;

struct _GfsGlLinear {
  /*< private >*/
  GfsGlScalar parent;
  GfsGlVarFunc * vf;
  GfsVariable * nx, * ny;
  
  /*< public >*/
  GString * expr;
  GfsVariable * use_scalar;
  gboolean reversed;
};

#define GFS_GL_LINEAR(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlLinear,\
					         gfs_gl_linear_class ())
#define GFS_IS_GL_LINEAR(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_linear_class ()))

GfsGlClass * gfs_gl_linear_class  (void);
GtsFile *    gfs_gl_linear_set    (GfsGlLinear * gl, 
				   const gchar * func);

/* GfsGlIsoline: Header */

typedef struct _GfsGlIsoline         GfsGlIsoline;

struct _GfsGlIsoline {
  /*< private >*/
  GfsGlLinear parent;
  GArray * levels;
  GfsVariable * used;
  GfsFrustum * f;
  gdouble val;

  /*< public >*/
  GfsVariable * min, * max;
  gchar * ls;
  gdouble n;
};

#define GFS_GL_ISOLINE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlIsoline,\
					         gfs_gl_isoline_class ())
#define GFS_IS_GL_ISOLINE(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_isoline_class ()))

GfsGlClass * gfs_gl_isoline_class             (void);
void         gfs_gl_isoline_set_levels        (GfsGlIsoline * gl, 
					       const gchar * levels);
void         gfs_gl_cell_traverse_visible_iso (GfsGl * gl,
					       GfsFrustum * f,
					       GfsVariable * min,
					       GfsVariable * max,
					       gdouble level,
					       FttCellTraverseFunc func,
					       gpointer data);

/* GfsGlVOF: Header */

typedef struct _GfsGlVOF         GfsGlVOF;

struct _GfsGlVOF {
  /*< private >*/
  GfsGlScalar parent;
  GfsGlVarFunc * vf;

  /*< public >*/
  GString * expr;
  GfsVariable * use_scalar;
  gboolean reversed, draw_edges, interpolate;
};

#define GFS_GL_VOF(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlVOF,\
					         gfs_gl_vof_class ())
#define GFS_IS_GL_VOF(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_vof_class ()))

GfsGlClass * gfs_gl_vof_class      (void);
GtsFile *    gfs_gl_vof_set        (GfsGlVOF * gl, 
				    const gchar * func);

/* GfsGlSolid: Header */

#if (!FTT_2D)
typedef struct _GfsGlSolid         GfsGlSolid;

struct _GfsGlSolid {
  /*< private >*/
  GfsGlScalar parent;
  gboolean needs_updating;

  /*< public >*/
  gboolean reversed;
  GfsVariable * p, * s, * use_scalar;
  GSList * solids;
};

#define GFS_GL_SOLID(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlSolid,\
					         gfs_gl_solid_class ())
#endif /* 3D */

#define GFS_IS_GL_SOLID(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_solid_class ()))

GfsGlClass * gfs_gl_solid_class  (void);
#if (!FTT_2D)
void         gfs_gl_solid_reset  (GfsGlSolid * gl);
#endif /* 3D */

/* GfsGlLevels: Header */

typedef struct _GfsGlLevels         GfsGlLevels;

struct _GfsGlLevels {
  /*< private >*/
  GfsGl2D parent;

  /*< public >*/
  GfsVariable * v;
};

#define GFS_GL_LEVELS(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlLevels,\
					         gfs_gl_levels_class ())
#define GFS_IS_GL_LEVELS(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_levels_class ()))

GfsGlClass * gfs_gl_levels_class  (void);

/* GfsGlVectors: Header */

typedef struct _GfsGlVectors         GfsGlVectors;

struct _GfsGlVectors {
  /*< private >*/
  GfsGlScalar parent;
  gdouble h, max;
  gboolean already_set;
  GfsGlVarFunc * vf[FTT_DIMENSION];

  /*< public >*/
  GString * expr[FTT_DIMENSION];
  GfsVariable * v[FTT_DIMENSION];
  gdouble scale;
  GtsColor c;
  gboolean use_scalar;
};

#define GFS_GL_VECTORS(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlVectors,\
					         gfs_gl_vectors_class ())
#define GFS_IS_GL_VECTORS(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_vectors_class ()))

GfsGlClass * gfs_gl_vectors_class        (void);
GtsFile *    gfs_gl_vectors_set          (GfsGlVectors * gl,
					  FttComponent c,
					  const gchar * func);

/* GfsGlStreamline: Header */

typedef struct _GfsGlStreamline         GfsGlStreamline;

struct _GfsGlStreamline {
  /*< private >*/
  GtsObject parent;
  GLuint list;

  /*< public >*/
  FttVector c;
  GList * l;
};

#define GFS_GL_STREAMLINE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlStreamline,\
					         gfs_gl_streamline_class ())
#define GFS_IS_GL_STREAMLINE(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_streamline_class ()))

GtsObjectClass * gfs_gl_streamline_class        (void);

/* GfsGlStreamlines: Header */

typedef struct _GfsGlStreamlines         GfsGlStreamlines;

struct _GfsGlStreamlines {
  /*< private >*/
  GfsGlVectors parent;
  GfsVariable * s;
  GList * stream, * selected, * candidate;

  /*< public >*/
  gboolean show_cells;
  gdouble dmin, radius;
};

#define GFS_GL_STREAMLINES(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlStreamlines,\
					         gfs_gl_streamlines_class ())
#define GFS_IS_GL_STREAMLINES(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_streamlines_class ()))

GfsGlClass * gfs_gl_streamlines_class           (void);
void         gfs_gl_streamlines_reset           (GfsGlStreamlines * gl);
void         gfs_gl_streamlines_reset_selected  (GfsGlStreamlines * gl);
GfsGlStreamline * gfs_gl_streamlines_add        (GfsGlStreamlines * gl, 
						 FttVector p);
gboolean     gfs_gl_streamlines_remove_selected (GfsGlStreamlines * gl);
void         gfs_gl_streamlines_update_display_lists (GfsGlStreamlines * gl);
gdouble      gfs_gl_streamlines_closest         (GfsGlStreamlines * gl, 
						 FttVector * p,
						 GtsPoint * closest);
void         gfs_gl_streamlines_evenly_spaced   (GfsGlStreamlines * gl,
						 gboolean (* callback) (GfsGlStreamlines *, 
									gpointer),
						 gpointer data);

/* GfsGlEllipses: Header */

typedef struct _GfsGlEllipses         GfsGlEllipses;

struct _GfsGlEllipses {
  /*< private >*/
  GfsGlScalar parent;
  gdouble h, max;
  gboolean already_set;
  GfsGlVarFunc * vf[4];

  /*< public >*/
  GString * expr[4];
  GfsVariable * v[4];
  gdouble scale;
  GtsColor c;
  gboolean use_scalar;
};

#define GFS_GL_ELLIPSES(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlEllipses,\
					         gfs_gl_ellipses_class ())
#define GFS_IS_GL_ELLIPSES(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_ellipses_class ()))

GfsGlClass * gfs_gl_ellipses_class        (void);
GtsFile *    gfs_gl_ellipses_set          (GfsGlEllipses * gl, 
					   guint i, 
					   const gchar * func);

/* GfsGlLocation: Header */

typedef struct _GfsGlLocation         GfsGlLocation;

struct _GfsGlLocation {
  /*< private >*/
  GfsGl parent;

  /*< public >*/
  gdouble size;
  gboolean label;
};

#define GFS_GL_LOCATION(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlLocation,\
					         gfs_gl_location_class ())

#define GFS_IS_GL_LOCATION(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_location_class ()))

GfsGlClass * gfs_gl_location_class  (void);

/* GfsGlHeight: Header */

#define GFS_IS_GL_HEIGHT(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_height_class ()))

GfsGlClass * gfs_gl_height_class  (void);

/* GfsGlLocate: Header */

typedef struct _GfsGlLocate         GfsGlLocate;

struct _GfsGlLocate {
  /*< private >*/
  GfsGl parent;

  /*< public >*/
  FttVector p;
};

#define GFS_GL_LOCATE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlLocate,\
					         gfs_gl_locate_class ())
#define GFS_IS_GL_LOCATE(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_locate_class ()))

GfsGlClass *  gfs_gl_locate_class             (void);

/* GfsGlPipes: Header */

#define GFS_IS_GL_PIPES(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_pipes_class ()))

GfsGlClass * gfs_gl_pipes_class  (void);

struct _GfsGlViewParams {
  gboolean do_init;
  gfloat beginx, beginy;
  gfloat dx, dy, tx, ty, sx, sy, sz;
  gfloat quat[4];
  gfloat dquat[4];
  gfloat fov;

  GtsColor bg;
  gfloat lc;
  gfloat timing, base_res, res, reactivity;
  gboolean motion, cp[6];
  gfloat lw; /* base line width for vector drawings */
};

void              gfs_gl_view_params_init    (GfsGlViewParams * p);
GfsGlViewParams * gfs_gl_view_params_new     (void);
void              gfs_gl_view_params_write   (GfsGlViewParams * p, 
					      FILE * fp);
void              gfs_gl_view_params_read    (GfsGlViewParams * p, 
					      GtsFile * fp);

typedef struct {
  GLint format, sort, options;
  guint width, height;
  FILE * fp;
  gfloat lw; /* base line width */
} GfsGl2PSParams;

typedef enum {
  GFSGL_PPM_OFFSCREEN = GL2PS_PGF + 1,
  GFSGL_PPM_SCREEN,
  GFSGL_SCREEN,
  GFSGL_GFSVIEW,
  GFSGL_GERRIS,
  GFSGL_GNUPLOT,
  GFSGL_OBJ,
  GFSGL_KML,
  GFSGL_PPM = GFSGL_PPM_OFFSCREEN
} GfsGlFormat;

void              gfs_gl2ps_params_init      (GfsGl2PSParams * p);
void              gfs_gl2ps_params_read      (GfsGl2PSParams * p,
					      GtsFile * fp);
void              gfs_gl2ps_params_write     (GfsGl2PSParams * p,
					      FILE * fp);
void              gfs_gl_write_image         (FILE * fp, 
					      const GLubyte * buffer, 
					      guint width, 
					      guint height);

/* GfsGlClipPlane: Header */

typedef struct _GfsGlClipPlane         GfsGlClipPlane;

struct _GfsGlClipPlane {
  /*< private >*/
  GfsGl2D parent;
  gint i;

  /*< public >*/
  gboolean disabled;
};

#define GFS_GL_CLIP_PLANE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlClipPlane,\
					         gfs_gl_clip_plane_class ())
#define GFS_IS_GL_CLIP_PLANE(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_clip_plane_class ()))

GfsGlClass * gfs_gl_clip_plane_class        (void);
void         gfs_gl_clip_plane_disable      (GfsGlClipPlane * gl);

/* GfsGlCutPlane: Header */

typedef struct _GfsGlCutPlane         GfsGlCutPlane;

struct _GfsGlCutPlane {
  /*< private >*/
  GfsGl2D parent;
  GList * list;

  /*< public >*/
};

#define GFS_GL_CUT_PLANE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlCutPlane,\
					         gfs_gl_cut_plane_class ())
#define GFS_IS_GL_CUT_PLANE(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_cut_plane_class ()))

GfsGlClass * gfs_gl_cut_plane_class        (void);

#if (!FTT_2D)

/* GfsGlIsosurface: Header */

typedef struct _GfsGlIsosurface         GfsGlIsosurface;

struct _GfsGlIsosurface {
  /*< private >*/
  GfsGlScalar parent;
  GfsGlVarFunc * vf;

  /*< public >*/
  GString * expr;
  GfsVariable * v, * min, * max, * p, * use_scalar;
  gdouble level, minv, maxv;
  gboolean reversed;
};

#define GFS_GL_ISOSURFACE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlIsosurface,\
					         gfs_gl_isosurface_class ())
#define GFS_IS_GL_ISOSURFACE(obj)         (gts_object_is_from_class (obj,\
						 gfs_gl_isosurface_class ()))

GfsGlClass * gfs_gl_isosurface_class        (void);
void         gfs_gl_isosurface_reset        (GfsGlIsosurface * gl);
GtsFile *    gfs_gl_isosurface_set          (GfsGlIsosurface * gl, 
					     const gchar * func);

#endif /* 3D */

gdouble      gfs_gl_domain_extent           (GfsDomain * domain,
					     GList * symmetries);

typedef struct {
  GLfloat * feedback;
} GfsGlFeedback;

GfsGlFeedback * gfs_gl_feedback_begin (guint buffersize);
gboolean        gfs_gl_feedback_end   (GfsGlFeedback * f,
				       GfsSimulation * sim,
				       FILE * fp,
				       GfsGlFormat format);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __GFSGL_H__ */
