/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001 National Institute of Water and Atmospheric Research
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

#ifndef __FLUID_H__
#define __FLUID_H__

#include <glib.h>
#include <gts.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "ftt.h"

typedef struct _GfsStateVector     GfsStateVector;
typedef struct _GfsSolidVector     GfsSolidVector;
typedef struct _GfsFaceStateVector GfsFaceStateVector;

struct _GfsFaceStateVector {
  gdouble un;
  gdouble v;
};

struct _GfsStateVector {
  /* temporary face variables */
  GfsFaceStateVector f[FTT_NEIGHBORS];

  /* solid boundaries */
  GfsSolidVector * solid;

  /* centered temporary variables */
  gdouble div, dp, res;
  gdouble g[FTT_DIMENSION];

  /* centered primitive variables */
  gdouble p;
  gdouble u, v;
#if (!FTT_2D)
  gdouble w;
#endif /* FTT_3D */
};

typedef enum {
  /* centered temporary variables */
  GFS_DIV = 0,
  GFS_DP,
  GFS_RES,
  GFS_GX,
  GFS_GY,
#if (!FTT_2D)
  GFS_GZ,
#endif /* FTT_3D */
  /* centered primitive variables */
  GFS_P,
  GFS_U, GFS_V,
#if (!FTT_2D)
  GFS_W,
#endif /* FTT_3D */
} GfsPermanentVariable;

struct _GfsSolidVector {
  gdouble s[FTT_NEIGHBORS];
  gdouble a, v, fv;
  FttCell * merged;
  FttVector cm, ca;
};

typedef enum {
  GFS_FLAG_USED =      1 <<  FTT_FLAG_USER,
  GFS_FLAG_BOUNDARY  = 1 << (FTT_FLAG_USER + 1),
  GFS_FLAG_DIRICHLET = 1 << (FTT_FLAG_USER + 2),
  GFS_FLAG_USER =            FTT_FLAG_USER + 3 /* user flags start here */
} GfsFlags;

/* GfsVariable: Header */

typedef struct _GfsVariable                GfsVariable;
typedef struct _GfsSurfaceGenericBc        GfsSurfaceGenericBc;

typedef void (* GfsVariableDerivedFunc)    (FttCell * cell, GfsVariable * v);
typedef void (* GfsVariableFineCoarseFunc) (FttCell * cell, GfsVariable * v);

struct _GfsVariable {
  /*< private >*/
  GtsObject parent;

  /*< public >*/
  guint i;
  gchar * name;
  gboolean centered;
  GfsVariableDerivedFunc derived;
  GfsVariableFineCoarseFunc fine_coarse;
  GtsContainer * sources;
  GfsSurfaceGenericBc * surface_bc;
  GfsVariable * next, * permanent;
  GtsObject * p;
};

typedef struct _GfsVariableClass    GfsVariableClass;

struct _GfsVariableClass {
  /*< private >*/
  GtsObjectClass parent_class;

  /*< public >*/
};

#define GFS_VARIABLE1(obj)            GTS_OBJECT_CAST (obj,\
					         GfsVariable,\
					         gfs_variable_class ())
#define GFS_VARIABLE_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsVariableClass,\
						 gfs_variable_class())
#define GFS_IS_VARIABLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_variable_class ()))
#define gfs_variable_parent(v)         ((v)->p)
#define gfs_variable_set_parent(v, pa) ((v)->p = GTS_OBJECT (pa))

GfsVariableClass * gfs_variable_class  (void);
GfsVariable *      gfs_variable_new    (GfsVariableClass * klass,
					GtsObject * parent,
					const gchar * name,
					gboolean centered,
					guint i);

/* Permanent variables: defined in fluid.c */
GTS_C_VAR GfsVariable * gfs_div, * gfs_dp, * gfs_res;
GTS_C_VAR GfsVariable * gfs_gx, * gfs_gy;
#if (!FTT_2D)
GTS_C_VAR GfsVariable * gfs_gz;
#endif /* FTT_3D */
GTS_C_VAR GfsVariable * gfs_centered_variables;
GTS_C_VAR GfsVariable * gfs_p;

/* Derived variables: defined in fluid.c */
GTS_C_VAR GfsVariable * gfs_derived_first, * gfs_derived_last;

#define GFS_STATE(cell)               ((GfsStateVector *) (cell)->data)
#define GFS_VARIABLE(cell, index)     ((&GFS_STATE (cell)->div)[index])
#define GFS_VELOCITY_INDEX(component) (GFS_U + (component))
#define GFS_GRADIENT_INDEX(component) (GFS_GX + (component))
#define GFS_VELOCITY_COMPONENT(index) ((index) - GFS_U)
#define GFS_GRADIENT_COMPONENT(index) ((index) - GFS_GX)

#define GFS_FACE_NORMAL_VELOCITY(fa)\
  (GFS_STATE ((fa)->cell)->f[(fa)->d].un)
#define GFS_FACE_NORMAL_VELOCITY_LEFT(fa)\
  (GFS_STATE ((fa)->cell)->f[(fa)->d].un)
#define GFS_FACE_NORMAL_VELOCITY_RIGHT(fa)\
  (GFS_STATE ((fa)->neighbor)->f[FTT_OPPOSITE_DIRECTION ((fa)->d)].un)

#define GFS_FACE_FRACTION(fa) (GFS_IS_MIXED ((fa)->cell) ?\
                               GFS_STATE ((fa)->cell)->solid->s[(fa)->d] : 1.)
#define GFS_FACE_FRACTION_LEFT(fa) GFS_FACE_FRACTION(fa)
#define GFS_FACE_FRACTION_RIGHT(fa) (GFS_IS_MIXED ((fa)->neighbor) ?\
                 GFS_STATE ((fa)->neighbor)->solid->s[FTT_OPPOSITE_DIRECTION ((fa)->d)] : 1.)

#define GFS_IS_FLUID(cell)      ((cell) != NULL &&\
                                 GFS_STATE (cell)->solid == NULL)
#define GFS_IS_MIXED(cell)      ((cell) != NULL &&\
                                 GFS_STATE (cell)->solid != NULL)
#define GFS_IS_SMALL(cell)      (GFS_IS_MIXED (cell) &&\
                                 GFS_STATE (cell)->solid->a < 0.5)
#define GFS_CELL_IS_BOUNDARY(cell) (((cell)->flags & GFS_FLAG_BOUNDARY) != 0)

void                  gfs_cell_cleanup              (FttCell * cell);
void                  gfs_cell_reset                (FttCell * cell, 
						     GfsVariable * v);
void                  gfs_get_from_above            (FttCell * cell, 
						     const GfsVariable * v);
void                  gfs_get_from_below_intensive  (FttCell * cell, 
						     const GfsVariable * v);
void                  gfs_get_from_below_extensive  (FttCell * cell, 
						     const GfsVariable * v);
gdouble               gfs_face_interpolated_value   (const FttCellFace * face,
						     guint v);
typedef gdouble    (* GfsCenterGradient)            (FttCell * cell,
						     FttComponent c,
						     guint v);
gdouble               gfs_center_gradient           (FttCell * cell,
						     FttComponent c,
						     guint v);
gdouble               gfs_center_van_leer_gradient  (FttCell * cell,
						     FttComponent c,
						     guint v);

typedef struct _GfsGradient GfsGradient;

struct _GfsGradient {
  gdouble a, b;
};

void                  gfs_face_gradient              (const FttCellFace * face,
						      GfsGradient * g,
						      guint v,
						      gint max_level);
void                  gfs_face_weighted_gradient     (const FttCellFace * face,
						      GfsGradient * g,
						      guint v,
						      gint max_level);
void                  gfs_face_gradient_flux         (const FttCellFace * face,
						      GfsGradient * g,
						      guint v,
						      gint max_level);
void                  gfs_cell_dirichlet_gradient    (FttCell * cell,
						      guint v,
						      gint max_level,
						      gdouble v0,
						      FttVector * grad);
gdouble               gfs_cell_dirichlet_gradient_flux (FttCell * cell,
							guint v,
							gint max_level,
							gdouble v0);
gdouble               gfs_cell_dirichlet_value         (FttCell * cell,
							GfsVariable * v,
							gint max_level);
void                  gfs_face_gradient_flux_centered(const FttCellFace * face,
						      GfsGradient * g,
						      guint v,
						      gint max_level);

void                  gfs_normal_divergence          (FttCell * cell);
void                  gfs_normal_divergence_2D       (FttCell * cell);
void                  gfs_divergence                 (FttCell * cell);
gdouble               gfs_vorticity_value            (FttCell * cell,
						      FttVector * lambda);
void                  gfs_vorticity                  (FttCell * cell,
						      GfsVariable * v);
void                  gfs_velocity_norm              (FttCell * cell,
						      GfsVariable * v);
void                  gfs_velocity_norm2             (FttCell * cell,
						      GfsVariable * v);
void                  gfs_velocity_lambda2           (FttCell * cell,
						      GfsVariable * v);
GtsRange              gfs_stats_variable             (FttCell * root, 
						      GfsVariable * v, 
						      FttTraverseFlags flags,
						      gint max_depth);

typedef struct _GfsNorm GfsNorm;

struct _GfsNorm {
  gdouble bias, first, second, infty, w;
};

void                  gfs_norm_init                 (GfsNorm * n);
void                  gfs_norm_reset                (GfsNorm * n);
void                  gfs_norm_add                  (GfsNorm * n, 
						     gdouble val,
						     gdouble weight);
void                  gfs_norm_update               (GfsNorm * n);

GfsNorm               gfs_norm_variable             (FttCell * root, 
						     GfsVariable * v, 
						     FttTraverseFlags flags,
						     gint max_depth);
  
void                  gfs_cell_traverse_mixed       (FttCell * root,
						     FttTraverseType order,
						     FttTraverseFlags flags,
						     FttCellTraverseFunc func,
						     gpointer data);
typedef void       (* FttCellTraverseCutFunc)       (FttCell * cell,
						     GtsSurface * s,
						     gpointer data);
void                  gfs_cell_traverse_cut         (FttCell * root,
						     GtsSurface * s,
						     FttTraverseType order,
						     FttTraverseFlags flags,
						     FttCellTraverseCutFunc func,
						     gpointer data);
gdouble               gfs_interpolate               (FttCell * cell,
						     FttVector p,
						     GfsVariable * v);
GfsVariable *         gfs_variable_list_copy        (GfsVariable * v,
						     GtsObject * parent);
void                  gfs_variable_list_destroy     (GfsVariable * v);
GfsVariable *         gfs_variable_from_name        (GfsVariable * variables,
						     const gchar * name);
GfsVariable *         gfs_variables_from_list       (GfsVariable * variables,
						     gchar * list,
						     gchar ** error);

void                  ftt_cell_refine_corners       (FttCell * cell,
						     FttCellInitFunc init,
						     gpointer data);
gdouble               gfs_center_curvature          (FttCell * cell,
						     FttComponent c,
						     guint v);
gdouble               gfs_streamline_curvature      (FttCell * cell);
gdouble               gfs_cell_laplacian            (FttCell * cell, 
						     GfsVariable * v);

typedef struct {
#if FTT_2D
  FttCell * c[7];
  gdouble w[7];
#else  /* 3D */
  FttCell * c[29];
  gdouble w[29];
#endif /* 3D */
  guint n;  
} GfsInterpolator;

void                  gfs_cell_corner_interpolator  (FttCell * cell,
						     FttDirection d[FTT_DIMENSION],
						     gint max_level,
						     gboolean centered,
						     GfsInterpolator * inter);
gdouble               gfs_cell_corner_value         (FttCell * cell,
						     FttDirection * d,
						     GfsVariable * v,
						     gint max_level);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FLUID_H__ */
