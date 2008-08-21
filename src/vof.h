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

#ifndef __VOF_H__
#define __VOF_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "advection.h"
#include "variable.h"

#define GFS_IS_FULL(f)             ((f) == 0. || (f) == 1.)

gdouble gfs_line_area              (FttVector * m, 
				    gdouble alpha);
void    gfs_line_center            (FttVector * m, 
				    gdouble alpha, 
				    gdouble a, 
				    FttVector * p);
gdouble gfs_line_area_center       (FttVector * m, 
				    gdouble alpha, 
				    FttVector * p);
gdouble gfs_line_alpha             (FttVector * m, 
				    gdouble c);
#if FTT_2D
#  define gfs_plane_volume         gfs_line_area
#  define gfs_plane_alpha          gfs_line_alpha
#  define gfs_plane_center         gfs_line_center
#  define gfs_plane_area_center     gfs_line_area_center
#else /* 3D */
gdouble gfs_plane_volume           (FttVector * m, 
				    gdouble alpha);
gdouble gfs_plane_alpha            (FttVector * m, 
				    gdouble c);
void    gfs_plane_center           (FttVector * m, 
				    gdouble alpha, 
				    gdouble a,
				    FttVector * p);
gdouble gfs_plane_area_center      (FttVector * m, 
				    gdouble alpha, 
				    FttVector * p);
#endif /* 3D */
void    gfs_youngs_gradient        (FttCell * cell, 
				    GfsVariable * v,
				    FttVector * g);

/* GfsVariableTracerVOF: header */

typedef struct _GfsVariableTracerVOF                GfsVariableTracerVOF;

struct _GfsVariableTracerVOF {
  /*< private >*/
  GfsVariableTracer parent;

  /*< public >*/
  GfsVariable * m[FTT_DIMENSION], * alpha;
};

#define GFS_VARIABLE_TRACER_VOF(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableTracerVOF,\
					           gfs_variable_tracer_vof_class ())
#define GFS_IS_VARIABLE_TRACER_VOF(obj)         (gts_object_is_from_class (obj,\
						   gfs_variable_tracer_vof_class ()))

GfsVariableClass * gfs_variable_tracer_vof_class  (void);

void     gfs_tracer_vof_advection  (GfsDomain * domain,
				    GfsAdvectionParams * par);
gdouble  gfs_vof_face_value        (const FttCellFace * face, 
				    GfsVariableTracerVOF * t);
guint    gfs_vof_facet             (FttCell * cell,
				    GfsVariableTracerVOF * t,
				    FttVector * p,
				    FttVector * m);
gdouble  gfs_vof_center            (FttCell * cell,
				    GfsVariableTracerVOF * t,
				    FttVector * p);
gdouble  gfs_vof_plane_interpolate (FttCell * cell,
				    FttVector * p,
				    guint level,
				    GfsVariableTracerVOF * t,
				    FttVector * m);
gdouble  gfs_vof_interpolate       (FttCell * cell,
				    FttVector * p,
				    guint level,
				    GfsVariableTracerVOF * t);
gdouble  gfs_height_curvature      (FttCell * cell, 
				    GfsVariableTracerVOF * t,
				    gdouble * kmax);
gdouble  gfs_fit_curvature         (FttCell * cell,
				    GfsVariableTracerVOF * t,
				    gdouble * kmax);
gdouble  gfs_vof_correctness       (FttCell * cell, 
				    GfsVariableTracerVOF * t);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __VOF_H__ */
