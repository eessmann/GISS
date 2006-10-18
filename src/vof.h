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

#define GFS_IS_FULL(f)             ((f) == 0. || (f) == 1.)

gdouble gfs_line_area              (FttVector * m, 
				    gdouble alpha);
void    gfs_line_center            (FttVector * m, 
				    gdouble alpha, 
				    gdouble a, 
				    FttVector * p);
gdouble gfs_line_alpha             (FttVector * m, 
				    gdouble c);
#if FTT_2D
#  define gfs_plane_volume         gfs_line_area
#  define gfs_plane_alpha          gfs_line_alpha
#else /* 3D */
gdouble gfs_plane_volume           (FttVector * m, 
				    gdouble alpha);
gdouble gfs_plane_alpha            (FttVector * m, 
				    gdouble c);
void    gfs_plane_center           (FttVector * m, 
				    gdouble alpha, 
				    gdouble a,
				    FttVector * p);
#endif /* 3D */
void    gfs_youngs_normal          (FttCell * cell, 
				    GfsVariable * v,
				    FttVector * n);
void    gfs_cell_vof_advection     (FttCell * cell,
				    FttComponent c,
				    GfsAdvectionParams * par);
void    gfs_tracer_vof_advection   (GfsDomain * domain,
				    GfsAdvectionParams * par);
void    gfs_vof_coarse_fine        (FttCell * parent, 
				    GfsVariable * v);
gboolean gfs_vof_plane             (FttCell * cell, 
				    GfsVariable * v,
				    FttVector * m, 
				    gdouble * alpha);
GSList * gfs_vof_facet             (FttCell * cell, 
				    GfsVariable * v);
gdouble  gfs_vof_interpolate       (FttCell * cell,
				    FttVector * p,
				    guint level,
				    GfsVariable * v);
gdouble  gfs_height_curvature      (FttCell * cell, 
				    GfsVariable * v);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __VOF_H__ */
