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

#ifndef __TIMESTEP_H__
#define __TIMESTEP_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "advection.h"
#include "poisson.h"

typedef struct _GfsMultilevelParams GfsMultilevelParams;

struct _GfsMultilevelParams {
  gdouble tolerance;
  guint nrelax;
  guint minlevel;
  guint nitermax;

  guint dimension;
  guint niter;
  GfsNorm residual_before, residual;
};

#include "variable.h"

void          gfs_multilevel_params_init      (GfsMultilevelParams * par);
void          gfs_multilevel_params_write     (GfsMultilevelParams * par, 
					       FILE * fp);
void          gfs_multilevel_params_read      (GfsMultilevelParams * par, 
					       GtsFile * fp);
void          gfs_correct_normal_velocities   (GfsDomain * domain,
					       guint dimension,
					       GfsVariable * p,
					       GfsVariable ** g,
					       gdouble dt);
void          gfs_mac_projection              (GfsDomain * domain,
					       GfsMultilevelParams * par,
					       GfsAdvectionParams * apar,
					       GfsVariable * p,
					       GfsVariable ** g);
void          gfs_correct_centered_velocities (GfsDomain * domain,
					       guint dimension,
					       GfsVariable ** g,
					       gdouble dt);
void          gfs_approximate_projection      (GfsDomain * domain,
					       GfsMultilevelParams * par,
					       GfsAdvectionParams * apar,
					       GfsVariable * p,
					       GfsVariable * res);
void          gfs_predicted_face_velocities   (GfsDomain * domain,
					       guint d,
					       GfsAdvectionParams * par);

void          gfs_diffusion                   (GfsDomain * domain,
					       GfsMultilevelParams * par,
					       GfsVariable * v,
					       GfsVariable * rhs, 
					       GfsVariable * dia);
void          gfs_centered_velocity_advection_diffusion (GfsDomain * domain,
							 guint dimension,
							 GfsAdvectionParams * par,
							 GfsMultilevelParams * dpar,
							 GfsVariable ** g);
void          gfs_tracer_advection_diffusion  (GfsDomain * domain,
					       GfsAdvectionParams * par,
					       GfsMultilevelParams * dpar,
					       GfsVariable * half);

/* GfsSurfaceGenericBc: Header */

struct _GfsSurfaceGenericBc {
  /*< private >*/
  GtsObject parent;

  /*< public >*/
  GfsVariable * v;  
};

typedef struct _GfsSurfaceGenericBcClass    GfsSurfaceGenericBcClass;

struct _GfsSurfaceGenericBcClass {
  /*< private >*/
  GtsObjectClass parent_class;

  /*< public >*/
  void (* bc) (FttCell *, GfsSurfaceGenericBc *);
};

#define GFS_SURFACE_GENERIC_BC(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSurfaceGenericBc,\
					         gfs_surface_generic_bc_class ())
#define GFS_SURFACE_GENERIC_BC_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsSurfaceGenericBcClass,\
						 gfs_surface_generic_bc_class())
#define GFS_IS_SURFACE_GENERIC_BC(obj)         (gts_object_is_from_class (obj,\
						 gfs_surface_generic_bc_class ()))

GfsSurfaceGenericBcClass * gfs_surface_generic_bc_class  (void);

/* GfsSurfaceBc: Header */

typedef struct _GfsSurfaceBc         GfsSurfaceBc;

struct _GfsSurfaceBc {
  /*< private >*/
  GfsSurfaceGenericBc parent;

  /*< public >*/
  GfsFunction * type, * val;
};

#define GFS_SURFACE_BC(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSurfaceBc,\
					         gfs_surface_bc_class ())
#define GFS_IS_SURFACE_BC(obj)         (gts_object_is_from_class (obj,\
						 gfs_surface_bc_class ()))

GfsSurfaceGenericBcClass * gfs_surface_bc_class  (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __TIMESTEP_H__ */
