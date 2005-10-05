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

#ifndef __SOURCE_H__
#define __SOURCE_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "event.h"

gdouble    gfs_variable_mac_source     (GfsVariable * v, 
					FttCell * cell);
void       gfs_domain_variable_centered_sources (GfsDomain * domain, 
						 GfsVariable * v,
						 GfsVariable * sv,
						 gdouble dt);
/* GfsSourceGeneric: Header */

typedef struct _GfsSourceGeneric         GfsSourceGeneric;

struct _GfsSourceGeneric {
  /*< private >*/
  GfsEvent parent;

  /*< public >*/
  GfsVariable * v;
};

typedef struct _GfsSourceGenericClass    GfsSourceGenericClass;

struct _GfsSourceGenericClass {
  /*< private >*/
  GfsEventClass parent_class;

  /*< public >*/
  gdouble (* mac_value)      (GfsSourceGeneric *, FttCell *, GfsVariable *);
  gdouble (* centered_value) (GfsSourceGeneric *, 
			      FttCell *, 
			      GfsVariable *);
};

#define GFS_SOURCE_GENERIC(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceGeneric,\
					         gfs_source_generic_class ())
#define GFS_SOURCE_GENERIC_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsSourceGenericClass,\
						 gfs_source_generic_class())
#define GFS_IS_SOURCE_GENERIC(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_generic_class ()))

GfsSourceGenericClass * gfs_source_generic_class  (void);

/* GfsSourceScalar: Header */

typedef struct _GfsSourceScalar         GfsSourceScalar;

struct _GfsSourceScalar {
  /*< private >*/
  GfsSourceGeneric parent;

  /*< public >*/
  GfsVariable * v;
};

#define GFS_SOURCE_SCALAR(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceScalar,\
					         gfs_source_scalar_class ())
#define GFS_IS_SOURCE_SCALAR(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_scalar_class ()))

GfsSourceGenericClass * gfs_source_scalar_class  (void);

/* GfsSourceVelocity: Header */

typedef struct _GfsSourceVelocity         GfsSourceVelocity;

struct _GfsSourceVelocity {
  /*< private >*/
  GfsSourceGeneric parent;

  /*< public >*/
  GfsVariable ** v;
};

#define GFS_SOURCE_VELOCITY(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceVelocity,\
					         gfs_source_velocity_class ())
#define GFS_IS_SOURCE_VELOCITY(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_velocity_class ()))

GfsSourceGenericClass * gfs_source_velocity_class  (void);

/* GfsSource: Header */

typedef struct _GfsSource         GfsSource;

struct _GfsSource {
  /*< private >*/
  GfsSourceScalar parent;

  /*< public >*/
  GfsFunction * intensity;
};

#define GFS_SOURCE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSource,\
					         gfs_source_class ())
#define GFS_IS_SOURCE(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_class ()))

GfsSourceGenericClass * gfs_source_class  (void);

/* GfsSourceControl: Header */

typedef struct _GfsSourceControl         GfsSourceControl;

struct _GfsSourceControl {
  /*< private >*/
  GfsSource parent;
  gdouble s;
};

#define GFS_SOURCE_CONTROL(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceControl,\
					         gfs_source_control_class ())
#define GFS_IS_SOURCE_CONTROL(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_control_class ()))

GfsSourceGenericClass * gfs_source_control_class  (void);

/* GfsDiffusion: Header */

typedef struct _GfsDiffusion         GfsDiffusion;

struct _GfsDiffusion {
  /*< private >*/
  GfsEvent parent;

  /*< public >*/
  GfsFunction * val;
  GfsVariable * mu;
};

typedef struct _GfsDiffusionClass    GfsDiffusionClass;

struct _GfsDiffusionClass {
  /*< private >*/
  GfsEventClass parent_class;

  /*< public >*/
  gdouble (* face)  (GfsDiffusion *, FttCellFace *);
  gdouble (* cell)  (GfsDiffusion *, FttCell *);
};

#define GFS_DIFFUSION(obj)            GTS_OBJECT_CAST (obj,\
					         GfsDiffusion,\
					         gfs_diffusion_class ())
#define GFS_DIFFUSION_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsDiffusionClass,\
						 gfs_diffusion_class())
#define GFS_IS_DIFFUSION(obj)         (gts_object_is_from_class (obj,\
						 gfs_diffusion_class ()))

GfsDiffusionClass * gfs_diffusion_class  (void);
gdouble             gfs_diffusion_face   (GfsDiffusion * d, 
					  FttCellFace * f);
gdouble             gfs_diffusion_cell   (GfsDiffusion * d, 
					  FttCell * cell);

/* GfsSourceDiffusion: Header */

struct _GfsSourceDiffusion {
  /*< private >*/
  GfsSourceScalar parent;

  /*< public >*/
  GfsDiffusion * D;
};

#define GFS_SOURCE_DIFFUSION(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceDiffusion,\
					         gfs_source_diffusion_class ())
#define GFS_IS_SOURCE_DIFFUSION(obj)         (gts_object_is_from_class (obj,\
						gfs_source_diffusion_class ()))

GfsSourceGenericClass *   gfs_source_diffusion_class  (void);
gdouble                   gfs_source_diffusion_face   (GfsSourceDiffusion * d, 
						       FttCellFace * f);
gdouble                   gfs_source_diffusion_cell   (GfsSourceDiffusion * d, 
						       FttCell * cell);

/* GfsSourceDiffusionExplicit: Header */

typedef struct _GfsSourceDiffusionExplicit         GfsSourceDiffusionExplicit;

struct _GfsSourceDiffusionExplicit {
  /*< private >*/
  GfsSourceDiffusion parent;

  /*< public >*/
  GfsVariable * s;
};

#define GFS_SOURCE_DIFFUSION_EXPLICIT(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceDiffusionExplicit,\
					         gfs_source_diffusion_explicit_class ())
#define GFS_IS_SOURCE_DIFFUSION_EXPLICIT(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_diffusion_explicit_class ()))

GfsSourceGenericClass * gfs_source_diffusion_explicit_class  (void);

/* GfsSourceViscosity: Header */

typedef struct _GfsSourceViscosity         GfsSourceViscosity;

struct _GfsSourceViscosity {
  /*< private >*/
  GfsSourceDiffusion parent;

  /*< public >*/
  GfsVariable ** v;
};

#define GFS_SOURCE_VISCOSITY(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceViscosity,\
					         gfs_source_viscosity_class ())
#define GFS_IS_SOURCE_VISCOSITY(obj) (gts_object_is_from_class (obj,\
				       gfs_source_viscosity_class ()))

GfsSourceGenericClass * gfs_source_viscosity_class  (void);

/* GfsSourceCoriolis: Header */

typedef struct _GfsSourceCoriolis         GfsSourceCoriolis;

struct _GfsSourceCoriolis {
  /*< private >*/
  GfsSourceVelocity parent;
  GfsVariable * u[2];

  /*< public >*/
  GfsFunction * omegaz, * drag;
};

#define GFS_SOURCE_CORIOLIS(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceCoriolis,\
					         gfs_source_coriolis_class ())
#define GFS_IS_SOURCE_CORIOLIS(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_coriolis_class ()))

GfsSourceGenericClass * gfs_source_coriolis_class    (void);
void                    gfs_source_coriolis_implicit (GfsDomain * domain,
						      gdouble dt);
GfsSourceCoriolis *     gfs_has_source_coriolis      (GfsDomain * domain);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __SOURCE_H__ */
