/* Gerris - The GNU Flow Solver			(-*-C-*-)
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
 * Author of the object: Gaurav Tomar
 */
#include "particle.h"
#include "config.h"

/* ParticleLagrange */
typedef struct _GfsParticulate GfsParticulate;
struct _GfsParticulate{
  GfsParticle parent;
  FttVector vel;
  gdouble mass, volume;
  FttVector force;
  GtsSListContainer * forces;
};

#define GFS_PARTICULATE(obj)            GTS_OBJECT_CAST (obj,		\
							 GfsParticulate, gfs_particulate_class ())
#define GFS_IS_PARTICULATE(obj)         (gts_object_is_from_class (obj, gfs_particulate_class ()))

GfsEventClass * gfs_particulate_class  (void);

/* ParticleList */
typedef struct _GfsParticleList GfsParticleList;
struct _GfsParticleList{
  GfsEventList parent;
  GtsSListContainer * forces;
};

#define GFS_PARTICLE_LIST(obj)            GTS_OBJECT_CAST (obj,		\
							   GfsParticleList, \
							   gfs_particle_list_class ())

#define GFS_IS_PARTICLE_LIST(obj)         (gts_object_is_from_class (obj, \
								     gfs_particle_list_class ()))

GfsEventClass * gfs_particle_list_class  (void);

/* Particle Force */
typedef struct _GfsParticleForce GfsParticleForce;
struct _GfsParticleForce{
  GtsSListContainee parent;
  FttVector (* force) (GfsParticle *p, GfsParticleForce *force);
};

#define PARTICLE_FORCE(obj)            GTS_OBJECT_CAST (obj,		\
							GfsParticleForce, \
							gfs_particle_force_class ())
#define GFS_IS_PARTICLE_FORCE(obj)         (gts_object_is_from_class (obj, \
								      gfs_particle_force_class ()))

GtsSListContaineeClass * gfs_particle_force_class  (void);

/* Lift Force */
typedef struct _ForceLift ForceLift;
struct _ForceLift{
  GfsParticleForce parent;
  GfsFunction * coefficient;
  GfsVariable *re_p, *u_rel, *v_rel, *w_rel, *pdia;
  GfsParticulate *p;
};

#define FORCE_LIFT(obj)            GTS_OBJECT_CAST (obj,		\
						    ForceLift,		\
						    gfs_force_lift_class ())
#define GFS_IS_FORCE_LIFT(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_lift_class ()))
GtsSListContaineeClass * gfs_force_lift_class  (void);

/* Drag Force */
typedef struct _ForceDrag ForceDrag;
struct _ForceDrag{
  GfsParticleForce parent;
  GfsFunction * coefficient;
  GfsVariable *re_p, *u_rel, *v_rel, *w_rel, *pdia;
  GfsParticulate *p;
};

#define FORCE_DRAG(obj)            GTS_OBJECT_CAST (obj,		\
						    ForceDrag,		\
						    gfs_force_drag_class ())
#define GFS_IS_FORCE_DRAG(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_drag_class ()))
GtsSListContaineeClass * gfs_force_drag_class  (void);

/* Buoy Force */
typedef struct _ForceBuoy ForceBuoy;
struct _ForceBuoy{
  GfsParticleForce parent;
  GfsParticulate *p;
};

#define FORCE_BUOY(obj)            GTS_OBJECT_CAST (obj,		\
						    ForceBuoy,		\
						    gfs_force_buoy_class ())
#define GFS_IS_FORCE_BUOY(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_buoy_class ()))
GtsSListContaineeClass * gfs_force_buoy_class  (void);
