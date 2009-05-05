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

#ifndef __MOVING_H__
#define __MOVING_H__

#include <gts.h>
#include "variable.h"
#include "utils.h"
#include "solid.h"
#include "advection.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef struct _GfsSolidMoving         GfsSolidMoving;

struct _GfsSolidMoving {
  /*< private >*/
  GfsSolid parent;

  /*< public >*/
  GfsFunction * level;
  gboolean active;
};

GfsEventClass * gfs_solid_moving_class (void);

#define GFS_SOLID_MOVING(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSolidMoving,\
					         gfs_solid_moving_class ())
#define GFS_IS_SOLID_MOVING(obj)         (gts_object_is_from_class (obj,\
						 gfs_solid_moving_class ()))

/* GfsMovingSimulation: Header */

typedef struct _GfsMovingSimulation         GfsMovingSimulation;

struct _GfsMovingSimulation {
  /*< private >*/
  GfsSimulation parent;

  /*< public >*/
  GfsVariable * old_solid, ** sold2;
};

#define GFS_MOVING_SIMULATION(obj)            GTS_OBJECT_CAST (obj,\
					         GfsMovingSimulation,\
					         gfs_moving_simulation_class ())
#define GFS_IS_MOVING_SIMULATION(obj)         (gts_object_is_from_class (obj,\
						 gfs_moving_simulation_class ()))

GfsSimulationClass * gfs_moving_simulation_class            (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MOVING_H__ */
