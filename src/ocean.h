/* Gerris - The GNU Flow Solver
 * Copyright (C) 2004 Stéphane Popinet
 * National Institute of Water and Atmospheric Research
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

#ifndef __OCEAN_H__
#define __OCEAN_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "simulation.h"
#include "source.h"

/* GfsOcean: Header */

typedef struct _GfsOcean         GfsOcean;

struct _GfsOcean {
  GfsSimulation parent;

  GPtrArray * layer;
  GfsDomain * toplayer;
};

#define GFS_OCEAN(obj)            GTS_OBJECT_CAST (obj,\
					         GfsOcean,\
					         gfs_ocean_class ())
#define GFS_IS_OCEAN(obj)         (gts_object_is_from_class (obj,\
						 gfs_ocean_class ()))

GfsSimulationClass * gfs_ocean_class          (void);
void                 gfs_hydrostatic_pressure (GfsDomain * domain,
					       GfsVariable * p,
					       GfsVariable * rho,
					       gdouble g);

/* GfsSourceHydrostatic: Header */

typedef struct _GfsSourceHydrostatic         GfsSourceHydrostatic;

struct _GfsSourceHydrostatic {
  /*< private >*/
  GfsSourceVelocity parent;
  GfsVariable * ph1;
  gboolean not_first;

  /*< public >*/
  GfsVariable * ph, * rho;
};

#define GFS_SOURCE_HYDROSTATIC(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceHydrostatic,\
					         gfs_source_hydrostatic_class ())
#define GFS_IS_SOURCE_HYDROSTATIC(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_hydrostatic_class ()))

GfsSourceGenericClass * gfs_source_hydrostatic_class    (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __OCEAN_H__ */
