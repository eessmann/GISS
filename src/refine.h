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

#ifndef __REFINE_H__
#define __REFINE_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "simulation.h"

/* GfsRefine: Header */

typedef struct _GfsRefine             GfsRefine;
typedef struct _GfsRefineClass        GfsRefineClass;

struct _GfsRefine {
  GtsSListContainee parent;

  GfsFunction * maxlevel;
};

struct _GfsRefineClass {
  GtsSListContaineeClass parent_class;

  void (* refine) (GfsRefine * refine, GfsSimulation * simulation);
};

#define GFS_REFINE(obj)            GTS_OBJECT_CAST (obj,\
					           GfsRefine,\
					           gfs_refine_class ())
#define GFS_REFINE_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						   GfsRefineClass,\
						   gfs_refine_class())
#define GFS_IS_REFINE(obj)         (gts_object_is_from_class (obj,\
						   gfs_refine_class ()))
     
GfsRefineClass * gfs_refine_class  (void);
GfsRefine *      gfs_refine_new    (GfsRefineClass * klass);

/* GfsRefineSolid: Header */

typedef struct _GfsRefineSolid         GfsRefineSolid;
typedef struct _GfsRefineSolidClass    GfsRefineSolidClass;

struct _GfsRefineSolid {
  GfsRefine parent;
};

struct _GfsRefineSolidClass {
  GfsRefineClass parent_class;
};

#define GFS_REFINE_SOLID(obj)            GTS_OBJECT_CAST (obj,\
					           GfsRefineSolid,\
					           gfs_refine_solid_class ())
#define GFS_REFINE_SOLID_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						   GfsRefineSolidClass,\
						   gfs_refine_solid_class())
#define GFS_IS_REFINE_SOLID(obj)         (gts_object_is_from_class (obj,\
						   gfs_refine_solid_class ()))
     
GfsRefineSolidClass * gfs_refine_solid_class  (void);

/* GfsRefineDistance: Header */

typedef struct _GfsRefineDistance         GfsRefineDistance;
typedef struct _GfsRefineDistanceClass    GfsRefineDistanceClass;

struct _GfsRefineDistance {
  GfsRefine parent;

  GtsSurface * surface;
  GNode * stree;
};

struct _GfsRefineDistanceClass {
  GfsRefineClass parent_class;
};

#define GFS_REFINE_DISTANCE(obj)            GTS_OBJECT_CAST (obj,\
					          GfsRefineDistance,\
					          gfs_refine_distance_class ())
#define GFS_REFINE_DISTANCE_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						  GfsRefineDistanceClass,\
						  gfs_refine_distance_class())
#define GFS_IS_REFINE_DISTANCE(obj)         (gts_object_is_from_class (obj,\
						  gfs_refine_distance_class ()))
     
GfsRefineDistanceClass * gfs_refine_distance_class  (void);
GfsRefineDistance *      gfs_refine_distance_new (GfsRefineDistanceClass * klass);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __REFINE_H__ */
