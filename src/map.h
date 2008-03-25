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

#ifndef __MAP_H__
#define __MAP_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "ftt.h"

/* GfsMap: Header */

typedef struct _GfsMap         GfsMap;

struct _GfsMap {
  /*< private >*/
  GtsSListContainee parent;

  /*< public >*/
};

typedef struct _GfsMapClass    GfsMapClass;

struct _GfsMapClass {
  /*< private >*/
  GtsSListContaineeClass parent_class;

  /*< public >*/
  void (* transform) (GfsMap * map, const FttVector * src, FttVector * dest);
  void (* inverse)   (GfsMap * map, const FttVector * src, FttVector * dest);
};

#define GFS_MAP(obj)            GTS_OBJECT_CAST (obj,\
					         GfsMap,\
					         gfs_map_class ())
#define GFS_MAP_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsMapClass,\
						 gfs_map_class())
#define GFS_IS_MAP(obj)         (gts_object_is_from_class (obj,\
						 gfs_map_class ()))

GfsMapClass * gfs_map_class      (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MAP_H__ */
