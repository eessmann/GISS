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

#ifndef __TENSION_H__
#define __TENSION_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "source.h"

/* GfsSourceTensionCSS: Header */

typedef struct _GfsSourceTensionCSS         GfsSourceTensionCSS;

struct _GfsSourceTensionCSS {
  /*< private >*/
  GfsSourceVelocity parent;
  GfsVariable * g[3];
  
  /*< public >*/
  GfsVariable * c, * t[FTT_DIMENSION];
  gdouble sigma;
};

#define GFS_SOURCE_TENSION_CSS(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceTensionCSS,\
					         gfs_source_tension_css_class ())
#define GFS_IS_SOURCE_TENSION_CSS(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_tension_css_class ()))

GfsSourceGenericClass * gfs_source_tension_css_class (void);

/* GfsSourceTension: Header */

typedef struct _GfsSourceTension         GfsSourceTension;

struct _GfsSourceTension {
  /*< private >*/
  GfsSourceVelocity parent;
  
  /*< public >*/
  GfsVariable * c, * k;
};

#define GFS_SOURCE_TENSION(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceTension,\
					         gfs_source_tension_class ())
#define GFS_IS_SOURCE_TENSION(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_tension_class ()))

GfsSourceGenericClass * gfs_source_tension_class (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __TENSION_H__ */
