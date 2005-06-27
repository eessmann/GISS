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

#ifndef __ADAPTIVE_H__
#define __ADAPTIVE_H__

#include "simulation.h"
#include "event.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void          gfs_cell_coarse_init          (FttCell * cell,
					     GfsDomain * domain);
void          gfs_cell_fine_init            (FttCell * cell,
					     GfsDomain * domain);
void          gfs_adapt_stats_init          (GfsAdaptStats * s);
void          gfs_adapt_stats_update        (GfsAdaptStats * s);
void          gfs_simulation_adapt          (GfsSimulation * simulation,
					     GfsAdaptStats * s);

/* GfsAdapt: Header */

typedef struct _GfsAdapt         GfsAdapt;

struct _GfsAdapt {
  /*< private >*/
  GfsEvent parent;
  gboolean active;

  /*< public >*/
  GfsFunction * minlevel, * maxlevel;
  guint mincells, maxcells;
  gdouble cmax, weight;
  GfsVariable * c;
  GtsKeyFunc cost;
};

#define GFS_ADAPT(obj)            GTS_OBJECT_CAST (obj,\
					         GfsAdapt,\
					         gfs_adapt_class ())
#define GFS_IS_ADAPT(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_class ()))

GfsEventClass * gfs_adapt_class  (void);

/* GfsAdaptVorticity: Header */

typedef struct _GfsAdaptVorticity         GfsAdaptVorticity;

struct _GfsAdaptVorticity {
  /*< private >*/
  GfsAdapt parent;
  GfsVariable ** u;
  gdouble maxa;

  /*< public >*/
};

#define GFS_ADAPT_VORTICITY(obj)            GTS_OBJECT_CAST (obj,\
					         GfsAdaptVorticity,\
					         gfs_adapt_vorticity_class ())
#define GFS_IS_ADAPT_VORTICITY(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_vorticity_class ()))

GfsEventClass * gfs_adapt_vorticity_class  (void);
 
/* GfsAdaptStreamlineCurvature: Header */

#define GFS_IS_ADAPT_STREAMLINE_CURVATURE(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_streamline_curvature_class ()))

GfsEventClass * gfs_adapt_streamline_curvature_class  (void);
 
/* GfsAdaptFunction: Header */

typedef struct _GfsAdaptFunction         GfsAdaptFunction;

struct _GfsAdaptFunction {
  /*< private >*/
  GfsAdapt parent;

  /*< public >*/
  GfsFunction * f;
};

#define GFS_ADAPT_FUNCTION(obj)            GTS_OBJECT_CAST (obj,\
					         GfsAdaptFunction,\
					         gfs_adapt_function_class ())
#define GFS_IS_ADAPT_FUNCTION(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_function_class ()))

GfsEventClass * gfs_adapt_function_class  (void);

/* GfsAdaptGradient: Header */

typedef struct _GfsAdaptGradient         GfsAdaptGradient;

struct _GfsAdaptGradient {
  /*< private >*/
  GfsAdapt parent;

  /*< public >*/
  GfsVariable * v;
};

#define GFS_ADAPT_GRADIENT(obj)            GTS_OBJECT_CAST (obj,\
					         GfsAdaptGradient,\
					         gfs_adapt_gradient_class ())
#define GFS_IS_ADAPT_GRADIENT(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_gradient_class ()))

GfsEventClass * gfs_adapt_gradient_class  (void);

/* GfsAdaptCurvature: Header */

#define GFS_IS_ADAPT_CURVATURE(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_curvature_class ()))

GfsEventClass * gfs_adapt_curvature_class  (void);

/* GfsAdaptNotBox: Header */

typedef struct _GfsAdaptNotBox         GfsAdaptNotBox;

struct _GfsAdaptNotBox {
  /*< private >*/
  GfsAdapt parent;
  FttVector p1, p2;

  /*< public >*/
  GfsBox * box;
};

#define GFS_ADAPT_NOT_BOX(obj)            GTS_OBJECT_CAST (obj,\
					         GfsAdaptNotBox,\
					         gfs_adapt_not_box_class ())
#define GFS_IS_ADAPT_NOT_BOX(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_not_box_class ()))

GfsEventClass *  gfs_adapt_not_box_class  (void);
GfsAdaptNotBox * gfs_adapt_not_box_new    (GfsEventClass * klass,
					   GfsBox * box);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __ADAPTIVE_H__ */
