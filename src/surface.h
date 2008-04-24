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

#ifndef __SURFACE_H__
#define __SURFACE_H__

#include <gts.h>
#include "ftt.h"

/* GfsGenericSurface: Header */

typedef GtsObject GfsGenericSurface;

typedef struct {
  GtsPoint * E, * D;
  gdouble x;
  guint n;
  gint inside;
} GfsSegment;

typedef struct _GfsGenericSurfaceClass    GfsGenericSurfaceClass;

struct _GfsGenericSurfaceClass {
  /*< private >*/
  GtsObjectClass parent_class;

  /*< public >*/
  GfsGenericSurface * (* cell_is_cut)          (FttCell * cell,
						GfsGenericSurface * s,
						gboolean flatten,
						gint maxlevel);
  guint               (* segment_intersection) (GfsGenericSurface * s,
						FttCell * cell,
						GfsSegment * I);
  void                (* segment_normal)       (GfsGenericSurface * s,
						FttCell * cell,
						GfsSegment * I,
						GtsVector n);
};

#define GFS_GENERIC_SURFACE(obj)            GTS_OBJECT_CAST (obj,\
					         GtsObject,\
					         gfs_generic_surface_class ())
#define GFS_GENERIC_SURFACE_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsGenericSurfaceClass,\
						 gfs_generic_surface_class())
#define GFS_IS_GENERIC_SURFACE(obj)         (gts_object_is_from_class (obj,\
						 gfs_generic_surface_class ()))

GfsGenericSurfaceClass * gfs_generic_surface_class  (void);
guint              gfs_surface_segment_intersection (GfsGenericSurface * s,
						     FttCell * cell,
						     GfsSegment * I);
void               gfs_surface_segment_normal       (GfsGenericSurface * s,
						     FttCell * cell,
						     GfsSegment * I,
						     GtsVector n);
GfsGenericSurface *      gfs_cell_is_cut (FttCell * cell,
					  GfsGenericSurface * s,
					  gboolean flatten,
					  gint maxlevel);
typedef void       (* FttCellTraverseCutFunc) (FttCell * cell,
					       GfsGenericSurface * s,
					       gpointer data);
void               gfs_cell_traverse_cut       (FttCell * root,
						GfsGenericSurface * s,
						FttTraverseType order,
						FttTraverseFlags flags,
						FttCellTraverseCutFunc func,
						gpointer data);
void               gfs_cell_traverse_cut_2D    (FttCell * root,
						GfsGenericSurface * s,
						FttTraverseType order,
						FttTraverseFlags flags,
						FttCellTraverseCutFunc func,
						gpointer data);
void               gfs_generic_surface_read    (GfsGenericSurface * s, 
						gpointer sim,
						GtsFile * fp);
void               gfs_generic_surface_write   (GfsGenericSurface * s,
						gpointer sim,
						FILE * fp);

/* GfsSurface: Header */

typedef struct _GfsSurface         GfsSurface;

struct _GfsSurface {
  /*< private >*/
  GtsObject parent;
  GtsVector rotate, scale, translate;
  gboolean flip;
  GfsFunction * f;
  GtsMatrix * m;

  /*< public >*/
  GtsSurface * s;
  gboolean twod;
};

#define GFS_SURFACE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSurface,\
					         gfs_surface_class ())
#define GFS_IS_SURFACE(obj)         (gts_object_is_from_class (obj,\
						 gfs_surface_class ()))

GfsGenericSurfaceClass *   gfs_surface_class          (void);
gdouble            gfs_surface_implicit_value (GfsSurface * s, 
					       GtsPoint p);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __SURFACE_H__ */
