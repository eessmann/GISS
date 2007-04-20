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

/* GfsSurface: Header */

typedef struct _GfsSurface         GfsSurface;

struct _GfsSurface {
  /*< private >*/
  GtsObject parent;
  GtsVector rotate, scale, translate;
  gdouble angle;
  gboolean flip;

  /*< public >*/
  GtsSurface * s;
  GfsFunction * f;
  GtsMatrix * m;
};

typedef struct {
  GtsPoint * E, * D;
  gdouble x;
  guint n;
  gint inside;
} GfsSegment;

#define GFS_SURFACE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSurface,\
					         gfs_surface_class ())
#define GFS_IS_SURFACE(obj)         (gts_object_is_from_class (obj,\
						 gfs_surface_class ()))

GtsObjectClass *   gfs_surface_class          (void);
void               gfs_surface_read           (GfsSurface * s, 
					       gpointer sim,
					       GtsFile * fp);
void               gfs_surface_write          (GfsSurface * s,
					       gpointer sim,
					       FILE * fp);
gdouble            gfs_surface_implicit_value (GfsSurface * s, 
					       GtsPoint p);
guint              gfs_surface_segment_intersection (GfsSurface * s,
						     GfsSegment * I);
GfsSurface *       gfs_cell_is_cut            (FttCell * cell,
					       GfsSurface * s,
					       gboolean flatten);
typedef void       (* FttCellTraverseCutFunc) (FttCell * cell,
					       GfsSurface * s,
					       gpointer data);
void               gfs_cell_traverse_cut       (FttCell * root,
						GfsSurface * s,
						FttTraverseType order,
						FttTraverseFlags flags,
						FttCellTraverseCutFunc func,
						gpointer data);
void               gfs_cell_traverse_cut_2D    (FttCell * root,
						GfsSurface * s,
						FttTraverseType order,
						FttTraverseFlags flags,
						FttCellTraverseCutFunc func,
						gpointer data);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __SURFACE_H__ */
