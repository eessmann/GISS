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

#ifndef __UTILS_H__
#define __UTILS_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <gmodule.h>
#include "ftt.h"

#define GFS_DOUBLE_TO_POINTER(d)     (*((gpointer *) &(d)))

/* GfsFunction: Header */

typedef struct _GfsFunction         GfsFunction;

struct _GfsFunction {
  /*< private >*/
  GtsObject parent;
  GModule * module;
  GString * expr;
  gdouble (* f) (gdouble x, gdouble y, gdouble z, gdouble t);
  gdouble val;

  /*< public >*/
};

typedef struct _GfsFunctionClass    GfsFunctionClass;

struct _GfsFunctionClass {
  GtsObjectClass parent_class;
};

#define GFS_FUNCTION(obj)            GTS_OBJECT_CAST (obj,\
					         GfsFunction,\
					         gfs_function_class ())
#define GFS_FUNCTION_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsFunctionClass,\
						 gfs_function_class())
#define GFS_IS_FUNCTION(obj)         (gts_object_is_from_class (obj,\
						 gfs_function_class ()))

GfsFunctionClass * gfs_function_class       (void);
GfsFunction *      gfs_function_new         (GfsFunctionClass * klass,
					     gdouble val);
gdouble            gfs_function_face_value  (GfsFunction * f, 
					     FttCellFace * fa,
					     gdouble t);
gdouble            gfs_function_value       (GfsFunction * f, 
					     FttVector * p,
					     gdouble t);
void               gfs_function_read        (GfsFunction * f, 
					     GtsFile * fp);
void               gfs_function_write       (GfsFunction * f, 
					     FILE * fp);
GtsObjectClass *   gfs_object_class_from_name (const gchar * name);

void               gfs_eigenvalues          (gdouble a[FTT_DIMENSION][FTT_DIMENSION],
					     gdouble d[FTT_DIMENSION],
					     gdouble v[FTT_DIMENSION][FTT_DIMENSION]);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __UTILS_H__ */

