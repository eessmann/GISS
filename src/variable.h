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

#ifndef __VARIABLE_H__
#define __VARIABLE_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef struct _GfsSurfaceGenericBc        GfsSurfaceGenericBc;

#include "advection.h"
#include "timestep.h"

/* GfsVariable: Header */

typedef void (* GfsVariableDerivedFunc)    (FttCell * cell, GfsVariable * v);
typedef void (* GfsVariableFineCoarseFunc) (FttCell * cell, GfsVariable * v);

struct _GfsVariable {
  /*< private >*/
  GtsObject parent;

  /*< public >*/
  guint i;
  gchar * name;
  gboolean centered;
  GfsVariableDerivedFunc derived;
  GfsVariableFineCoarseFunc fine_coarse;
  GtsContainer * sources;
  GfsSurfaceGenericBc * surface_bc;
  GfsVariable * next, * permanent;
  GtsObject * p;
};

typedef struct _GfsVariableClass    GfsVariableClass;

struct _GfsVariableClass {
  /*< private >*/
  GtsObjectClass parent_class;

  /*< public >*/
};

#define GFS_VARIABLE1(obj)            GTS_OBJECT_CAST (obj,\
					         GfsVariable,\
					         gfs_variable_class ())
#define GFS_VARIABLE_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsVariableClass,\
						 gfs_variable_class())
#define GFS_IS_VARIABLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_variable_class ()))
#define gfs_variable_parent(v)         ((v)->p)
#define gfs_variable_set_parent(v, pa) ((v)->p = GTS_OBJECT (pa))

GfsVariableClass *    gfs_variable_class            (void);
GfsVariable *         gfs_variable_new              (GfsVariableClass * klass,
						     GtsObject * parent,
						     const gchar * name,
						     gboolean centered,
						     guint i);
GfsVariable *         gfs_variable_list_copy        (GfsVariable * v,
						     GtsObject * parent);
void                  gfs_variable_list_destroy     (GfsVariable * v);
GfsVariable *         gfs_variable_from_name        (GfsVariable * variables,
						     const gchar * name);
GfsVariable *         gfs_variables_from_list       (GfsVariable * variables,
						     gchar * list,
						     gchar ** error);

/* GfsVariableTracer: header */

typedef struct _GfsVariableTracer                GfsVariableTracer;

struct _GfsVariableTracer {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsAdvectionParams advection;
  GfsMultilevelParams diffusion;
};

#define GFS_VARIABLE_TRACER(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableTracer,\
					           gfs_variable_tracer_class ())
#define GFS_IS_VARIABLE_TRACER(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_tracer_class ()))

GfsVariableClass * gfs_variable_tracer_class  (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __VARIABLE_H__ */
