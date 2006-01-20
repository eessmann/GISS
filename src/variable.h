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

#include "timestep.h"
#include "event.h"

/* GfsVariable: Header */

typedef void (* GfsVariableFineCoarseFunc) (FttCell * cell, GfsVariable * v);

struct _GfsVariable {
  /*< private >*/
  GfsEvent parent;

  /*< public >*/
  guint i;
  FttComponent component;
  gchar * name;
  gboolean centered;
  GfsVariableFineCoarseFunc fine_coarse, coarse_fine;
  GtsContainer * sources;
  GfsSurfaceGenericBc * surface_bc;
  GfsDomain * domain;
};

typedef struct _GfsVariableClass    GfsVariableClass;

struct _GfsVariableClass {
  /*< private >*/
  GfsEventClass parent_class;

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

GfsVariableClass *    gfs_variable_class            (void);
GfsVariable *         gfs_variable_new              (GfsVariableClass * klass,
						     GfsDomain * domain,
						     const gchar * name);
#define               gfs_temporary_variable(d)     (gfs_variable_new (gfs_variable_class (),\
                                                                      (d), NULL))
GfsVariable *         gfs_variable_from_name        (GSList * i,
						     const gchar * name);
GSList *              gfs_variables_from_list       (GSList * i,
						     gchar * list,
						     gchar ** error);
#define gfs_variable_set_vector(v, c)  ((v)->component = (c))

/* GfsVariableTracer: header */

typedef struct _GfsVariableTracer                GfsVariableTracer;

struct _GfsVariableTracer {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsAdvectionParams advection;
};

#define GFS_VARIABLE_TRACER(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableTracer,\
					           gfs_variable_tracer_class ())
#define GFS_IS_VARIABLE_TRACER(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_tracer_class ()))

GfsVariableClass * gfs_variable_tracer_class  (void);

/* GfsVariableResidual: header */

#define GFS_IS_VARIABLE_RESIDUAL(obj)         (gts_object_is_from_class (obj,\
					       gfs_variable_residual_class ()))

GfsVariableClass * gfs_variable_residual_class  (void);

/* GfsVariableFiltered: header */

typedef struct _GfsVariableFiltered                GfsVariableFiltered;

struct _GfsVariableFiltered {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsVariable * v;
  guint niter;
};

#define GFS_VARIABLE_FILTERED(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableFiltered,\
					           gfs_variable_filtered_class ())
#define GFS_IS_VARIABLE_FILTERED(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_filtered_class ()))

GfsVariableClass * gfs_variable_filtered_class  (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __VARIABLE_H__ */
