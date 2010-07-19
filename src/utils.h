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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#  ifdef HAVE_MPI
#    include <mpi.h>

# define gfs_all_reduce(domain, p, type, op) {			        \
    if ((domain)->pid >= 0) {						\
      union { int a; float b; double c;} global;			\
      MPI_Allreduce (&(p), &global, 1, type, op, MPI_COMM_WORLD);	\
      memcpy (&(p), &global, sizeof (p));				\
    }									\
  }

# define gfs_error(pid, ...) {	                                        \
    int rank;                                                           \
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);                              \
    if (rank == (pid))		             				\
      fprintf (stderr, __VA_ARGS__);                                    \
    else if ((pid) < 0) {						\
      int size;                                                         \
      MPI_Comm_size (MPI_COMM_WORLD, &size);                            \
      if (size > 1) {			                                \
        char name[MPI_MAX_PROCESSOR_NAME];                              \
        MPI_Get_processor_name (name, &size);                           \
        fprintf (stderr, "PE %d (%s): ", rank, name);			\
      }                                                                 \
      fprintf (stderr, __VA_ARGS__);                                    \
    }                                                                   \
  }

#  else /* doesn't HAVE_MPI */
    /* gfs_all_reduce() defaults to nothing without MPI */
#    define gfs_all_reduce(domain, p, type, op)
#    define gfs_error(pid, ...) fprintf(stderr, __VA_ARGS__)
#  endif /* doesn't HAVE_MPI */
#endif /* HAVE_CONFIG_H */

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

#define GFS_DOUBLE_TO_POINTER(d)     (*((gpointer *) &(d)))
#define GFS_NODATA                   G_MAXDOUBLE
#define GFS_HAS_DATA(cell,v)         (GFS_VALUE (cell, v) != GFS_NODATA)

gboolean gfs_char_in_string (char c, const char * s);
gchar *  gfs_file_statement (GtsFile * fp);

/* GfsGlobal: Header */

typedef struct _GfsGlobal         GfsGlobal;

#define GFS_GLOBAL(obj)            GTS_OBJECT_CAST (obj,\
					         GfsGlobal,\
					         gfs_global_class ())
#define GFS_IS_GLOBAL(obj)         (gts_object_is_from_class (obj,\
						 gfs_global_class ()))

GtsObjectClass * gfs_global_class  (void);

/* GfsFunction: Header */

typedef struct _GfsFunction         GfsFunction;

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
GfsFunction *      gfs_function_new_from_variable (GfsFunctionClass * klass, 
						   GfsVariable * v);
void               gfs_function_set_units   (GfsFunction * f,
					     gdouble units);
gchar *            gfs_function_description (GfsFunction * f,
					     gboolean truncate);
gdouble            gfs_function_face_value  (GfsFunction * f,
					     FttCellFace * fa);
gdouble            gfs_function_value       (GfsFunction * f,
					     FttCell * cell);
void               gfs_function_set_constant_value (GfsFunction * f, 
						    gdouble val);
gdouble            gfs_function_get_constant_value (GfsFunction * f);
GfsVariable *      gfs_function_get_variable (GfsFunction * f);
void               gfs_function_read        (GfsFunction * f, 
					     gpointer domain,
					     GtsFile * fp);
void               gfs_function_write       (GfsFunction * f, 
					     FILE * fp);
GString *          gfs_function_expression  (GtsFile * fp, 
					     gboolean * is_expression);

/* GfsFunctionSpatial: Header */

#define GFS_IS_FUNCTION_SPATIAL(obj)         (gts_object_is_from_class (obj,\
					      gfs_function_spatial_class ()))

typedef gdouble (* GfsFunctionSpatialFunc)    (double x, double y, double z, double t);

GfsFunctionClass * gfs_function_spatial_class (void);
GfsFunction *      gfs_function_spatial_new   (GfsFunctionClass * klass, 
					       GfsFunctionSpatialFunc func);
gdouble            gfs_function_spatial_value (GfsFunction * f, const FttVector * p);

/* GfsFunctionMap: Header */

#define GFS_IS_FUNCTION_MAP(obj)         (gts_object_is_from_class (obj,\
					      gfs_function_map_class ()))

GfsFunctionClass * gfs_function_map_class (void);

/* GfsFunctionConstant: Header */

#define GFS_IS_FUNCTION_CONSTANT(obj)         (gts_object_is_from_class (obj,\
					       gfs_function_constant_class ()))

GfsFunctionClass * gfs_function_constant_class (void);
gdouble            gfs_read_constant           (GtsFile * fp,
						gpointer domain);

GtsObjectClass *   gfs_object_class_from_name (const gchar * name);

void               gfs_eigenvalues          (gdouble a[FTT_DIMENSION][FTT_DIMENSION],
					     gdouble d[FTT_DIMENSION],
					     gdouble v[FTT_DIMENSION][FTT_DIMENSION]);
gdouble            gfs_matrix_inverse       (gdouble ** m, 
					     guint n,
					     gdouble pivmin);
gpointer           gfs_matrix_new           (guint n, 
					     guint p,
					     guint size);
void               gfs_matrix_free          (gpointer m);

typedef struct {
  gboolean started;
  glong start, stop;
} GfsClock;

GfsClock *         gfs_clock_new            (void);
void               gfs_clock_start          (GfsClock * t);
void               gfs_clock_stop           (GfsClock * t);
gdouble            gfs_clock_elapsed        (GfsClock * t);
void               gfs_clock_destroy        (GfsClock * t);

FILE *             gfs_union_open           (FILE * fp, 
					     int rank);
void               gfs_union_close          (FILE * fp, 
					     int rank, 
					     FILE * fpp);

/* GfsFormat: Header */

typedef struct _GfsFormat GfsFormat;

typedef enum {
  ITER,
  TIME,
  PID,
  NONE
} GfsFormatType;

struct _GfsFormat {
  gchar * s;
  GfsFormatType t;
};

GfsFormat *        gfs_format_new          (gchar * s,
			                    guint len, 
					    GfsFormatType t);

void               gfs_format_destroy      (GfsFormat * f);

gchar *            gfs_format_string       (GSList * list, 
					    gint pid, 
					    guint niter,
					    gdouble time);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __UTILS_H__ */

