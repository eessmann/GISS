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

#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "event.h"

/* GfsOutput: Header */

typedef struct _GfsOutput         GfsOutput;
typedef struct _GfsOutputClass    GfsOutputClass;
typedef struct _GfsOutputFile     GfsOutputFile;

struct _GfsOutput {
  GfsEvent parent;

  GfsOutputFile * file;
  gchar * format;
  GSList * formats;
  gboolean dynamic, first_call;
};

struct _GfsOutputClass {
  GfsEventClass parent_class;
};

#define GFS_OUTPUT(obj)            GTS_OBJECT_CAST (obj,\
					           GfsOutput,\
					           gfs_output_class ())
#define GFS_OUTPUT_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						   GfsOutputClass,\
						   gfs_output_class())
#define GFS_IS_OUTPUT(obj)         (gts_object_is_from_class (obj,\
						   gfs_output_class ()))
     
GfsOutputClass * gfs_output_class  (void);
void             gfs_output_mute   (GfsOutput * output);

struct _GfsOutputFile {
  guint refcount;
  gchar * name;
  FILE * fp;
};

GfsOutputFile * gfs_output_file_open    (const gchar * name,
					 const gchar * mode);
void            gfs_output_file_close   (GfsOutputFile * file);

/* GfsOutputTime: Header */

GfsOutputClass * gfs_output_time_class  (void);

/* GfsOutputProgress: Header */

GfsOutputClass * gfs_output_progress_class  (void);

/* GfsOutputProjectionStats: Header */

GfsOutputClass * gfs_output_projection_stats_class  (void);

/* GfsOutputDiffusionStats: Header */

GfsOutputClass * gfs_output_diffusion_stats_class  (void);

/* GfsOutputSolidStats: Header */

GfsOutputClass * gfs_output_solid_stats_class  (void);

/* GfsOutputAdaptStats: Header */

GfsOutputClass * gfs_output_adapt_stats_class  (void);

/* GfsOutputTiming: Header */

GfsOutputClass * gfs_output_timing_class (void);

/* GfsOutputBalance: Header */

GfsOutputClass * gfs_output_balance_class  (void);

/* GfsOutputSolidForce: Header */

GfsOutputClass * gfs_output_solid_force_class (void);

/* GfsOutputLocation: Header */

typedef struct _GfsOutputLocation         GfsOutputLocation;

struct _GfsOutputLocation {
  /*< private >*/
  GfsOutput parent;

  /*< public >*/
  FttVector p;
};

#define GFS_OUTPUT_LOCATION(obj)            GTS_OBJECT_CAST (obj,\
					         GfsOutputLocation,\
					         gfs_output_location_class ())
#define GFS_IS_OUTPUT_LOCATION(obj)         (gts_object_is_from_class (obj,\
						 gfs_output_location_class ()))

GfsOutputClass * gfs_output_location_class  (void);

/* GfsOutputSimulation: Header */

typedef struct _GfsOutputSimulation         GfsOutputSimulation;

struct _GfsOutputSimulation {
  GfsOutput parent;

  gint max_depth;
  GfsVariable * var;
};

#define GFS_OUTPUT_SIMULATION(obj)            GTS_OBJECT_CAST (obj,\
					     GfsOutputSimulation,\
					     gfs_output_simulation_class ())
#define GFS_OUTPUT_SIMULATION_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
					     GfsOutputClass,\
					     gfs_output_simulation_class())
#define GFS_IS_OUTPUT_SIMULATION(obj)         (gts_object_is_from_class (obj,\
					     gfs_output_simulation_class ()))
     
GfsOutputClass * gfs_output_simulation_class  (void);

/* GfsOutputBoundaries: Header */

GfsOutputClass * gfs_output_boundaries_class  (void);

/* GfsOutputScalar: Header */

typedef struct _GfsOutputScalar         GfsOutputScalar;

struct _GfsOutputScalar {
  /*< private >*/
  GfsOutput parent;

  gboolean autoscale;
  
  /*< public >*/
  GfsVariable * v;
  gdouble min, max;
  gint maxlevel;
  GtsBBox * box;
};

#define GFS_OUTPUT_SCALAR(obj)            GTS_OBJECT_CAST (obj,\
					         GfsOutputScalar,\
					         gfs_output_scalar_class ())
#define GFS_IS_OUTPUT_SCALAR(obj)         (gts_object_is_from_class (obj,\
						 gfs_output_scalar_class ()))

GfsOutputClass * gfs_output_scalar_class  (void);

/* GfsOutputScalarNorm: Header */

GfsOutputClass * gfs_output_scalar_norm_class  (void);

/* GfsOutputScalarStats: Header */

GfsOutputClass * gfs_output_scalar_stats_class  (void);

/* GfsOutputScalarSum: Header */

GfsOutputClass * gfs_output_scalar_sum_class  (void);

/* GfsOutputEnergy: Header */

GfsOutputClass * gfs_output_energy_class  (void);

/* GfsOutputErrorNorm: Header */

typedef struct _GfsOutputErrorNorm        GfsOutputErrorNorm;

struct _GfsOutputErrorNorm {
  /*< private >*/
  GfsOutputScalar parent;
  GfsVariable * v;
  
  /*< public >*/
  GfsFunction * s;
  gboolean unbiased;
};

#define GFS_OUTPUT_ERROR_NORM(obj)            GTS_OBJECT_CAST (obj,\
					         GfsOutputErrorNorm,\
					         gfs_output_error_norm_class ())
#define GFS_IS_OUTPUT_ERROR_NORM(obj)         (gts_object_is_from_class (obj,\
						 gfs_output_error_norm_class ()))

GfsOutputClass * gfs_output_error_norm_class  (void);

/* GfsOutputCorrelation: Header */

GfsOutputClass * gfs_output_correlation_class  (void);

/* GfsOutputSquares: Header */

#define GFS_IS_OUTPUT_SQUARES(obj)         (gts_object_is_from_class (obj,\
						 gfs_output_squares_class ()))

GfsOutputClass * gfs_output_squares_class  (void);

/* GfsOutputStreamline: Header */

typedef struct _GfsOutputStreamline         GfsOutputStreamline;

struct _GfsOutputStreamline {
  /*< private >*/
  GfsOutputScalar parent;

  /*< public >*/
  FttVector p;
};

#define GFS_OUTPUT_STREAMLINE(obj)         GTS_OBJECT_CAST (obj,\
					       GfsOutputStreamline,\
					       gfs_output_streamline_class ())
#define GFS_IS_OUTPUT_STREAMLINE(obj)     (gts_object_is_from_class (obj,\
					       gfs_output_streamline_class ()))

GfsOutputClass * gfs_output_streamline_class  (void);

/* GfsOutputStreakline: Header */

typedef struct _GfsOutputStreakline         GfsOutputStreakline;

struct _GfsOutputStreakline {
  /*< private >*/
  GfsOutputScalar parent;

  gboolean started;
  gdouble ds;

  /*< public >*/
  GtsPoint * p;
  GSList * streak;
};

#define GFS_OUTPUT_STREAKLINE(obj)         GTS_OBJECT_CAST (obj,\
					       GfsOutputStreakline,\
					       gfs_output_streakline_class ())
#define GFS_IS_OUTPUT_STREAKLINE(obj)     (gts_object_is_from_class (obj,\
					       gfs_output_streakline_class ()))

GfsOutputClass * gfs_output_streakline_class  (void);

/* GfsOutputPPM: Header */

#define GFS_IS_OUTPUT_PPM(obj)         (gts_object_is_from_class (obj,\
						 gfs_output_ppm_class ()))

GfsOutputClass * gfs_output_ppm_class  (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __OUTPUT_H__ */
