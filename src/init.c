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

#include "config.h"

#ifdef HAVE_FENV_H
# define _GNU_SOURCE
# include <fenv.h>
#endif /* HAVE_FENV_H */

#include <stdlib.h>
#include <locale.h>

#include "boundary.h"
#include "mpi_boundary.h"
#include "init.h"
#include "refine.h"
#include "output.h"
#include "adaptive.h"
#include "source.h"
#include "tension.h"
#include "ocean.h"
#include "wave.h"
#include "levelset.h"
#include "vof.h"
#include "solid.h"
#include "moving.h"
#include "river.h"
#include "balance.h"
#include "map.h"
#include "metric.h"
#include "particle.h"

#include "modules.h"

#ifdef HAVE_MPI
# include <mpi.h>
#endif /* HAVE_MPI */

static void gfs_log (const gchar * log_domain,
		     GLogLevelFlags log_level,
		     const gchar * message)
{
  int type = 0;
  gchar * pe;
  const gchar stype[][10] = {
    "ERROR", "CRITICAL", "WARNING", "MESSAGE", "INFO", "DEBUG"
  };

#ifdef HAVE_MPI
  int rank = -1;
  MPI_Comm_size (MPI_COMM_WORLD, &rank);
  if (rank > 1)
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  else
    rank = -1;
  if (rank >= 0) {
    char name[MPI_MAX_PROCESSOR_NAME];
    int length;
    MPI_Get_processor_name (name, &length);
    pe = g_strdup_printf ("PE %d (%s): ", rank, name);
  }
  else
#endif /* HAVE_MPI */
    pe = g_strdup ("");

  switch (log_level & G_LOG_LEVEL_MASK) {
  case G_LOG_LEVEL_ERROR:    type = 0; break;
  case G_LOG_LEVEL_CRITICAL: type = 1; break;
  case G_LOG_LEVEL_WARNING:  type = 2; break;
  case G_LOG_LEVEL_MESSAGE:  type = 3; break;
  case G_LOG_LEVEL_INFO:     type = 4; break;
  case G_LOG_LEVEL_DEBUG:    type = 5; break;
  default:
    g_assert_not_reached ();
  }
  fprintf (stderr, "\n%s-%s **: %s%s\n\n", 
	   log_domain, stype[type], pe, message);
  g_free (pe);
}

/**
 * gfs_classes:
 *
 * Returns: a pointer to a NULL-terminated array of all the classes
 * usable in Gerris parameter files.
 */
GtsObjectClass ** gfs_classes (void)
{
  static GtsObjectClass ** classes = NULL;
  if (classes == NULL) { gpointer klass[] = {

  gfs_global_class (),
  gfs_simulation_class (),
    gfs_ocean_class (),
    gfs_advection_class (),
    gfs_poisson_class (),
    gfs_simulation_moving_class (),
    gfs_axi_class (),
    gfs_wave_class (),
    gfs_river_class (),

  gfs_surface_bc_class (),

  gfs_box_class (),

  gfs_gedge_class (),

  gfs_bc_dirichlet_class (),
  gfs_bc_subcritical_class (),
  gfs_bc_neumann_class (),
  gfs_bc_navier_class (),
  gfs_bc_flather_class (),

  gfs_boundary_class (),
    gfs_boundary_inflow_constant_class (),
    gfs_boundary_outflow_class (),
    gfs_boundary_gradient_class (),
    gfs_boundary_periodic_class (),
      gfs_boundary_mpi_class (),

  gfs_refine_class (),
    gfs_refine_solid_class (),
    gfs_refine_surface_class (),
      gfs_refine_distance_class (),
      gfs_refine_height_class (),

  gfs_event_class (),
    gfs_variable_class (),
      gfs_variable_boolean_class (),
      gfs_variable_tracer_class (),
        gfs_variable_tracer_vof_class (),
      gfs_variable_residual_class (),
      gfs_variable_filtered_class (),
      gfs_variable_diagonal_class (),
      gfs_variable_function_class (),
#if FTT_2D
        gfs_variable_stream_function_class (),
#endif /* FTT_2D */
      gfs_variable_curvature_class (),
        gfs_variable_position_class (),
      gfs_variable_distance_class (),

    gfs_solid_class (),
      gfs_solid_moving_class(),

    gfs_init_class (),
    gfs_init_mask_class (),
    gfs_init_flow_constant_class (),
    gfs_init_fraction_class (),
#if FTT_2D
    gfs_init_vorticity_class (),
#endif /* FTT_2D */
    gfs_init_wave_class (),

    gfs_metric_lon_lat_class (),
    gfs_metric_cubed_class (),

    gfs_adapt_class (),
      gfs_adapt_vorticity_class (),
      gfs_adapt_streamline_curvature_class (),
      gfs_adapt_function_class (),
      gfs_adapt_gradient_class (),
        gfs_adapt_error_class (),

    gfs_event_sum_class (),
      gfs_event_sum_direction_class (),
    gfs_event_harmonic_class (),
    gfs_event_stop_class (),
    gfs_event_script_class (),
    gfs_event_balance_class (),
    gfs_source_generic_class (),
      gfs_source_scalar_class (),
        gfs_source_class (),
          gfs_source_control_class (),
            gfs_source_control_field_class (),
          gfs_source_flux_class (),
        gfs_source_diffusion_class (),
          gfs_source_diffusion_explicit_class (),
      gfs_source_velocity_class (),
        gfs_source_viscosity_class (),
          gfs_source_viscosity_explicit_class (),
        gfs_source_friction_class (),
        gfs_source_coriolis_class (),
          gfs_source_tension_class (),
          gfs_source_tension_css_class (),
#if !FTT_2D
        gfs_source_hydrostatic_class (),
#endif /* 2D3 or 3D */
    gfs_remove_droplets_class (),
    gfs_remove_ponds_class (),
    gfs_event_filter_class (),
    gfs_event_list_class (),
   
    gfs_output_class (),
      gfs_output_time_class (),
      gfs_output_progress_class (),
      gfs_output_projection_stats_class (),
      gfs_output_diffusion_stats_class (),
      gfs_output_solid_stats_class (),
      gfs_output_adapt_stats_class (),
      gfs_output_timing_class (),
      gfs_output_balance_class (),
      gfs_output_solid_force_class (),
      gfs_output_location_class (),
        gfs_output_particle_class (),
      gfs_output_simulation_class (),
      gfs_output_boundaries_class (),

      gfs_output_scalar_class (),
        gfs_output_scalar_norm_class (),
        gfs_output_scalar_stats_class (),
        gfs_output_scalar_sum_class (),
        gfs_output_scalar_maxima_class (),
        gfs_output_scalar_histogram_class (),
        gfs_output_droplet_sums_class (),
        gfs_output_error_norm_class (),
          gfs_output_correlation_class (),
	gfs_output_squares_class (),
	gfs_output_streamline_class (),
        gfs_output_ppm_class (),  

  gfs_map_class (),
    gfs_map_function_class (),

  gfs_particle_class (),

  NULL};

    guint n = 0;
    gpointer * c = klass;
    while (*(c++)) n++;
    classes = g_malloc ((n + 1)*sizeof (gpointer));
    memcpy (classes, klass, (n + 1)*sizeof (gpointer));
  }
  return classes;
}

typedef void (* AtExitFunc) (void);

/**
 * gfs_init:
 * @argc: a pointer on the number of command line arguments passed to
 * the program.
 * @argv: a pointer on the command line arguments passed to the
 * program.
 *
 * Initializes the Gerris library. This function must be called before
 * any other function of the library.
 */
void gfs_init (int * argc, char *** argv)
{
  static gboolean initialized = FALSE;

  if (initialized)
    return;

  if (!setlocale (LC_ALL, "POSIX"))
    g_warning ("cannot set locale to POSIX");

#ifdef HAVE_MPI
  MPI_Initialized (&initialized);
  if (!initialized) {
    if (!argc || !argv) {
      int argc1 = 1;
      char ** argv1;

      argv1 = g_malloc (sizeof (char *));
      argv1[0] = g_strdup ("gfs_init");
      MPI_Init (&argc1, &argv1);
      g_free (argv1[0]); g_free (argv1);
    }
    else
      MPI_Init (argc, argv);
    MPI_Errhandler_set (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
  }
  atexit ((AtExitFunc) MPI_Finalize);
#endif /* HAVE_MPI */
  initialized = TRUE;

#ifdef FE_NOMASK_ENV
# ifdef FE_DIVBYZERO
  feenableexcept (FE_DIVBYZERO);
# endif /* FE_DIVBYZERO */
#endif /* FE_NO_MASK_ENV */

  g_log_set_handler (G_LOG_DOMAIN,
		     G_LOG_LEVEL_ERROR |
		     G_LOG_LEVEL_CRITICAL |
		     G_LOG_LEVEL_WARNING |
		     G_LOG_LEVEL_MESSAGE |
		     G_LOG_LEVEL_INFO |
		     G_LOG_LEVEL_DEBUG |
		     G_LOG_FLAG_FATAL |
		     G_LOG_FLAG_RECURSION,
		     (GLogFunc) gfs_log, NULL);

  /* Instantiates classes before reading any domain or simulation file */
  gfs_classes ();

  /* If modules are not supported, calls modules init functions */
#include "modules.c"
}
