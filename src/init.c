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

#ifdef HAVE_FPU_SETCW
# include <fpu_control.h>
  static fpu_control_t fpu_trap_exceptions = 
       _FPU_IEEE & ~(_FPU_MASK_ZM /*| _FPU_MASK_IM | _FPU_MASK_OM*/);
#endif /* HAVE_FPU_SETCW */

#include <stdlib.h>
#include "boundary.h"
#include "init.h"
#include "refine.h"
#include "output.h"
#include "adaptive.h"
#include "source.h"
#include "tension.h"
#include "ocean.h"

#include "modules.h"

#ifdef HAVE_MPI
# include <mpi.h>
# include "mpi_boundary.h"
#endif /* HAVE_MPI */

static void gfs_log (const gchar * log_domain,
		     GLogLevelFlags log_level,
		     const gchar * message)
{
  int rank = -1, type = 0;
  gchar pe[10];
  const gchar stype[][10] = {
    "ERROR", "CRITICAL", "WARNING", "MESSAGE", "INFO", "DEBUG"
  };

#ifdef HAVE_MPI
  MPI_Comm_size (MPI_COMM_WORLD, &rank);
  if (rank > 1)
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  else
    rank = -1;
#endif /* HAVE_MPI */
  if (rank >= 0)
    sprintf (pe, "PE %d: ", rank);
  else
    pe[0] = '\0';

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
}

static void cell_vorticity (FttCell * cell, GfsVariable * v)
{
  GFS_VARIABLE (cell, v->i) = gfs_vorticity (cell);
}

static void cell_velocity_norm (FttCell * cell, GfsVariable * v)
{
  GFS_VARIABLE (cell, v->i) = gfs_velocity_norm (cell);
}

static void cell_velocity_norm2 (FttCell * cell, GfsVariable * v)
{
  GFS_VARIABLE (cell, v->i) = gfs_velocity_norm2 (cell);
}

static void cell_level (FttCell * cell, GfsVariable * v)
{
  GFS_VARIABLE (cell, v->i) = ftt_cell_level (cell);
}

static void cell_fraction (FttCell * cell, GfsVariable * v)
{
  GFS_VARIABLE (cell, v->i) = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
}

static void cell_lambda2 (FttCell * cell, GfsVariable * v)
{
  gdouble size = ftt_cell_size (cell);
  GFS_VARIABLE (cell, v->i) /= gfs_velocity_lambda2 (cell)/(size*size);
}

static void cell_curvature (FttCell * cell, GfsVariable * v)
{
  GFS_VARIABLE (cell, v->i) = gfs_streamline_curvature (cell)/ftt_cell_size (cell);
}

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
  guint i = 0;
  GfsVariable * v;
  static gboolean initialized = FALSE;

  if (initialized)
    return;

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
    atexit ((void (*)(void)) MPI_Finalize);
  }
#endif /* HAVE_MPI */
  initialized = TRUE;

#ifdef HAVE_FPU_SETCW
   _FPU_SETCW (fpu_trap_exceptions);
#endif /* HAVE_FPU_SETCW */

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

  /* Initialize permanent variables */
  gfs_div = v = gfs_variable_new (gfs_variable_class (), NULL, NULL, FALSE, i++);
  v->permanent = NULL;
  g_assert (v->i == GFS_DIV);
  gfs_dp = v = v->next = gfs_variable_new (gfs_variable_class (), NULL, NULL, TRUE, i++);
  v->permanent = NULL;
  g_assert (v->i  == GFS_DP);
  gfs_res = v = v->next = gfs_variable_new (gfs_variable_class (), NULL, NULL, FALSE, i++);
  v->permanent = NULL;
  g_assert (v->i == GFS_RES);
  gfs_gx = v = v->next = gfs_variable_new (gfs_variable_class (), NULL, NULL, FALSE, i++);
  v->permanent = NULL;
  g_assert (v->i  == GFS_GX);
  gfs_gy = v = v->next = gfs_variable_new (gfs_variable_class (), NULL, NULL, FALSE, i++);
  v->permanent = NULL;
  g_assert (v->i  == GFS_GY);
#if (!FTT_2D)
  gfs_gz = v = v->next = gfs_variable_new (gfs_variable_class (), NULL, NULL, FALSE, i++);
  v->permanent = NULL;
  g_assert (v->i  == GFS_GZ);
#endif /* FTT_3D */
  gfs_centered_variables = gfs_p = v = v->next = 
    gfs_variable_new (gfs_variable_class (), NULL, "P", TRUE, i++);
  g_assert (v->i   == GFS_P);
  v = v->next = gfs_variable_new (gfs_variable_class (), NULL, "U", FALSE, i++);
  g_assert (v->i == GFS_U);
  v = v->next = gfs_variable_new (gfs_variable_class (), NULL, "V", FALSE, i++);
  g_assert (v->i == GFS_V);
#if (!FTT_2D)
  v = v->next = gfs_variable_new (gfs_variable_class (), NULL, "W", FALSE, i++);
  g_assert (v->i == GFS_W);
#endif /* FTT_3D */

  /* Initializes derived variables */
  gfs_derived_first = v = 
    gfs_variable_new (gfs_variable_class (), NULL, "Vorticity", FALSE, GFS_DIV);
  v->derived = cell_vorticity;
  v = v->next = gfs_variable_new (gfs_variable_class (), NULL, "Divergence", FALSE, GFS_DIV);
  v->derived = (GfsVariableDerivedFunc) gfs_divergence;
  v = v->next = gfs_variable_new (gfs_variable_class (), NULL, "Velocity", FALSE, GFS_DIV);
  v->derived = cell_velocity_norm;
  v = v->next = gfs_variable_new (gfs_variable_class (), NULL, "Velocity2", FALSE, GFS_DIV);
  v->derived = cell_velocity_norm2;
  v = v->next = gfs_variable_new (gfs_variable_class (), NULL, "Level", FALSE, GFS_DIV);
  v->derived = cell_level;
  v = v->next = gfs_variable_new (gfs_variable_class (), NULL, "A", FALSE, GFS_DIV);
  v->derived = cell_fraction;
  v = v->next = gfs_variable_new (gfs_variable_class (), NULL, "Lambda2", FALSE, GFS_DIV);
  v->derived = cell_lambda2;
  v = v->next = gfs_variable_new (gfs_variable_class (), NULL, "Curvature", FALSE, GFS_DIV);
  v->derived = cell_curvature;
  gfs_derived_last = v;
  v->next = gfs_centered_variables;

  /* Instantiates classes before reading any domain or simulation file */
  gfs_simulation_class ();
    gfs_ocean_class ();
    gfs_advection_class ();

  gfs_variable_class ();
    gfs_variable_tracer_class ();
    gfs_variable_residual_class ();

  gfs_surface_bc_class ();

  gfs_box_class ();
    gfs_box_not_adapt_class ();
  gfs_gedge_class ();

  gfs_bc_dirichlet_class ();
  gfs_bc_neumann_class ();

  gfs_boundary_class ();
    gfs_boundary_inflow_constant_class ();
    gfs_boundary_outflow_class ();
#ifdef HAVE_MPI
    gfs_boundary_mpi_class ();
#endif /* HAVE_MPI */

  gfs_refine_class ();
    gfs_refine_solid_class ();
    gfs_refine_surface_class ();
      gfs_refine_distance_class ();
      gfs_refine_height_class ();

  gfs_event_class ();
    gfs_init_class ();
    gfs_init_flow_constant_class ();
    gfs_init_fraction_class ();
#if FTT_2D
    gfs_init_vorticity_class ();
#endif /* FTT_2D */
    gfs_adapt_class ();
      gfs_adapt_vorticity_class ();
      gfs_adapt_streamline_curvature_class ();
      gfs_adapt_function_class ();
      gfs_adapt_gradient_class ();
        gfs_adapt_curvature_class ();
      gfs_adapt_not_box_class ();
    gfs_event_sum_class ();
      gfs_event_sum2_class ();
    gfs_event_stop_class ();
    gfs_event_script_class ();
    gfs_source_generic_class ();
      gfs_source_class ();
        gfs_source_control_class ();
      gfs_source_coriolis_class ();
      gfs_source_hydrostatic_class ();
      gfs_source_diffusion_class ();
        gfs_source_diffusion_explicit_class ();
        gfs_source_viscosity_class ();
            gfs_source_vector_class ();
        gfs_source_tension_class ();
    gfs_remove_droplets_class ();
    gfs_remove_ponds_class ();
   
    gfs_output_class ();
      gfs_output_time_class ();
      gfs_output_progress_class ();
      gfs_output_projection_stats_class ();
      gfs_output_diffusion_stats_class ();
      gfs_output_solid_stats_class ();
      gfs_output_adapt_stats_class ();
      gfs_output_timing_class ();
      gfs_output_balance_class ();
      gfs_output_solid_force_class ();
      gfs_output_location_class ();
      gfs_output_simulation_class ();
      gfs_output_boundaries_class ();
      gfs_output_energy_class ();
      gfs_output_particle_class ();

      gfs_output_scalar_class ();
        gfs_output_scalar_norm_class ();
        gfs_output_scalar_stats_class ();
        gfs_output_scalar_sum_class ();
        gfs_output_scalar_histogram_class ();
        gfs_output_error_norm_class ();
          gfs_output_correlation_class ();
	gfs_output_squares_class ();
	gfs_output_streamline_class ();
        gfs_output_ppm_class ();

  /* If modules are not supported, calls modules init functions */
#include "modules.c"
}
