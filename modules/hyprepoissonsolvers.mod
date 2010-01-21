/* Gerris - The GNU Flow Solver                       (-*-C-*-)
 * Copyright (C) 2001-2008 National Institute of Water and Atmospheric Research
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

#include <math.h>
#include <stdlib.h>
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "gfs.h"

typedef enum {
  HYPRE_BOOMER_AMG,
  HYPRE_PCG
} HypreSolverType;

typedef enum {
  HYPRE_AMG_PRECOND,
  HYPRE_PARASAILS_PRECOND,
  NO_PRECOND
} HyprePrecondType;

typedef struct _HypreProblem HypreProblem;
typedef struct _HypreSolverParams HypreSolverParams;

struct _HypreSolverParams {
  HypreSolverType solver_type;
  HyprePrecondType precond_type;
  gint relax_type;
  gint coarsening_type;
  gint cycle_type;
  gint nlevel;
  gboolean verbose;
  gint ncyclemax;
  gint ncyclemin;
  gint niter;
  gint nrelax;
  gdouble tolerance;
};

/* Parameters to the projection schemes are stored in proj_hp */
HypreSolverParams proj_hp;

struct _HypreProblem {
  HYPRE_IJMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_IJVector b;
  HYPRE_ParVector par_b;
  HYPRE_IJVector x;
  HYPRE_ParVector par_x;
};

static void hypre_solver_write (HypreSolverParams * par,FILE * fp)
{
 
  if ( par == NULL) {   
    fputs ("{\n",fp);
    fputs ("  solver_type      = boomer_amg\n", fp);
    fputs ("  relax_type       = sor-j-forward\n", fp);
    fputs ("  precond_type     = none\n", fp);
    fputs ("  coarsening_type  = hmis\n", fp);
    fputs ("  cycle_type       = 1\n", fp);
    fputs ("  nlevel           = 11\n", fp);
    fputs ("}\n", fp);
  }
  else {
    g_return_if_fail (par != NULL);
    g_return_if_fail (fp != NULL);
    
    fprintf (fp,"{\n");
    switch (par->solver_type) {
    case HYPRE_BOOMER_AMG: fputs ("  solver_type      = boomer_amg\n", fp); break;
    case HYPRE_PCG:        fputs ("  solver_type      = pcg\n", fp); break;
    }
    
    switch (par->relax_type) {
    case 0: fputs ("  relax_type       = jacobi\n", fp); break;
    case 1: fputs ("  relax_type       = gauss_seidel\n", fp); break;
    case 3: fputs ("  relax_type       = sor-j-forward\n", fp); break;
    case 4: fputs ("  relax_type       = sor-j-backward\n", fp); break;
    case 5: fputs ("  relax_type       = gs-j\n", fp); break;
    case 6: fputs ("  relax_type       = ssor-j\n", fp); break;
    case 7: fputs ("  relax_type       = matvec-jacobi\n", fp); break;
    case 9: fputs ("  relax_type       = direct\n", fp); break;
    }

    switch (par->precond_type) {
    case HYPRE_AMG_PRECOND: fputs       ("  precond_type     = amg\n", fp); break;
    case HYPRE_PARASAILS_PRECOND: fputs ("  precond_type     = parasails\n", fp); break;
    case NO_PRECOND: fputs              ("  precond_type     = none\n", fp); break; 
    }
    
    switch (par->coarsening_type) {
    case 0: fputs  ("  coarsening_type  = cljp\n", fp); break;
    case 3: fputs  ("  coarsening_type  = ruge_stueben\n", fp); break;
    case 6: fputs  ("  coarsening_type  = falgout\n", fp); break;
    case 8: fputs  ("  coarsening_type  = pmis\n", fp); break;
    case 10: fputs ("  coarsening_type  = hmis\n", fp); break;
    case 21: fputs ("  coarsening_type  = cgc\n", fp); break;
    case 22: fputs ("  coarsening_type  = cgc_e\n", fp); break;
    }

    fprintf (fp,"  cycle_type       = %i\n", par->cycle_type);
    
    fprintf (fp,"  nlevel           = %i\n", par->nlevel);
    
    fprintf (fp,"  verbose          = %i\n", par->verbose);

    fprintf (fp,"  ncyclemax        = %i\n", par->ncyclemax);

    fprintf (fp,"  ncyclemin        = %i\n", par->ncyclemin);

    fprintf (fp,"  nrelax           = %i\n", par->nrelax);

    fprintf (fp,"  tolerance        = %e\n", par->tolerance);

    fputc ('}', fp);
  }
}

static void hypre_solver_read (HypreSolverParams * par, GtsFile * fp)
{
  GtsFileVariable var[] = {
    {GTS_STRING,  "relax_type",      TRUE},
    {GTS_STRING,  "solver_type",     TRUE},
    {GTS_STRING,  "precond_type",    TRUE},
    {GTS_STRING,  "coarsening_type", TRUE},
    {GTS_INT,     "cycle_type",      TRUE},
    {GTS_INT,     "nlevel",          TRUE},
    {GTS_INT,     "verbose",         TRUE},
    {GTS_INT,     "ncyclemax",       TRUE},
    {GTS_INT,     "ncyclemin",       TRUE},
    {GTS_DOUBLE,  "tolerance",       TRUE},
    {GTS_INT,     "nrelax",          TRUE},
    {GTS_NONE}
  };
  gchar * solver_type = NULL, * relax_type = NULL, * coarsening_type = NULL;
  gchar * precond_type = NULL;

  /* Default hypre solver parameters */
  /* Boomer AMG is the default solver */
  par->solver_type = HYPRE_BOOMER_AMG;
  par->precond_type = NO_PRECOND;
  par->relax_type = 5;
  par->coarsening_type = 22;
  par->cycle_type = 1;
  par->nlevel = 11;
  par->verbose = FALSE;
  par->ncyclemax = 100;
  par->ncyclemin = 1;
  par->niter = 0;
  par->nrelax = 4;
  par->tolerance = 1e-5;

  g_assert (par != NULL);
  g_assert (fp != NULL);

  var[0].data = &relax_type;
  var[1].data = &solver_type;
  var[2].data = &precond_type;
  var[3].data = &coarsening_type;
  var[4].data = &par->cycle_type;
  var[5].data = &par->nlevel;
  var[6].data = &par->verbose;
  var[7].data = &par->ncyclemax;
  var[8].data = &par->ncyclemin;
  var[9].data = &par->tolerance;
  var[10].data= &par->nrelax;

  gts_file_assign_variables (fp, var);

  if (solver_type) {
    if (!strcmp (solver_type, "boomer_amg"))
      par->solver_type = HYPRE_BOOMER_AMG;
    else if (!strcmp (solver_type, "pcg"))
      par->solver_type = HYPRE_PCG;
    else if (fp->type != GTS_ERROR)
      gts_file_variable_error (fp, var, "solver_type",
			       "unknown Hypre Solver `%s'", solver_type);
    g_free (solver_type);
  }

  /* The default is no preconditioner */
  if (precond_type) {
    if (!strcmp (precond_type, "amg"))
      par->precond_type = HYPRE_AMG_PRECOND;
    else if (!strcmp (precond_type, "parasails"))
      par->precond_type = HYPRE_PARASAILS_PRECOND;
    else if (!strcmp (precond_type, "none"))
      par->precond_type = NO_PRECOND;
    else if (fp->type != GTS_ERROR)
      gts_file_variable_error (fp, var, "precond_type",
			       "unknown Hypre Preconditioner `%s'", precond_type);
    g_free (precond_type);
  }

  if (relax_type) {
    if (!strcmp (relax_type, "jacobi"))
      par->relax_type = 0;
    else if (!strcmp (relax_type, "gauss_seidel"))
      par->relax_type = 1;
    else if (!strcmp (relax_type, "sor-j-forward"))
      par->relax_type = 3;
    else if (!strcmp (relax_type, "sor-j-backward"))
      par->relax_type = 4;
    else if (!strcmp (relax_type, "gs-j"))
      par->relax_type = 5;
    else if (!strcmp (relax_type, "ssor-j"))
      par->relax_type = 6;
    else if (!strcmp (relax_type, "matvec-jacobi"))
      par->relax_type = 7;
    else if (!strcmp (relax_type, "direct"))
      par->relax_type = 9;
    else if (fp->type != GTS_ERROR)
      gts_file_variable_error (fp, var, "relax_type",
			       "unknown Hypre Relax `%s'", relax_type);
    g_free (relax_type);
  }

  if (coarsening_type) {
    if (par->solver_type != HYPRE_BOOMER_AMG) {
      printf("Warning *** Coarsening algorithms are only for the Boomer AMG Solver !!\n");
      printf("Warning *** None will be used with the selected solver.\n");
    }

    if (!strcmp (coarsening_type, "cljp"))
      par->coarsening_type = 0;
    else if (!strcmp (coarsening_type, "ruge_stueben"))
      par->coarsening_type = 3;
    else if (!strcmp (coarsening_type, "falgout"))
      par->coarsening_type = 6;
    else if (!strcmp (coarsening_type, "pmis"))
      par->coarsening_type = 8;
    else if (!strcmp (coarsening_type, "hmis"))
      par->coarsening_type = 10;
    else if (!strcmp (coarsening_type, "cgc"))
      par->coarsening_type = 21;
    else if (!strcmp (coarsening_type, "cgc_e"))
      par->coarsening_type = 22;
    else if (fp->type != GTS_ERROR)
      gts_file_variable_error (fp, var, "coarsening_type",
			       "unknown Hypre Coarsening `%s'", coarsening_type);
    g_free (coarsening_type);
  }

  if (par->cycle_type != 1 && par->cycle_type != 2 && par->cycle_type != 3 && par->cycle_type != 4 && par->cycle_type != 5 && par->cycle_type != 6 && par->cycle_type != 7 && par->cycle_type != 8 && par->cycle_type != 11 && par->cycle_type != 12 && par->cycle_type != 13 && par->cycle_type != 14)
    gts_file_variable_error (fp, var, "cycle_type",
			       "unknown Cycle Type `%i'", par->cycle_type);
    
  if (par->nlevel < 1)
    gts_file_variable_error (fp, var, "nlevel",
			       "error in hypre solver parameter nlevel < 0.");

  if (par->verbose != 1 && par->verbose != 0)
    gts_file_variable_error (fp, var, "verbose",
			     "error in hypre solver parameter verbose != 0 or 1: `%i'", par->verbose);

  if (par->ncyclemax < 1)
    gts_file_variable_error (fp, var, "verbose",
			     "error in hypre solver parameter ncyclemax can't be < 1  `%i'", par->ncyclemax);

  if (par->ncyclemin < 1)
    gts_file_variable_error (fp, var, "verbose",
			     "error in hypre solver parameter ncyclemin can't be < 1  `%i'", par->ncyclemin);

  if (par->tolerance <= 0.) {
    gts_file_variable_error (fp, var, "tolerance",
			     "tolerance `%g' must be strictly positive",
			     par->tolerance);
    return;
  }

  if (par->nrelax <= 0.) {
    gts_file_variable_error (fp, var, "nrelax",
			     "nrelax `%i' must be strictly positive",
			     par->nrelax);
    return;
  }
}

static void hypre_projection_params_read (GfsMultilevelParams * par, GtsFile * fp)
{
  hypre_solver_read (&proj_hp, fp);
}

static void hypre_projection_params_write (GfsMultilevelParams * par, FILE * fp)
{
  hypre_solver_write (&proj_hp, fp);
}

static void hypre_approx_params_read (GfsMultilevelParams * par, GtsFile * fp)
{
  printf("No distinction is made between approx and mac projection\n");
  printf(" in the hyprepoissonsolvers module \n");
  printf("ApproxProjectionParams won't get read \n");
  printf("Use ProjectionParams to specify parameters");
}

/* Initialize module */
void          gfs_module_read     (GtsFile * fp, GfsSimulation * sim);

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{

  return NULL;
}

void gfs_module_read (GtsFile * fp, GfsSimulation * sim)
{
  sim->approx_projection_params.read = hypre_approx_params_read;
  sim->approx_projection_params.write =  hypre_projection_params_write;

  sim->projection_params.read = hypre_projection_params_read;
  sim->projection_params.write = hypre_projection_params_write;
}
