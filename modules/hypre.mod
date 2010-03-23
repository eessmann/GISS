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
#include <_hypre_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>

#include "variable.h"
#include "poisson.h"

/*#define DEBUG*/

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
  gint maxsize;
};

/***********************************************/
/*     Boomer Algebraic Multigrid Solver       */
/***********************************************/
static void call_AMG_Boomer_solver (GfsDomain * domain, GfsMultilevelParams * par,
				    HypreProblem * hp)
{
  HYPRE_Solver solver;
  int num_iterations;
  double final_res_norm;

  gfs_domain_timer_start (domain, "Hypre: AMG_Boomer_solver");

  if (proj_hp.nlevel == 0)
    proj_hp.nlevel = gfs_domain_depth (domain);
  
  /* Create solver */
  HYPRE_BoomerAMGCreate(&solver);

  if (proj_hp.verbose)
    HYPRE_BoomerAMGSetPrintLevel(solver, 3);  /* print solve info + parameters */
  HYPRE_BoomerAMGSetCoarsenType(solver, proj_hp.coarsening_type);
  HYPRE_BoomerAMGSetRelaxType(solver, proj_hp.relax_type);
  HYPRE_BoomerAMGSetCycleType(solver, proj_hp.cycle_type);
  HYPRE_BoomerAMGSetNumSweeps(solver, proj_hp.nrelax);     /* Sweeps on each level */
  HYPRE_BoomerAMGSetMaxLevels(solver, proj_hp.nlevel);  /* maximum number of levels */
  HYPRE_BoomerAMGSetTol(solver, proj_hp.tolerance);        /* conv. tolerance */
  HYPRE_BoomerAMGSetMaxIter(solver, proj_hp.ncyclemax); /* maximum number of iterations */
  HYPRE_BoomerAMGSetMinIter(solver, proj_hp.ncyclemin); /* minimum number of iterations */


  /* Now setup and solve! */
  HYPRE_BoomerAMGSetup(solver, hp->parcsr_A, hp->par_b, hp->par_x);
  HYPRE_BoomerAMGSolve(solver, hp->parcsr_A, hp->par_b, hp->par_x);

  HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
  par->niter = num_iterations;

  /* Prints informations on the residual */
  if (proj_hp.verbose) {
  
    
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

    printf("\n");
    printf("Iterations = %d\n", num_iterations);
    printf("Final Relative Residual Norm = %e\n", final_res_norm);
    printf("\n");
  }
  
  /* Destroy solver */
  HYPRE_BoomerAMGDestroy(solver);
  gfs_domain_timer_stop (domain, "Hypre: AMG_Boomer_solver");
}

/******************************************/
/*       PreConjugateGradient Solver      */
/******************************************/
static void call_PCG_solver (GfsDomain * domain, GfsMultilevelParams * par,
			     HypreProblem * hp)
{
  HYPRE_Solver solver;
  gfs_domain_timer_start (domain, "Hypre: PCG_Solver");

  /* Create solver */
  HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
  
  /* Set some parameters (See Reference Manual for more parameters) */
  if (proj_hp.verbose)
    HYPRE_BoomerAMGSetPrintLevel(solver, 3);  /* print solve info + parameters */
  HYPRE_PCGSetMaxIter(solver, par->nitermax); /* max iterations */
  HYPRE_PCGSetTol(solver, par->tolerance); /* conv. tolerance */
  HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
  HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

  HYPRE_ParCSRPCGSetup(solver, hp->parcsr_A, hp->par_b, hp->par_x);
  HYPRE_ParCSRPCGSolve(solver, hp->parcsr_A, hp->par_b, hp->par_x);

  /*  Run info - needed logging turned on */
  if (proj_hp.verbose) {
    int num_iterations;
    double final_res_norm;
    HYPRE_PCGGetNumIterations(solver, &num_iterations);
    HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    
    printf("\n");
  	  printf("Iterations = %d\n", num_iterations);
  	  printf("Final Relative Residual Norm = %e\n", final_res_norm);
  	  printf("\n");
  }
  
  /* Destroy solver */
  HYPRE_ParCSRPCGDestroy(solver);
  gfs_domain_timer_stop (domain, "Hypre: PCG_Solver");
}



static void hypre_problem_new (HypreProblem * hp, GfsDomain * domain,
			       gdouble size)
{
  gfs_domain_timer_start (domain, "HYPRE: Solver setup");
  
  /* Create the matrix.*/
  HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, size-1, 0, size-1, &hp->A);

  /* Create the vectors rhs and solution.*/
  HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, size-1, &hp->b);
  HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, size-1, &hp->x);
  
  /* Choose a parallel csr format storage */
  HYPRE_IJMatrixSetObjectType(hp->A, HYPRE_PARCSR);
  HYPRE_IJVectorSetObjectType(hp->b, HYPRE_PARCSR);
  HYPRE_IJVectorSetObjectType(hp->x, HYPRE_PARCSR);

  /* Initialize before setting coefficients */
  HYPRE_IJMatrixInitialize(hp->A);
  HYPRE_IJVectorInitialize(hp->b);
  HYPRE_IJVectorInitialize(hp->x);
}

static void hypre_problem_destroy (HypreProblem * hp)
{
  g_assert (hp->A);

  HYPRE_IJMatrixDestroy(hp->A);
  HYPRE_IJVectorDestroy(hp->b);
  HYPRE_IJVectorDestroy(hp->x);
}

static void extract_stencil (GArray * stencil, HypreProblem * hp)
{
  double values[hp->maxsize];
  int cols[hp->maxsize];
  gint i, index;
  
  for (i = 0; i < stencil->len; i++) {
    GfsStencilElement * tmp = &g_array_index (stencil, GfsStencilElement, i);
    if (i == 0)
      index = tmp->cell_id;
    cols[i] = tmp->cell_id;
    values[i] = tmp->cell_coeff;
  }

  HYPRE_IJMatrixSetValues(hp->A, 1, &i, &index, cols, values);
}

static void hypre_problem_init (HypreProblem * hp, GfsLinearProblem * lp,
				GfsDomain * domain)
{
  double *rhs_values, *x_values;
  int    *rows;
  gint i;

  /* Now go through my local rows and set the matrix entries.*/
  rhs_values = malloc(lp->rhs->len * sizeof(double));
  x_values = malloc(lp->lhs->len * sizeof(double));
  rows = malloc(lp->lhs->len * sizeof(int));
    
  hp->maxsize = lp->maxsize;
  g_ptr_array_foreach (lp->LP, (GFunc) extract_stencil, hp);

  for (i = 0; i < lp->rhs->len; i++) {
    rhs_values[i] = g_array_index (lp->rhs, gdouble, i);
    x_values[i] = g_array_index (lp->lhs, gdouble, i);
    rows[i] = i;
  }
  
  HYPRE_IJVectorSetValues(hp->b, lp->rhs->len, rows, rhs_values );
  HYPRE_IJVectorSetValues(hp->x, lp->lhs->len, rows, x_values);

  /* Assemble after setting the coefficients */
  HYPRE_IJMatrixAssemble(hp->A);
  HYPRE_IJVectorAssemble(hp->b);
  HYPRE_IJVectorAssemble(hp->x);

  /* Get the parcsr matrix object to use */
  HYPRE_IJMatrixGetObject(hp->A, (void**) &hp->parcsr_A);
  HYPRE_IJVectorGetObject(hp->b, (void **) &hp->par_b);
  HYPRE_IJVectorGetObject(hp->x, (void **) &hp->par_x);

#ifdef DEBUG
  HYPRE_IJMatrixPrint(hp->A, "Aij.dat");
  HYPRE_IJVectorPrint(hp->x, "xi.dat");
  HYPRE_IJVectorPrint(hp->b, "bi.dat");
#endif

  free(x_values);
  free(rhs_values);
  free(rows);
  gfs_domain_timer_stop (domain, "HYPRE: Solver setup");
}

static void hypre_problem_copy (HypreProblem * hp, GfsLinearProblem * lp)
{
  double *x_values;
  int    *rows;
  gint i;

  /* Copy the solution to the GfsLinearProblem structure */
  x_values = malloc( lp->lhs->len * sizeof(double));
  rows = malloc( lp->lhs->len * sizeof(int));
    
  for (i=0; i< lp->lhs->len; i++) {
    x_values[i] =  0.;
    rows[i] = i;
  }
    
  HYPRE_IJVectorGetValues(hp->x, lp->lhs->len, rows, x_values);
    
  for (i = 0; i < lp->lhs->len; i++)
    g_array_index (lp->lhs, gdouble, i) = x_values[i];

  free(x_values);
  free(rows);
}

static void copy_poisson_solution (FttCell * cell, GfsLinearProblem * lp)
{
  GFS_VALUE (cell, lp->lhs_v) = g_array_index (lp->lhs, gdouble,
					       (gint) GFS_VALUE (cell, lp->id));
}

static void copy_poisson_problem_solution_to_simulation_tree (GfsDomain * domain,
							      GfsLinearProblem * lp)
{
  
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) copy_poisson_solution, lp);
}


static void solve_poisson_problem_using_hypre (GfsDomain * domain,
					       GfsLinearProblem * lp,
					       GfsMultilevelParams * par)
{
  HypreProblem hp;
     
  hypre_problem_new (&hp, domain, lp->rhs->len);
  hypre_problem_init (&hp, lp, domain);
  
  /* Choose a solver and solve the system */
  if (proj_hp.solver_type == HYPRE_BOOMER_AMG)
    call_AMG_Boomer_solver (domain, par, &hp);
  else if (proj_hp.solver_type == HYPRE_PCG)
    call_PCG_solver (domain, par, &hp);
  else
    g_assert_not_reached();
  
  hypre_problem_copy (&hp, lp);
  hypre_problem_destroy (&hp);
}

static void correct (FttCell * cell, gpointer * data)
{
  GfsVariable * u = data[0];
  GfsVariable * dp = data[1];
  GFS_VALUE (cell, u) += GFS_VALUE (cell, dp);
}

/**
 * gfs_hypre_poisson_solve:
 * @domain: the domain over which the poisson problem is solved.
 * @par: the parameters of the poisson problem.
 * @lhs: the variable to use as left-hand side.
 * @rhs: the variable to use as right-hand side.
 * @res: the variable to store the residual
 * @dia: the diagonal weight.
 * @dt: the length of the time-step.
 *
 * Solves the poisson problem over domain using one of the solvers
 * of the HYPRE library.
 *
 * First the poisson problem is extracted as a GfsLinearProblem, which
 * is then fed to the HYPRE library.
 *
 * The solution is then copied to the quadtree.
 */
static void gfs_hypre_poisson_solve (GfsDomain * domain,
				     GfsMultilevelParams * par,
				     GfsVariable * lhs,
				     GfsVariable * rhs,
				     GfsVariable * res,
				     GfsVariable * dia,
				     gdouble dt)
{
  GfsVariable * dp = gfs_temporary_variable (domain);
  gpointer data[2];

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, dp);

  /* calculates the intial residual and its norm */
  gfs_residual (domain, par->dimension, FTT_TRAVERSE_LEAFS, -1, lhs, rhs, dia, res);
  par->residual_before = par->residual = 
    gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, dt, res);

  GfsLinearProblem * lp = gfs_get_poisson_problem (domain, par, res, dp, dia, -1);
 
  solve_poisson_problem_using_hypre (domain, lp, par);

  copy_poisson_problem_solution_to_simulation_tree (domain, lp);

  gfs_linear_problem_destroy (lp);

  /* correct on leaf cells */
  data[0] = lhs;
  data[1] = dp;
  gfs_traverse_and_bc (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) correct, data,
		       lhs, lhs);
  /* compute new residual on leaf cells */
  gfs_residual (domain, par->dimension, FTT_TRAVERSE_LEAFS, -1, lhs, rhs, dia, res);
  par->residual = gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, dt, res);
  gts_object_destroy (GTS_OBJECT (dp));
}

static void hypre_solver_write (HypreSolverParams * par,FILE * fp)
{
 
  if ( par == NULL) {
    fputs ("{\n",fp);
    fputs ("  solver_type      = boomer_amg\n", fp);
    fputs ("  relax_type       = sor-j-forward\n", fp);
    fputs ("  precond_type     = none\n", fp);
    fputs ("  coarsening_type  = cgc_e\n", fp);
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
  par->nlevel = 0;
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
      printf("Warning *** Coarsening algorithms are only for the BoomerAMG Solver !!\n");
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

  if ( par->cycle_type < 1 || (par->cycle_type > 8 && par->cycle_type < 11) ||
       par->cycle_type > 14)
    gts_file_variable_error (fp, var, "cycle_type",
			     "unknown Cycle Type `%i'", par->cycle_type);
  
  if (par->nlevel < 0)
    gts_file_variable_error (fp, var, "nlevel",
			     "error in hypre solver parameter nlevel < 0.");

  if (par->verbose != 1 && par->verbose != 0)
    gts_file_variable_error (fp, var, "verbose",
			     "error in hypre solver parameter verbose != 0 or 1: `%i'",
			     par->verbose);
  
  if (par->ncyclemax < 1)
    gts_file_variable_error (fp, var, "verbose",
			     "error in hypre solver parameter ncyclemax can't be < 1  `%i'",
			     par->ncyclemax);

  if (par->ncyclemin < 1)
    gts_file_variable_error (fp, var, "verbose",
			     "error in hypre solver parameter ncyclemin can't be < 1  `%i'",
			     par->ncyclemin);

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

/* Initialize module */
void          gfs_module_read     (GtsFile * fp, GfsSimulation * sim);
void          gfs_module_write    (FILE * fp);

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  return NULL;
}

void gfs_module_read (GtsFile * fp, GfsSimulation * sim)
{
  g_return_if_fail (fp != NULL);

  hypre_solver_read (&proj_hp, fp);
  
  /* initialise the poisson cycle hook */
  sim->approx_projection_params.poisson_solve = gfs_hypre_poisson_solve;
  sim->projection_params.poisson_solve = gfs_hypre_poisson_solve;
}

void gfs_module_write (FILE * fp)
{
  g_return_if_fail (fp != NULL);

  hypre_solver_write (&proj_hp, fp);
}
