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
#include "mpi_boundary.h"

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
};

/* Parameters for the projection schemes are stored in proj_hp */
HypreSolverParams proj_hp;

struct _HypreProblem {
  HYPRE_IJMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_IJVector b;
  HYPRE_ParVector par_b;
  HYPRE_IJVector x;
  HYPRE_ParVector par_x;
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
  HYPRE_BoomerAMGSetNumSweeps(solver, par->nrelax);     /* Sweeps on each level */
  HYPRE_BoomerAMGSetMaxLevels(solver, proj_hp.nlevel);  /* maximum number of levels */
  HYPRE_BoomerAMGSetTol(solver, par->tolerance);        /* conv. tolerance */
  HYPRE_BoomerAMGSetMaxIter(solver, par->nitermax); /* maximum number of iterations */
  HYPRE_BoomerAMGSetMinIter(solver, par->nitermin); /* minimum number of iterations */

  /* Now setup and solve! */
  HYPRE_BoomerAMGSetup(solver, hp->parcsr_A, hp->par_b, hp->par_x);
  HYPRE_BoomerAMGSolve(solver, hp->parcsr_A, hp->par_b, hp->par_x);

  HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
  par->niter = num_iterations;

  /* Prints informations on the residual */
  if (proj_hp.verbose && domain->pid <= 0) {      
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
  if (proj_hp.verbose && domain->pid <= 0) {
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
			       gint size, gint istart)
{
  gfs_domain_timer_start (domain, "HYPRE: Solver setup");
  
  /* Create the matrix.*/
  HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 
		       istart, istart + size - 1, istart, istart + size - 1, 
		       &hp->A);

  /* Create the vectors rhs and solution.*/
  HYPRE_IJVectorCreate(MPI_COMM_WORLD, istart, istart + size-1, &hp->b);
  HYPRE_IJVectorCreate(MPI_COMM_WORLD, istart, istart + size-1, &hp->x);
  
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

static void extract_stencil (GfsStencil * stencil, HypreProblem * hp)
{
  int ncols = stencil->id->len;
  int rows = g_array_index (stencil->id, int, 0);
  HYPRE_IJMatrixSetValues (hp->A, 1, &ncols, &rows, 
			   (int *) stencil->id->data, (double *) stencil->coeff->data);
}

static void hypre_problem_init (HypreProblem * hp, GfsLinearProblem * lp,
				GfsDomain * domain)
{
  double *rhs_values, *x_values;
  int    *rows;
  gint i;

  /* Now go through my local rows and set the matrix entries.*/
  rhs_values = (double *) lp->rhs->data;
  x_values = (double *) lp->lhs->data;
  rows = g_malloc (lp->lhs->len*sizeof(int));
  g_ptr_array_foreach (lp->LP, (GFunc) extract_stencil, hp);

  for (i = 0; i < lp->rhs->len; i++)
    rows[i] = lp->istart + i;
  
  HYPRE_IJVectorSetValues(hp->b, lp->rhs->len, rows, rhs_values );
  HYPRE_IJVectorSetValues(hp->x, lp->lhs->len, rows, x_values);

  /* Assemble after setting the coefficients */
  HYPRE_IJMatrixAssemble(hp->A);
  HYPRE_IJVectorAssemble(hp->b);
  HYPRE_IJVectorAssemble(hp->x);

  /* Get the parcsr matrix object to use */
  HYPRE_IJMatrixGetObject(hp->A, (void **) &hp->parcsr_A);
  HYPRE_IJVectorGetObject(hp->b, (void **) &hp->par_b);
  HYPRE_IJVectorGetObject(hp->x, (void **) &hp->par_x);

#ifdef DEBUG
  HYPRE_IJMatrixPrint(hp->A, "Aij.dat");
  HYPRE_IJVectorPrint(hp->x, "xi.dat");
  HYPRE_IJVectorPrint(hp->b, "bi.dat");
#endif

  free(rows);
  gfs_domain_timer_stop (domain, "HYPRE: Solver setup");
}

static void hypre_problem_copy (HypreProblem * hp, GfsLinearProblem * lp)
{
  double *x_values;
  int    *rows;
  gint i;

  /* Copy the solution to the GfsLinearProblem structure */
  x_values = g_malloc (lp->lhs->len*sizeof (double));
  rows = g_malloc (lp->lhs->len*sizeof (int));
    
  for (i = 0; i < lp->lhs->len; i++) {
    x_values[i] = 0.;
    rows[i] = i + lp->istart;
  }
    
  HYPRE_IJVectorGetValues(hp->x, lp->lhs->len, rows, x_values);
    
  for (i = 0; i < lp->lhs->len; i++)
    g_array_index (lp->lhs, gdouble, i) = x_values[i];

  free(x_values);
  free(rows);
}

typedef struct {
  GfsLinearProblem * lp;
  GfsVariable * lhs;
} CopyParams;

static void copy_poisson_solution (FttCell * cell, CopyParams * p)
{
   GFS_VALUE (cell, p->lhs) =
     g_array_index (p->lp->lhs, gdouble, (int) GFS_VALUE (cell, p->lp->id) - p->lp->istart);
}

static void bc_copy_solution (FttCellFace * f, CopyParams * p)
{
  GFS_VALUE (f->cell, p->lhs) =
    g_array_index (p->lp->lhs, gdouble, (int) GFS_VALUE (f->cell, p->lp->id) - p->lp->istart);
}

static void box_copy_poisson_solution_bc (GfsBox * box, CopyParams * p)
{ 
  FttDirection d;
  
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      if(!GFS_IS_BOUNDARY_MPI(box->neighbor[d])) {
	GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
	ftt_face_traverse_boundary (b->root, b->d,
				    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				    (FttFaceTraverseFunc) bc_copy_solution, p);
      }
    }
}

static void solve_poisson_problem_using_hypre (GfsDomain * domain,
					       GfsLinearProblem * lp,
					       GfsMultilevelParams * par)
{
  HypreProblem hp;
  gdouble tolerance = par->tolerance, res0 = 0.;
  gdouble size = ftt_level_size (gfs_domain_depth (domain));
  gint i, len = lp->rhs->len;

  for (i=0;i< lp->rhs->len;i++)
    res0 += pow(g_array_index (lp->rhs, gdouble, i),2.);
  res0 = sqrt(res0);

  if (domain->pid >= 0) {
    gfs_all_reduce (domain, res0, MPI_DOUBLE, MPI_SUM);
    gfs_all_reduce (domain, len, MPI_INT, MPI_SUM);
  }

  /* Tolerance has to be rescaled to account for the different of method */
  /* used by Hypre to computed the norm of the residual */
  par->tolerance *= sqrt(((gdouble) len))/res0;
 
  hypre_problem_new (&hp, domain, lp->rhs->len, lp->istart);
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
  par->tolerance = tolerance;
}

static void correct (FttCell * cell, gpointer * data)
{
  GfsVariable * u = data[0];
  GfsVariable * dp = data[1];
  GFS_VALUE (cell, u) += GFS_VALUE (cell, dp);
}

static void correct_bc (FttCellFace * f, gpointer * data)
{
  GfsVariable * u = data[0];
  GfsVariable * dp = data[1];

  GFS_VALUE (f->cell, u) += GFS_VALUE (f->cell, dp);
}

static void box_correct_bc (GfsBox * box, gpointer data)
{ 
  FttDirection d;
  
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      if(!GFS_IS_BOUNDARY_MPI(box->neighbor[d])) {
	GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
	ftt_face_traverse_boundary (b->root, b->d,
				    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				    (FttFaceTraverseFunc) correct_bc, data);
      }
    }
}

static void hypre_poisson_solve (GfsDomain * domain,
				 GfsMultilevelParams * par,
				 GfsVariable * lhs,
				 GfsVariable * rhs,
				 GfsVariable * res,
				 GfsVariable * dia,
				 gdouble dt)
{
  /* calculates the initial residual and its norm */
  gfs_residual (domain, par->dimension, FTT_TRAVERSE_LEAFS, -1, lhs, rhs, dia, res);
  par->residual_before = par->residual = 
    gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, dt, res);

  if (par->nitermax > 0) {
    GfsVariable * dp = gfs_temporary_variable (domain);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, dp);
    GfsLinearProblem * lp = gfs_get_poisson_problem (domain, res, dp, dia, -1, lhs);
 
    solve_poisson_problem_using_hypre (domain, lp, par);

    CopyParams p = { lp, dp };
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
    			      (FttCellTraverseFunc) copy_poisson_solution, &p);
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_copy_poisson_solution_bc, &p);

    gfs_linear_problem_destroy (lp);

    /* correct on leaf cells */
    gpointer data[2];
    data[0] = lhs;
    data[1] = dp;

    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) correct, data);
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_correct_bc, data);

    if ( domain->pid >= 0 )
      gfs_domain_copy_bc (domain, FTT_TRAVERSE_LEAFS, -1, lhs, lhs);

    /* compute new residual on leaf cells */
    gfs_residual (domain, par->dimension, FTT_TRAVERSE_LEAFS, -1, lhs, rhs, dia, res);
    par->residual = gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, dt, res);
    gts_object_destroy (GTS_OBJECT (dp));
  }
}

static void hypre_solver_write (HypreSolverParams * par,FILE * fp)
{
  fprintf (fp," {\n");
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

  fputc ('}', fp);
}

static void hypre_solver_read (HypreSolverParams * par, GtsFile * fp)
{
  gchar * solver_type = NULL, * relax_type = NULL, * coarsening_type = NULL;
  gchar * precond_type = NULL;
  GtsFileVariable var[] = {
    {GTS_STRING,  "relax_type",      TRUE, &relax_type},
    {GTS_STRING,  "solver_type",     TRUE, &solver_type},
    {GTS_STRING,  "precond_type",    TRUE, &precond_type},
    {GTS_STRING,  "coarsening_type", TRUE, &coarsening_type},
    {GTS_INT,     "cycle_type",      TRUE, &par->cycle_type},
    {GTS_INT,     "nlevel",          TRUE, &par->nlevel},
    {GTS_INT,     "verbose",         TRUE, &par->verbose},
    {GTS_NONE}
  };

  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (solver_type) {
    if (!strcmp (solver_type, "boomer_amg"))
      par->solver_type = HYPRE_BOOMER_AMG;
    else if (!strcmp (solver_type, "pcg"))
      par->solver_type = HYPRE_PCG;
    else
      gts_file_variable_error (fp, var, "solver_type", "unknown solver type `%s'", solver_type);
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
    else
      gts_file_variable_error (fp, var, "precond_type", 
			       "unknown preconditioner `%s'", precond_type);
    g_free (precond_type);
  }
  if (fp->type == GTS_ERROR)
    return;

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
    else
      gts_file_variable_error (fp, var, "relax_type", "unknown relax type `%s'", relax_type);
    g_free (relax_type);
  }
  if (fp->type == GTS_ERROR)
    return;

  if (coarsening_type) {
    if (par->solver_type != HYPRE_BOOMER_AMG)
      g_warning ("coarsening algorithms are only for the BoomerAMG Solver !!\n"
		 "none will be used with the selected solver");

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
    else
      gts_file_variable_error (fp, var, "coarsening_type",
			       "unknown coarsening type `%s'", coarsening_type);
    g_free (coarsening_type);
  }
  if (fp->type == GTS_ERROR)
    return;

  if (par->cycle_type < 1 || (par->cycle_type > 8 && par->cycle_type < 11) ||
      par->cycle_type > 14)
    gts_file_variable_error (fp, var, "cycle_type",
			     "unknown cycle type `%i'", par->cycle_type);
  else if (par->nlevel < 0)
    gts_file_variable_error (fp, var, "nlevel", "nlevel cannot be < 0");
}

/* Initialize module */
void          gfs_module_read     (GtsFile * fp, GfsSimulation * sim);
void          gfs_module_write    (FILE * fp);

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "hypre";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  return NULL;
}

void gfs_module_read (GtsFile * fp, GfsSimulation * sim)
{
  g_return_if_fail (fp != NULL);

  /* Default hypre solver parameters */
  /* Boomer AMG is the default solver */
  proj_hp.solver_type = HYPRE_BOOMER_AMG;
  proj_hp.precond_type = NO_PRECOND;
  proj_hp.relax_type = 5;
  proj_hp.coarsening_type = 22;
  proj_hp.cycle_type = 1;
  proj_hp.nlevel = 0;
  proj_hp.verbose = FALSE;

  if (fp->type == '{')
    hypre_solver_read (&proj_hp, fp);
  
  /* initialise the poisson cycle hook */
  sim->approx_projection_params.poisson_solve = hypre_poisson_solve;
  sim->projection_params.poisson_solve = hypre_poisson_solve;
}

void gfs_module_write (FILE * fp)
{
  g_return_if_fail (fp != NULL);

  hypre_solver_write (&proj_hp, fp);
}
