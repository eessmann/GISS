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

#include <math.h>
#include "poisson.h"
#include "solid.h"
#include "source.h"
#include "tension.h"

/**
 * gfs_multilevel_params_write:
 * @par: the multilevel parameters.
 * @fp: a file pointer.
 *
 * Writes in @fp a text representation of the multilevel parameters
 * @par.  
 */
void gfs_multilevel_params_write (GfsMultilevelParams * par, FILE * fp)
{
  g_return_if_fail (par != NULL);
  g_return_if_fail (fp != NULL);

  fprintf (fp,
           "{\n"
	   "  tolerance = %g\n"
	   "  nrelax    = %u\n"
           "  erelax    = %u\n"
	   "  minlevel  = %u\n"
	   "  nitermax  = %u\n"
	   "  nitermin  = %u\n"
	   "  weighted  = %d\n"
	   "  beta      = %g\n",
	   par->tolerance,
	   par->nrelax,
	   par->erelax,
	   par->minlevel,
	   par->nitermax,
	   par->nitermin,
	   par->weighted,
	   par->beta);
  if (par->omega != 1.)
    fprintf (fp, "  omega     = %g\n", par->omega);
  fputc ('}', fp);
}

void gfs_multilevel_params_init (GfsMultilevelParams * par)
{
  g_return_if_fail (par != NULL);

  par->tolerance = 1e-3;
  par->nrelax    = 4;
  par->erelax    = 1;
  par->minlevel  = 0;
  par->nitermax  = 100;
  par->nitermin  = 1;

  par->dimension = FTT_DIMENSION;
  par->weighted = FALSE;
  par->beta = 0.5;
  par->omega = 1.;
}

void gfs_multilevel_params_read (GfsMultilevelParams * par, GtsFile * fp)
{
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "tolerance", TRUE},
    {GTS_UINT,   "nrelax",    TRUE},
    {GTS_UINT,   "erelax",    TRUE},
    {GTS_UINT,   "minlevel",  TRUE},
    {GTS_UINT,   "nitermax",  TRUE},
    {GTS_UINT,   "nitermin",  TRUE},
    {GTS_INT,    "weighted",  TRUE},
    {GTS_DOUBLE, "beta",      TRUE},
    {GTS_DOUBLE, "omega",     TRUE},
    {GTS_NONE}
  };

  g_return_if_fail (par != NULL);
  g_return_if_fail (fp != NULL);

  var[0].data = &par->tolerance;
  var[1].data = &par->nrelax;
  var[2].data = &par->erelax;
  var[3].data = &par->minlevel;
  var[4].data = &par->nitermax;
  var[5].data = &par->nitermin;
  var[6].data = &par->weighted;
  var[7].data = &par->beta;
  var[8].data = &par->omega;

  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (par->tolerance <= 0.) {
    gts_file_variable_error (fp, var, "tolerance",
			     "tolerance `%g' must be strictly positive",
			     par->tolerance);
    return;
  }
  if (par->nrelax == 0)
    gts_file_variable_error (fp, var, "nrelax", "nrelax must be non zero");
  if (par->erelax == 0)
    gts_file_variable_error (fp, var, "erelax", "erelax must be non zero");
  if (par->beta < 0.5 || par->beta > 1.)
    gts_file_variable_error (fp, var, "beta", "beta must be in [0.5,1]");
}

static gdouble rate (gdouble a, gdouble b, guint n)
{
  if (a > 0. && b > 0. && n > 0)
    return exp (log (b/a)/n);
  return 0.;
}

/**
 * gfs_multilevel_params_stats_write:
 * @par: the multilevel parameters.
 * @fp: a file pointer.
 *
 * Writes in @fp the statistics contained in @p.
 */
void gfs_multilevel_params_stats_write (GfsMultilevelParams * par,
					FILE * fp)
{
  g_return_if_fail (par != NULL);
  g_return_if_fail (fp != NULL);

  fprintf (fp,
	   "    niter: %4d\n"
	   "    residual.bias:   % 10.3e % 10.3e\n"
	   "    residual.first:  % 10.3e % 10.3e %6.2g\n"
	   "    residual.second: % 10.3e % 10.3e %6.2g\n"
	   "    residual.infty:  % 10.3e % 10.3e %6.2g\n",
	   par->niter,
	   par->residual_before.bias,
	   par->residual.bias,
	   par->residual_before.first,
	   par->residual.first,
	   rate (par->residual.first,
		 par->residual_before.first,
		 par->niter),
	   par->residual_before.second,
	   par->residual.second,
	   rate (par->residual.second,
		 par->residual_before.second,
		 par->niter),
	   par->residual_before.infty,
	   par->residual.infty,
	   rate (par->residual.infty,
		 par->residual_before.infty,
		 par->niter));
}

      /** Methods for GfsLinearProblem **/
/**
 * gfs_linear_problem_new:
 *
 * Creates a new GfsLinearProblem.
 * Returns a pointer on the new GfsLinearProblem. 
 */
GfsLinearProblem * gfs_linear_problem_new ()
{
  return g_malloc (sizeof (GfsLinearProblem));
}

/**
 * gfs_linear_problem_init:
 * @lp: a pointer on a GfsLinearProblem.
 * 
 * Initialises a GfsLinearProblem. Creates the structure
 * to store the problem LP, the right hand and left hand
 * side vectors rhs and lhs.
 * Initialises the maximum size for a stencil to 0.
 * And the number of diagonal element to 0.
 */
void gfs_linear_problem_init (GfsLinearProblem * lp)
{
  lp->rhs = g_array_new (FALSE, FALSE, sizeof (gdouble));
  lp->lhs = g_array_new (FALSE, FALSE, sizeof (gdouble));
  lp->LP = g_ptr_array_new ();
  lp->maxsize = 0;
  lp->nleafs = 0;
}

/**
 * gfs_linear_problem_add_stencil:
 * @lp: a pointer on a GfsLinearProblem.
 * @stencil: a pointer on a GfsStencil. 
 *
 * Adds a stencil to the linear problem.
 * If the stencil is larger than the previous ones
 * lp->maxsize is updated.
 */
void gfs_linear_problem_add_stencil (GfsLinearProblem * lp, GfsStencil * stencil)
{
  g_assert (stencil != NULL);

  g_ptr_array_add (lp->LP, stencil);
  
  if (stencil->data->len > lp->maxsize)
    lp->maxsize = stencil->data->len;
}

static void destroy_stencil (GfsStencil * stencil)
{
  gfs_stencil_destroy (stencil);
}

/**
 * gfs_linear_problem_destroy:
 * @lp: a pointer on a GfsLinearProblem.
 * 
 * Destroys a GfsLinearProblem.
 */
void gfs_linear_problem_destroy (GfsLinearProblem * lp)
{
  gts_object_destroy (GTS_OBJECT (lp->id));

  g_array_free (lp->rhs, TRUE);  
  g_array_free (lp->lhs, TRUE);
  
  g_ptr_array_foreach (lp->LP, (GFunc) destroy_stencil, NULL);
  g_ptr_array_free (lp->LP, TRUE);
}

/*******************************************************************/

static void relax_coeff_stencil (FttCell * cell, GfsLinearProblem * lp)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;
  
  GfsStencil * stencil = gfs_stencil_new ();
  stencil->id = lp->id;
  stencil->u = lp->u;

  gfs_stencil_add_element (stencil, (gint) GFS_VALUE (cell, lp->id), 0.);

  g.a = GFS_VALUE (cell, lp->dia);
  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = neighbor.c[f.d];
    if (f.neighbor) {
      gfs_face_weighted_gradient_stencil (&f, &ng, lp->maxlevel, stencil);
      g.a += ng.a;
    }
  }

  if (g.a > 0.)
    gfs_stencil_add_element (stencil, (gint) GFS_VALUE (cell, lp->id),  -g.a);
  else {
    gfs_stencil_reinit (stencil);
    gfs_stencil_add_element (stencil, (gint) GFS_VALUE (cell, lp->id), 1.);
    g_array_index (lp->rhs, gdouble, (gint) GFS_VALUE (cell, lp->id)) = 0.;
  }

  gfs_linear_problem_add_stencil (lp, stencil);
}

static void leafs_numbering (FttCell * cell, GfsLinearProblem * lp) {
 
  GFS_VALUE (cell, lp->id) = (gdouble) lp->nleafs;
  g_array_append_val (lp->lhs, GFS_VALUE (cell, lp->lhs_v));
  g_array_append_val (lp->rhs, GFS_VALUE (cell, lp->rhs_v));
  lp->nleafs++;
}

static void bc_number (FttCellFace * f, GfsLinearProblem * lp)
{
  GFS_VALUE(f->cell, lp->id) = (gdouble) lp->nleafs;
  g_array_append_val (lp->lhs, GFS_VALUE (f->cell, lp->lhs_v));
  g_array_append_val (lp->rhs, GFS_VALUE (f->cell, lp->rhs_v));
  lp->nleafs++;
}

static void bc_leafs_numbering (GfsBox * box, GfsLinearProblem * lp) {
 
  FttDirection d;
  
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      
      b->type = GFS_BOUNDARY_CENTER_VARIABLE;
      ftt_face_traverse_boundary (b->root, b->d,
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, lp->maxlevel,
				  (FttFaceTraverseFunc) bc_number, lp);   
    }
}

void gfs_get_poisson_problem (GfsDomain * domain, 
			      GfsVariable * rhs, GfsVariable * lhs, GfsVariable * dia,
			      guint dimension,
			      GfsLinearProblem * lp)
{
  GfsVariable * id = gfs_temporary_variable (domain);

  gfs_domain_timer_start (domain, "get_poisson_problem");

  lp->rhs_v = rhs;
  lp->lhs_v = lhs;
  lp->id = id;
  
  /* Cell numbering */
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, lp->maxlevel,
			    (FttCellTraverseFunc) leafs_numbering, lp);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) bc_leafs_numbering, lp);
  /* End - Cell numbering */
 
  /* Creates stencils on the fly */
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS,
			    lp->maxlevel, (FttCellTraverseFunc)  relax_coeff_stencil, lp);

  gfs_domain_homogeneous_bc_stencil (domain, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS,
				     lp->maxlevel, lp->dp, lp->u, lp);
  
  /* ov / v good names ?*/
  
  /*End - Creates stencils on the fly */
  gfs_domain_timer_stop (domain, "get_poisson_problem");
}

typedef struct {
  guint u, rhs, dia, res;
  gint maxlevel;
  gdouble beta, omega;
  FttComponent component;
  guint axi;
} RelaxParams;

static void relax (FttCell * cell, RelaxParams * p)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;

  g.a = GFS_VARIABLE (cell, p->dia);
  g.b = 0.;
  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = neighbor.c[f.d];
    if (f.neighbor) {
      gfs_face_weighted_gradient (&f, &ng, p->u, p->maxlevel);
      g.a += ng.a;
      g.b += ng.b;
    }
  }
  if (g.a > 0.)
    GFS_VARIABLE (cell, p->u) = (g.b - GFS_VARIABLE (cell, p->rhs))/g.a;
  else
    GFS_VARIABLE (cell, p->u) = 0.;
}

static void relax2D (FttCell * cell, RelaxParams * p)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;

  g.a = GFS_VARIABLE (cell, p->dia);
  g.b = 0.;
  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS_2D; f.d++) {
    f.neighbor = neighbor.c[f.d];
    if (f.neighbor) {
      gfs_face_weighted_gradient_2D (&f, &ng, p->u, p->maxlevel);
      g.a += ng.a;
      g.b += ng.b;
    }
  }
  if (g.a > 0.)
    GFS_VARIABLE (cell, p->u) = 
      (1. - p->omega)*GFS_VARIABLE (cell, p->u) 
      + p->omega*(g.b - GFS_VARIABLE (cell, p->rhs))/g.a;
  else
    GFS_VARIABLE (cell, p->u) = 0.;
}

/**
 * gfs_relax:
 * @domain: the domain to relax.
 * @d: number of dimensions (2 or 3).
 * @max_depth: the maximum depth of the domain to relax.
 * @u: the variable to use as left-hand side.
 * @rhs: the variable to use as right-hand side.
 * @dia: the diagonal weight.
 *
 * Apply one pass of a Jacobi relaxation to all the leaf cells of
 * @domain with a level inferior or equal to @max_depth and to all the
 * cells at level @max_depth. The relaxation should converge (if the
 * right-hand-side @rhs verifies the solvability conditions) toward
 * the solution of a Poisson equation for @u at the maximum depth.
 */
void gfs_relax (GfsDomain * domain,
		guint d,
		gint max_depth,
		gdouble omega,
		GfsVariable * u,
		GfsVariable * rhs,
		GfsVariable * dia)
{
  RelaxParams p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (d > 1 && d <= 3);
  g_return_if_fail (u != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (dia != NULL);

  
  p.u = u->i;
  p.rhs = rhs->i;
  p.dia = dia->i;
  p.maxlevel = max_depth;
  p.omega = omega;
  
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, 
			    FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS,
			    max_depth,
			    (FttCellTraverseFunc) (d == 2 ? relax2D : relax), &p);
}

static void residual_set (FttCell * cell, RelaxParams * p)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;

  g.a = GFS_VARIABLE (cell, p->dia);
  g.b = 0.;
  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = neighbor.c[f.d];
    if (f.neighbor) {
      gfs_face_weighted_gradient (&f, &ng, p->u, p->maxlevel);
      g.a += ng.a;
      g.b += ng.b;
    }
  }
  GFS_VARIABLE (cell, p->res) = GFS_VARIABLE (cell, p->rhs) - 
    (g.b - GFS_VARIABLE (cell, p->u)*g.a);
}

static void residual_set2D (FttCell * cell, RelaxParams * p)
{
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;

  g.a = GFS_VARIABLE (cell, p->dia);
  g.b = 0.;
  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS_2D; f.d++) {
    f.neighbor = neighbor.c[f.d];
    if (f.neighbor) {
      gfs_face_weighted_gradient_2D (&f, &ng, p->u, p->maxlevel);
      g.a += ng.a;
      g.b += ng.b;
    }
  }
  GFS_VARIABLE (cell, p->res) = GFS_VARIABLE (cell, p->rhs) - 
    (g.b - GFS_VARIABLE (cell, p->u)*g.a);
}

/**
 * gfs_residual:
 * @domain: a domain.
 * @d: number of dimensions (2 or 3).
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 * @u: the variable to use as left-hand side.
 * @rhs: the variable to use as right-hand side.
 * @dia: the diagonal weight.
 * @res: the variable to use to store the residual.
 *
 * For each cell of @domain, computes the sum of the residual over
 * the volume of the cell for a Poisson equation with @u as
 * left-hand-side and @rhs as right-hand-side. Stores the result in
 * @res.  
 */
void gfs_residual (GfsDomain * domain,
		   guint d,
		   FttTraverseFlags flags,
		   gint max_depth,
		   GfsVariable * u, GfsVariable * rhs, GfsVariable * dia,
		   GfsVariable * res)
{
  RelaxParams p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (d > 1 && d <= 3);
  g_return_if_fail (u != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (dia != NULL);
  g_return_if_fail (res != NULL);

  p.u = u->i;
  p.rhs = rhs->i;
  p.dia = dia->i;
  p.res = res->i;
  p.maxlevel = max_depth;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth,
			    (FttCellTraverseFunc) (d == 2 ? residual_set2D : residual_set), &p);
}

static void reset_coeff (FttCell * cell)
{
  FttDirection d;
  GfsFaceStateVector * f = GFS_STATE (cell)->f;
  
  for (d = 0; d < FTT_NEIGHBORS; d++)
    f[d].v = 0.;
}

typedef struct {
  gdouble lambda2[FTT_DIMENSION];
  GfsFunction * alpha;
  GfsDomain * domain;
} PoissonCoeff;

static void poisson_coeff (FttCellFace * face,
			   PoissonCoeff * p)
{
  gdouble alpha = p->alpha ? gfs_function_face_value (p->alpha, face) : 1.;
  gdouble v = p->lambda2[face->d/2]*alpha*gfs_domain_face_fraction (p->domain, face);

  if (alpha <= 0.) {
    FttVector p;
    ftt_face_pos (face, &p);
    g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	   "alpha is negative (%g) at face (%g,%g,%g).\n"
	   "Please check your definition.",
	   alpha, p.x, p.y, p.z);
  }
  GFS_STATE (face->cell)->f[face->d].v = v;
  
  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v = v;
    break;
  case FTT_FINE_COARSE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v +=
      v/FTT_CELLS_DIRECTION (face->d);
    break;
  default:
    g_assert_not_reached ();
  }
}

static void face_coeff_from_below (FttCell * cell)
{
  FttDirection d;
  GfsFaceStateVector * f = GFS_STATE (cell)->f;
  guint neighbors = 0;

  for (d = 0; d < FTT_NEIGHBORS; d++) {
    FttCellChildren child;
    guint i, n;

    f[d].v = 0.;
    n = ftt_cell_children_direction (cell, d, &child);
    for (i = 0; i < n; i++)
      if (child.c[i])
	f[d].v += GFS_STATE (child.c[i])->f[d].v;
    f[d].v /= n;

    FttCell * neighbor;
    if (f[d].v > 0. && (neighbor = ftt_cell_neighbor (cell, d)) && !GFS_CELL_IS_BOUNDARY (neighbor))
      neighbors++;
  }

  if (neighbors == 1)
    for (d = 0; d < FTT_NEIGHBORS; d++)
      f[d].v = 0.;
}

/**
 * gfs_poisson_coefficients:
 * @domain: a #GfsDomain.
 * @alpha: the inverse of density or %NULL.
 *
 * Initializes the face coefficients for the Poisson equation
 * $\nabla\cdot\alpha\nabla p=\dots$.
 *
 * If @alpha is %NULL, it is taken to be unity.
 */
void gfs_poisson_coefficients (GfsDomain * domain,
			       GfsFunction * alpha)
{
  PoissonCoeff p;
  FttComponent i;

  g_return_if_fail (domain != NULL);

  for (i = 0; i < FTT_DIMENSION; i++) {
    gdouble lambda = (&domain->lambda.x)[i];

    p.lambda2[i] = lambda*lambda;
  }
  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) reset_coeff, NULL);
  p.alpha = alpha;
  p.domain = domain;
  gfs_domain_face_traverse (domain, FTT_XYZ, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) poisson_coeff, &p);
  gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) face_coeff_from_below, NULL);
}

static void tension_coeff (FttCellFace * face, gpointer * data)
{
  gdouble * lambda2 = data[0];
  GfsSourceTensionGeneric * t = data[1];
  GfsVariable * kappa = GFS_SOURCE_TENSION (data[1])->k;
  gdouble alpha = data[2] ? gfs_function_face_value (data[2], face) : 1.;
  gdouble v = lambda2[face->d/2]*alpha*gfs_domain_face_fraction (kappa->domain, face)*
    gfs_function_face_value (t->sigma, face);
  gdouble k1 = GFS_VARIABLE (face->cell, kappa->i);
  gdouble k2 = GFS_VARIABLE (face->neighbor, kappa->i);
#if 0
  gdouble c1 = GFS_VARIABLE (face->cell, t->c->i);
  gdouble c2 = GFS_VARIABLE (face->neighbor, t->c->i);
  gdouble w1 = c1*(1. - c1);
  gdouble w2 = c2*(1. - c2);

  if (w1 + w2 > 0.)
    v *= (w1*k1 + w2*k2)/(w1 + w2);
  else
#endif
  {
    if (k1 < G_MAXDOUBLE) {
      if (k2 < G_MAXDOUBLE)
	v *= (k1 + k2)/2.;
      else
	v *= k1;
    }
    else if (k2 < G_MAXDOUBLE)
      v *= k2;
    else
      v = 1e6;
  }
  g_assert (v <= 1e6);

  if (alpha <= 0.) {
    FttVector p;
    ftt_face_pos (face, &p);
    g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	   "alpha is negative (%g) at face (%g,%g,%g).\n"
	   "Please check your definition.",
	   alpha, p.x, p.y, p.z);
  }
  GFS_STATE (face->cell)->f[face->d].v = v;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v = v;
    break;
  case FTT_FINE_COARSE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v = G_MAXDOUBLE;
    break;
  default:
    g_assert_not_reached ();
  }
}

/**
 * gfs_source_tension_coefficients:
 * @s: a #GfsSourceTension.
 * @domain: a #GfsDomain.
 * @alpha: the inverse of density or %NULL.
 *
 * Initializes the face coefficients with the surface tension term
 * (interface curvature times surface tension coefficient).
 *
 * If @alpha is %NULL, it is taken to be unity.
 */
void gfs_source_tension_coefficients (GfsSourceTension * s,
				      GfsDomain * domain,
				      GfsFunction * alpha)
{
  gdouble lambda2[FTT_DIMENSION];
  gpointer data[3];
  FttComponent i;

  g_return_if_fail (s != NULL);
  g_return_if_fail (domain != NULL);

  for (i = 0; i < FTT_DIMENSION; i++) {
    gdouble lambda = (&domain->lambda.x)[i];

    lambda2[i] = lambda*lambda;
  }
  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) reset_coeff, NULL);
  data[0] = lambda2;
  data[1] = s;
  data[2] = alpha;
  gfs_domain_face_traverse (domain, FTT_XYZ, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) tension_coeff, data);
}

static void correct (FttCell * cell, gpointer * data)
{
  GfsVariable * u = data[0];
  GfsVariable * dp = data[1];
  GFS_VALUE (cell, u) += GFS_VALUE (cell, dp);
}

static void get_from_above (FttCell * parent, GfsVariable * v)
{
  guint level = ftt_cell_level (parent);
  FttCellNeighbors n;
  FttCellChildren child;
  FttComponent c;
  FttVector h;
  guint i;

  ftt_cell_neighbors (parent, &n);
  for (c = 0; c < FTT_DIMENSION; c++) {
    FttCellFace f;
    GfsGradient g;
    gdouble g1, g2;
    
    f.cell = parent;
    f.d = 2*c;
    f.neighbor = n.c[f.d];
    gfs_face_gradient (&f, &g, v->i, level);
    g1 = g.b - g.a*GFS_VARIABLE (parent, v->i);
    f.d = 2*c + 1;
    f.neighbor = n.c[f.d];
    gfs_face_gradient (&f, &g, v->i, level);
    g2 = g.b - g.a*GFS_VARIABLE (parent, v->i);
    (&h.x)[c] = (g1 - g2)/2.;
  }

  ftt_cell_children (parent, &child);
  for (i = 0; i < FTT_CELLS; i++) 
    if (child.c[i]) {
      FttVector p;
      
      GFS_VARIABLE (child.c[i], v->i) = GFS_VARIABLE (parent, v->i);
      ftt_cell_relative_pos (child.c[i], &p);
      for (c = 0; c < FTT_DIMENSION; c++)
	GFS_VARIABLE (child.c[i], v->i) += (&p.x)[c]*(&h.x)[c];
    }
}

static void get_from_below_3D (FttCell * cell, const GfsVariable * v)
{
  gdouble val = 0.;
  guint i;
  FttCellChildren child;

  ftt_cell_children (cell, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i])
      val += GFS_VARIABLE (child.c[i], v->i);
  GFS_VARIABLE (cell, v->i) = val/2.;
}

static void get_from_below_2D (FttCell * cell, const GfsVariable * v)
{
  gdouble val = 0.;
  guint i;
  FttCellChildren child;

  ftt_cell_children (cell, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i])
      val += GFS_VARIABLE (child.c[i], v->i);
  GFS_VARIABLE (cell, v->i) = val;
}

typedef struct {
  GfsVariable * s, * r, * u, * v;
  gdouble srs, rs2, beta;
} MRSData;

static void compute_beta (FttCell * cell, MRSData * data)
{
  gdouble rs = GFS_VALUE (cell, data->r) - GFS_VALUE (cell, data->s);
  data->rs2 += rs*rs;
  data->srs -= GFS_VALUE (cell, data->s)*rs;
}

static void update_sv (FttCell * cell, MRSData * data)
{
  GFS_VALUE (cell, data->s) += data->beta*(GFS_VALUE (cell, data->r) - GFS_VALUE (cell, data->s));
  GFS_VALUE (cell, data->v) += data->beta*(GFS_VALUE (cell, data->u) - GFS_VALUE (cell, data->v));
  GFS_VALUE (cell, data->r) = GFS_VALUE (cell, data->s);
  GFS_VALUE (cell, data->u) = GFS_VALUE (cell, data->v);
}

static void relax_loop (GfsDomain * domain, 
			GfsVariable * dp, GfsVariable * u, 
			RelaxParams * q, guint nrelax,
			guint dimension)
{
  guint n;

  gfs_domain_homogeneous_bc (domain,
			     FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, q->maxlevel, 
			     dp, u);
  for (n = 0; n < nrelax - 1; n++)
    gfs_traverse_and_homogeneous_bc (domain, FTT_PRE_ORDER, 
				     FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, q->maxlevel,
				     (FttCellTraverseFunc) (dimension == 2 ? relax2D : relax), q,
				     dp, u);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, 
			    FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, q->maxlevel,
			    (FttCellTraverseFunc) (dimension == 2 ? relax2D : relax), q);
}

/**
 * gfs_poisson_cycle:
 * @domain: the domain on which to solve the Poisson equation.
 * @p: the #GfsMultilevelParams.
 * @u: the variable to use as left-hand side.
 * @rhs: the variable to use as right-hand side.
 * @dia: the diagonal weight.
 * @res: the residual.
 *
 * Apply one multigrid iteration to the Poisson equation defined by @u
 * and @rhs.
 *
 * The initial value of @res on the leaves of @root must be set to
 * the residual of the Poisson equation (using gfs_residual()).
 *
 * The face coefficients must be set using gfs_poisson_coefficients().
 *
 * The values of @u on the leaf cells are updated as well as the values
 * of @res (i.e. the cell tree is ready for another iteration).
 *
 * Returns a gboolean : TRUE if the function uses its own convergence
 * criterion and satisfied it. FALSE otherwise.
 */
void gfs_poisson_cycle (GfsDomain * domain,
			   GfsMultilevelParams * p,
			   GfsVariable * u,
			   GfsVariable * rhs,
			   GfsVariable * dia,
			   GfsVariable * res)
{
  guint l, nrelax, minlevel;
  GfsVariable * dp;
  gpointer data[2];
  
  g_return_if_fail (domain != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (p->dimension > 1 && p->dimension <= 3);
  g_return_if_fail (u != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (dia != NULL);
  g_return_if_fail (res != NULL);

  dp = gfs_temporary_variable (domain);
  minlevel = MAX (domain->rootlevel, p->minlevel);

  /* compute residual on non-leafs cells */
  gfs_domain_cell_traverse (domain, 
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    p->dimension == 2 ? (FttCellTraverseFunc) get_from_below_2D : 
			    (FttCellTraverseFunc) get_from_below_3D,
			    res);

  /* relax top level */
  nrelax = p->nrelax;
  for (l = minlevel; l < p->depth; l++)
    nrelax *= p->erelax;

  RelaxParams q;
  q.u = dp->i;
  q.rhs = res->i;
  q.dia = dia->i;
  q.maxlevel = minlevel;
  q.omega = p->omega;
  
  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, q.maxlevel,
			    (FttCellTraverseFunc) gfs_cell_reset, dp);
  relax_loop (domain, dp, u, &q, nrelax, p->dimension);
  nrelax /= p->erelax;

  /* relax from top to bottom */
  for (q.maxlevel = minlevel + 1; q.maxlevel <= p->depth; q.maxlevel++, nrelax /= p->erelax) {
    /* get initial guess from coarser grid */ 
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_NON_LEAFS, 
			      q.maxlevel - 1,
			      (FttCellTraverseFunc) get_from_above, dp);
    relax_loop (domain, dp, u, &q, nrelax, p->dimension);
  }
  /* correct on leaf cells */
  data[0] = u;
  data[1] = dp;
  gfs_traverse_and_bc (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) correct, data,
		       u, u);
  /* compute new residual on leaf cells */
  gfs_residual (domain, p->dimension, FTT_TRAVERSE_LEAFS, -1, u, rhs, dia, res);

  gts_object_destroy (GTS_OBJECT (dp));
}

/**
 * gfs_poisson_solver:
 * @domain: the domain over which the poisson problem is solved.
 * @par: the parameters of the poisson problem.
 * @lhs: the variable to use as left-hand side.
 * @rhs: the variable to use as right-hand side.
 * @res: the varaible to store the residual
 * @dia: the diagonal weight.
 * @dt: the length of the time-step.
 *
 * Solve the poisson problem over domain using the Gerris' native
 * multigrid poisson solver.
 */
void gfs_poisson_solve (GfsDomain * domain, GfsMultilevelParams * par,
			GfsVariable * lhs, GfsVariable * rhs, GfsVariable * res,
			GfsVariable * dia, gdouble dt)
{
  gfs_domain_timer_start (domain, "poisson_solve");

  /* calculates the intial residual and its norm */
  gfs_residual (domain, par->dimension, FTT_TRAVERSE_LEAFS, -1, lhs, rhs, dia, res);
  par->residual_before = par->residual = 
    gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, dt, res);

  gdouble res_max_before = par->residual.infty;

  while (par->niter < par->nitermin ||
	 (par->residual.infty > par->tolerance && par->niter < par->nitermax)) {

    /* Does one iteration */
    gfs_poisson_cycle (domain, par, lhs, rhs, dia, res);
    
    par->residual = gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, dt, res);

    if (par->residual.infty == res_max_before) /* convergence has stopped!! */
      break;
    if (par->residual.infty > res_max_before/1.1 && par->minlevel < par->depth)
      par->minlevel++;
    res_max_before = par->residual.infty;
    par->niter++;
  }
  gfs_domain_timer_stop (domain, "poisson_solve");
}

typedef struct {
  GfsSourceDiffusion * d;
  gdouble lambda2[FTT_DIMENSION];
  gdouble dt;
  GfsVariable * rhoc, * axi;
  GfsFunction * alpha;
  GfsDomain * domain;
} DiffusionCoeff;

static void diffusion_coef (FttCellFace * face, DiffusionCoeff * c)
{
  gdouble v = 
    c->lambda2[face->d/2]*c->dt*
    gfs_source_diffusion_face (c->d, face)*
    gfs_domain_face_fraction (c->domain, face);

  GFS_STATE (face->cell)->f[face->d].v = v;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v = v;
    break;
  case FTT_FINE_COARSE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v +=
      v/FTT_CELLS_DIRECTION (face->d);
    break;
  default:
    g_assert_not_reached ();
  }
}

static void diffusion_mixed_coef (FttCell * cell, DiffusionCoeff * c)
{
  reset_coeff (cell);
  if (GFS_IS_MIXED (cell))
    GFS_STATE (cell)->solid->v = 
      c->dt*gfs_domain_solid_metric (c->domain, cell)*gfs_source_diffusion_cell (c->d, cell);
  if (c->rhoc) {
    gdouble rho = c->alpha ? 1./gfs_function_value (c->alpha, cell) : 1.;
    if (rho <= 0.) {
      FttVector p;
      ftt_cell_pos (cell, &p);
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	     "density is negative (%g) at cell (%g,%g,%g).\n"
	     "Please check your definition of alpha.",
	     rho, p.x, p.y, p.z);
    }
    gdouble f = gfs_domain_cell_fraction (c->domain, cell);
    GFS_VALUE (cell, c->rhoc) = rho*f;

    if (c->axi) {
      FttVector p;
      gfs_cell_cm (cell, &p);
      GFS_VALUE (cell, c->axi) = 2.*c->dt*gfs_source_diffusion_cell (c->d, cell)/(rho*p.y*p.y);
    }
  }
}

/**
 * gfs_diffusion_coefficients:
 * @domain: a #GfsDomain.
 * @d: a #GfsSourceDiffusion.
 * @dt: the time-step.
 * @rhoc: where to store the mass.
 * @axi: where to store the axisymmetric term (or %NULL).
 * @alpha: the inverse of density or %NULL.
 * @beta: the implicitness parameter (0.5 Crank-Nicholson, 1. backward Euler).
 *
 * Initializes the face coefficients for the diffusion equation.
 */
void gfs_diffusion_coefficients (GfsDomain * domain,
				 GfsSourceDiffusion * d,
				 gdouble dt,
				 GfsVariable * rhoc,
				 GfsVariable * axi,
				 GfsFunction * alpha,
				 gdouble beta)
{
  DiffusionCoeff coef;
  FttComponent i;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (d != NULL);
  g_return_if_fail (beta >= 0.5 && beta <= 1.);

  for (i = 0; i < FTT_DIMENSION; i++) {
    gdouble lambda = (&domain->lambda.x)[i];

    coef.lambda2[i] = lambda*lambda;
  }
  coef.d = d;
  coef.dt = beta*dt;
  coef.rhoc = rhoc;
  coef.alpha = alpha;
  coef.domain = domain;
  coef.axi = axi;
  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) diffusion_mixed_coef, &coef);
  gfs_domain_face_traverse (domain, FTT_XYZ, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) diffusion_coef, &coef);
  gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) face_coeff_from_below, 
			    NULL);
}

static void diffusion_rhs (FttCell * cell, RelaxParams * p)
{
  gdouble f, h, val;
  FttCellNeighbors neighbor;
  FttCellFace face;
  
  if (GFS_IS_MIXED (cell)) {
    if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
      f = gfs_cell_dirichlet_gradient_flux (cell, p->u, -1, GFS_STATE (cell)->solid->fv);
    else
      f = GFS_STATE (cell)->solid->fv;
  }
  else
    f = 0.; /* Neumann condition by default */
  h = ftt_cell_size (cell);
  val = GFS_VARIABLE (cell, p->u);
  face.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++) {
    GfsGradient g;

    face.neighbor = neighbor.c[face.d];
    gfs_face_gradient_flux (&face, &g, p->u, -1);
    if (face.d/2 == p->component) {
      g.a *= 2.;
      g.b *= 2.;
    }
    f += g.b - g.a*val;
  }
  GFS_VARIABLE (cell, p->rhs) += p->beta*f/(h*h*GFS_VARIABLE (cell, p->dia));
  if (p->axi)
    GFS_VARIABLE (cell, p->rhs) -= val*p->beta*GFS_VARIABLE (cell, p->axi);
}

/**
 * gfs_diffusion_rhs:
 * @domain: a #GfsDomain.
 * @v: a #GfsVariable.
 * @rhs: a #GfsVariable.
 * @rhoc: the mass.
 * @axi: the axisymmetric term.
 * @beta: the implicitness parameter (0.5 Crank-Nicholson, 1. backward Euler).
 *
 * Adds to the @rhs variable of @cell the right-hand side of the
 * diffusion equation for variable @v.
 *
 * The diffusion coefficients must have been already set using
 * gfs_diffusion_coefficients().
 */
void gfs_diffusion_rhs (GfsDomain * domain, 
			GfsVariable * v, GfsVariable * rhs, 
			GfsVariable * rhoc, GfsVariable * axi,
			gdouble beta)
{
  RelaxParams p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (rhoc != NULL);
  g_return_if_fail (beta >= 0.5 && beta <= 1.);

  p.u = v->i;
  p.rhs = rhs->i;
  p.dia = rhoc->i;
  p.beta = (1. - beta)/beta;
  p.component = GFS_IS_AXI (domain) ? v->component : FTT_DIMENSION;
  p.axi = axi ? axi->i : FALSE;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) diffusion_rhs, &p);
}

static void diffusion_relax (FttCell * cell, RelaxParams * p)
{
  gdouble a;
  GfsGradient g = { 0., 0. };
  gdouble h = ftt_cell_size (cell);
  FttCellNeighbors neighbor;
  FttCellFace face;

  a = GFS_VARIABLE (cell, p->dia);
  if (GFS_IS_MIXED (cell) && ((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
    g.b = gfs_cell_dirichlet_gradient_flux (cell, p->u, p->maxlevel, 0.);

  face.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++) {
    GfsGradient ng;

    face.neighbor = neighbor.c[face.d];
    gfs_face_gradient_flux (&face, &ng, p->u, p->maxlevel);
    if (face.d/2 == p->component) {
      ng.a *= 2.;
      ng.b *= 2.;
    }
    g.a += ng.a;
    g.b += ng.b;
  }
  a *= h*h;
  g_assert (a > 0.);
  g.a = 1. + g.a/a;
  if (p->axi)
    g.a += GFS_VARIABLE (cell, p->axi);
  g.b = GFS_VARIABLE (cell, p->res) + g.b/a;
  g_assert (g.a > 0.);
  GFS_VARIABLE (cell, p->u) = g.b/g.a;
}

static void diffusion_residual (FttCell * cell, RelaxParams * p)
{
  gdouble a;
  GfsGradient g = { 0., 0. };
  gdouble h;
  FttCellNeighbors neighbor;
  FttCellFace face;

  h = ftt_cell_size (cell);
  a = GFS_VARIABLE (cell, p->dia);
  if (GFS_IS_MIXED (cell)) {
    if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
      g.b = gfs_cell_dirichlet_gradient_flux (cell, p->u, -1, GFS_STATE (cell)->solid->fv);
    else
      g.b = GFS_STATE (cell)->solid->fv;
  }

  face.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++) {
    GfsGradient ng;

    face.neighbor = neighbor.c[face.d];
    gfs_face_gradient_flux (&face, &ng, p->u, -1);
    if (face.d/2 == p->component) {
      ng.a *= 2.;
      ng.b *= 2.;
    }
    g.a += ng.a;
    g.b += ng.b;
  }
  a *= h*h;
  g_assert (a > 0.);
  g.a = 1. + g.a/a;
  if (p->axi)
    g.a += GFS_VARIABLE (cell, p->axi);
  g.b = GFS_VARIABLE (cell, p->rhs) + g.b/a;
  GFS_VARIABLE (cell, p->res) = g.b - g.a*GFS_VARIABLE (cell, p->u);
}

/**
 * gfs_diffusion_residual:
 * @domain: a #GfsDomain.
 * @u: the variable to use as left-hand side.
 * @rhs: the right-hand side.
 * @rhoc: the mass.
 * @axi: the axisymmetric term.
 * @res: the residual.
 *
 * Sets the @res variable of each leaf cell of @domain to the residual
 * of the diffusion equation for @v.
 *
 * The diffusion coefficients must have been set using
 * gfs_diffusion_coefficients() and the right-hand side using
 * gfs_diffusion_rhs().
 */
void gfs_diffusion_residual (GfsDomain * domain,
			     GfsVariable * u,
			     GfsVariable * rhs,
			     GfsVariable * rhoc,
			     GfsVariable * axi,
			     GfsVariable * res)
{
  RelaxParams p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (u != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (rhoc != NULL);
  g_return_if_fail (res != NULL);

  p.u = u->i;
  p.rhs = rhs->i;
  p.dia = rhoc->i;
  p.res = res->i;
  p.component = GFS_IS_AXI (domain) ? u->component : FTT_DIMENSION;
  p.axi = axi ? axi->i : FALSE;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) diffusion_residual, &p);
}

static void diffusion_relax_loop (GfsDomain * domain, 
				  GfsVariable * dp, GfsVariable * u,
				  RelaxParams * p, guint nrelax)
{
  guint n;
  gfs_domain_homogeneous_bc (domain, 
			     FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, p->maxlevel,
			     dp, u);
  for (n = 0; n < nrelax - 1; n++)
    gfs_traverse_and_homogeneous_bc (domain, FTT_PRE_ORDER, 
				     FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, p->maxlevel,
				     (FttCellTraverseFunc) diffusion_relax, p,
				     dp, u);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, 
			    FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, p->maxlevel,
			    (FttCellTraverseFunc) diffusion_relax, p);
}

/**
 * gfs_diffusion_cycle:
 * @domain: the domain on which to solve the diffusion equation.
 * @levelmin: the top level of the multigrid hierarchy.
 * @depth: the total depth of the domain.
 * @nrelax: the number of relaxations to apply at each level.
 * @u: the variable to use as left-hand side.
 * @rhs: the right-hand side.
 * @rhoc: the mass.
 * @axi: the axisymmetric term.
 * @res: the residual.
 *
 * Apply one multigrid iteration to the diffusion equation for @u.
 *
 * The initial value of @res on the leaves of @root must be set to
 * the residual of the diffusion equation using gfs_diffusion_residual().
 *
 * The diffusion coefficients must be set using gfs_diffusion_coefficients().
 *
 * The values of @u on the leaf cells are updated as well as the values
 * of @res (i.e. the cell tree is ready for another iteration).
 */
void gfs_diffusion_cycle (GfsDomain * domain,
			  guint levelmin,
			  guint depth,
			  guint nrelax,
			  GfsVariable * u,
			  GfsVariable * rhs,
			  GfsVariable * rhoc,
			  GfsVariable * axi,
			  GfsVariable * res)
{
  GfsVariable * dp;
  RelaxParams p;
  gpointer data[2];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (u != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (rhoc != NULL);
  g_return_if_fail (res != NULL);

  dp = gfs_temporary_variable (domain);
  dp->component = u->component;

  /* compute residual on non-leafs cells */
  gfs_domain_cell_traverse (domain, 
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_get_from_below_intensive, res);

  /* relax top level */
  p.maxlevel = levelmin;
  p.u = dp->i;
  p.res = res->i;
  p.dia = rhoc->i;
  p.component = GFS_IS_AXI (domain) ? u->component : FTT_DIMENSION;
  p.axi = axi ? axi->i : FALSE;

  gfs_domain_cell_traverse (domain, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, levelmin,
			    (FttCellTraverseFunc) gfs_cell_reset, dp);
  diffusion_relax_loop (domain, dp, u, &p, 10*nrelax);
  /* relax from top to bottom */
  for (p.maxlevel = levelmin + 1; p.maxlevel <= depth; p.maxlevel++) {
    /* get initial guess from coarser grid */ 
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_NON_LEAFS,
			      p.maxlevel - 1,
			      (FttCellTraverseFunc) get_from_above, dp);
    diffusion_relax_loop (domain, dp, u, &p, nrelax);
  }
  /* correct on leaf cells */
  data[0] = u;
  data[1] = dp;
  gfs_traverse_and_bc (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) correct, data,
		       u, u);
  /* compute new residual on leaf cells */
  gfs_diffusion_residual (domain, u, rhs, rhoc, axi, res);

  gts_object_destroy (GTS_OBJECT (dp));
}

