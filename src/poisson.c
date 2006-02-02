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
	   "  beta      = %g\n"
	   "}",
	   par->tolerance,
	   par->nrelax,
	   par->erelax,
	   par->minlevel,
	   par->nitermax,
	   par->nitermin,
	   par->weighted,
	   par->beta);
}

void gfs_multilevel_params_init (GfsMultilevelParams * par)
{
  g_return_if_fail (par != NULL);

  par->tolerance = 1e-3;
  par->nrelax    = 4;
  par->erelax    = 1;
  par->minlevel  = 0;
  par->nitermax  = 100;
  par->nitermin  = 0;

  par->dimension = FTT_DIMENSION;
  par->weighted = FALSE;
  par->beta = 0.5;
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

#define FACE_GRADIENT(x, y, z, u) (gfs_face_weighted_gradient (x, y, z, u))
//#define FACE_GRADIENT(x, y, z, u) (gfs_face_gradient_flux_centered (x, y, z, u))

void gfs_face_gradient_flux_centered (const FttCellFace * face,
				      GfsGradient * g,
				      guint v,
				      gint max_level)
{
  guint level;

  g_return_if_fail (face != NULL);

  g->a = g->b = 0.;
  if (face->neighbor == NULL)
    return;

  level = ftt_cell_level (face->cell);
  if (ftt_cell_level (face->neighbor) < level) {
    g_assert_not_implemented ();
#if 0
    /* neighbor is at a shallower level */
    Gradient gcf;
    gdouble w = GFS_STATE (face->cell)->f[face->d].v;

    gcf = gradient_flux_fine_coarse (face, v, max_level);
    g->a = w*gcf.a;
    g->b = w*(gcf.b*GFS_VARIABLE (face->neighbor, v) + gcf.c);
#endif
  }
  else {
    if (level == max_level || FTT_CELL_IS_LEAF (face->neighbor)) {
      /* neighbor is at the same level */
      gdouble w = GFS_STATE (face->cell)->f[face->d].v;

      if (!GFS_IS_MIXED (face->cell) || !GFS_IS_MIXED (face->neighbor)) {
	g->a = w;
	g->b = w*GFS_VARIABLE (face->neighbor, v);
      }
      else {
	FttComponent c = face->d/2;
	FttComponent cp = FTT_ORTHOGONAL_COMPONENT (c);
	FttCell * n[4];
	gdouble ye = (1. - GFS_FACE_FRACTION (face))/2.;
	GfsSolidVector * s = GFS_STATE (face->cell)->solid;

	n[0] = face->cell; n[1] = face->neighbor;
	if ((s->s[2*cp] == 1. && s->s[2*cp + 1] < 1.) ||
	    (s->s[2*cp + 1] == 0. && 
	     s->s[2*cp] < 1. && s->s[2*cp] > 0.)) {
	  n[2] = ftt_cell_neighbor (n[0], 2*cp);
	  n[3] = ftt_cell_neighbor (n[1], 2*cp);
	}
	else if ((s->s[2*cp + 1] == 1. && s->s[2*cp] < 1.) ||
		 (s->s[2*cp] == 0. && 
		  s->s[2*cp + 1] < 1. && s->s[2*cp + 1] > 0.)) {
	  n[2] = ftt_cell_neighbor (n[0], 2*cp + 1);
	  n[3] = ftt_cell_neighbor (n[1], 2*cp + 1);
	}
	else {
	  g->a = w;
	  g->b = w*GFS_VARIABLE (n[1], v);
	  return;
	}

	//	g_assert (n[2] && n[3]);
	if (n[2] && n[3]) {
	  g->a = w*(1. - ye);
	  g->b = w*(1. - ye)*GFS_VARIABLE (n[1], v) + 
	    w*ye*(GFS_VARIABLE (n[3], v) - GFS_VARIABLE (n[2], v));
	}
	else {
	  g->a = w;
	  g->b = w*GFS_VARIABLE (n[1], v);
	}
      }
    }
    else {
    g_assert_not_implemented ();
#if 0
      /* neighbor is at a deeper level */
      FttCellChildren children;
      FttCellFace f;
      guint i, n;
      
      f.d = FTT_OPPOSITE_DIRECTION (face->d);
      n = ftt_cell_children_direction (face->neighbor, f.d, &children);
      f.neighbor = face->cell;
      for (i = 0; i < n; i++) {
	Gradient gcf;
	gdouble w;

	f.cell = children.c[i];
	w = GFS_STATE (f.cell)->f[f.d].v;
	
	/* check for mixed cell refinement violation (topology.fig) */
	g_assert (f.cell);
	
	gcf = gradient_fine_coarse (&f, v, max_level);
	g->a += w*gcf.b;
	g->b += w*(gcf.a*GFS_VARIABLE (f.cell, v) - gcf.c);
      }
#endif
#if (!FTT_2D && !FTT_2D3)
      g->a /= 2.;
      g->b /= 2.;
#endif /* not 2D and not 2D3 */
    }
  }
}

typedef struct {
  guint u, rhs, dia, res;
  gint maxlevel;
  gdouble beta;
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
      FACE_GRADIENT (&f, &ng, p->u, p->maxlevel);
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
      FACE_GRADIENT (&f, &ng, p->u, p->maxlevel);
      g.a += ng.a;
      g.b += ng.b;
    }
  }
  if (g.a > 0.)
    GFS_VARIABLE (cell, p->u) = (g.b - GFS_VARIABLE (cell, p->rhs))/g.a;
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
      FACE_GRADIENT (&f, &ng, p->u, p->maxlevel);
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
      FACE_GRADIENT (&f, &ng, p->u, p->maxlevel);
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

static void poisson_coeff (FttCellFace * face, gdouble * lambda2)
{
  GfsStateVector * s = GFS_STATE (face->cell);
  gdouble v = lambda2[face->d/2];

  if (GFS_IS_MIXED (face->cell))
    v *= s->solid->s[face->d];
  s->f[face->d].v = v;

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

static void reset_alpha_coeff (FttCell * cell, gpointer * data)
{
  FttDirection d;
  GfsFaceStateVector * f = GFS_STATE (cell)->f;
  GfsFunction * alpha = data[0];
  GfsVariable * a = data[1];
  
  for (d = 0; d < FTT_NEIGHBORS; d++)
    f[d].v = 0.;
  GFS_VARIABLE (cell, a->i) = gfs_function_value (alpha, cell);
}

static void poisson_alpha_coeff (FttCellFace * face,
				 gpointer * data)
{
  gdouble * lambda2 = data[0];
  GfsVariable * alpha = data[1];
  gdouble v = lambda2[face->d/2];
  GfsStateVector * s = GFS_STATE (face->cell);

  if (GFS_IS_MIXED (face->cell))
    v *= s->solid->s[face->d];
  v *= gfs_face_interpolated_value (face, alpha->i);
  s->f[face->d].v = v;

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

  for (d = 0; d < FTT_NEIGHBORS; d++) {
    FttCellChildren child;
    guint i, n;

    f[d].v = 0.;
    n = ftt_cell_children_direction (cell, d, &child);
    for (i = 0; i < n; i++)
      if (child.c[i])
	f[d].v += GFS_STATE (child.c[i])->f[d].v;
    f[d].v /= n;
  }
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
  gdouble lambda2[FTT_DIMENSION];
  FttComponent i;

  g_return_if_fail (domain != NULL);

  for (i = 0; i < FTT_DIMENSION; i++) {
    gdouble lambda = (&domain->lambda.x)[i];

    lambda2[i] = lambda*lambda;
  }
  if (alpha == NULL) {
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) reset_coeff, NULL);
    gfs_domain_face_traverse (domain, FTT_XYZ, 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) poisson_coeff, lambda2);
  }
  else {
    gpointer data[2];

    data[0] = alpha;
    data[1] = gfs_temporary_variable (domain);
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) reset_alpha_coeff, data);
    data[0] = lambda2;
    gfs_domain_face_traverse (domain, FTT_XYZ, 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) poisson_alpha_coeff, 
			      data);
    gts_object_destroy (data[1]);
  }
  gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) face_coeff_from_below, NULL);
}

static void correct (FttCell * cell, gpointer * data)
{
  GfsVariable * u = data[0];
  GfsVariable * dp = data[1];
  GFS_VARIABLE (cell, u->i) += GFS_VARIABLE (cell, dp->i);
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
 */
void gfs_poisson_cycle (GfsDomain * domain,
			GfsMultilevelParams * p,
			GfsVariable * u,
			GfsVariable * rhs,
			GfsVariable * dia,
			GfsVariable * res)
{
  guint n, l, nrelax, minlevel;
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
			    (FttCellTraverseFunc) gfs_get_from_below_extensive, res);

  /* relax top level */
  gfs_domain_cell_traverse (domain, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, minlevel,
			    (FttCellTraverseFunc) gfs_cell_reset, dp);
  nrelax = p->nrelax;
  for (l = minlevel; l < p->depth; l++)
    nrelax *= p->erelax;
  for (n = 0; n < nrelax; n++) {
    gfs_domain_homogeneous_bc (domain,
			       FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS,
			       minlevel, dp, u);
    gfs_relax (domain, p->dimension, minlevel, dp, res, dia);
  }
  nrelax /= p->erelax;

  /* relax from top to bottom */
  for (l = minlevel + 1; l <= p->depth; l++, nrelax /= p->erelax) {
    /* get initial guess from coarser grid */ 
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_NON_LEAFS, l - 1,
			      (FttCellTraverseFunc) get_from_above, dp);
    for (n = 0; n < nrelax; n++) {
      gfs_domain_homogeneous_bc (domain, 
				 FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS,
				 l, dp, u);
      gfs_relax (domain, p->dimension, l, dp, res, dia);
    }
  }
  /* correct on leaf cells */
  data[0] = u;
  data[1] = dp;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) correct, data);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, u);
  /* compute new residual on leaf cells */
  gfs_residual (domain, p->dimension, FTT_TRAVERSE_LEAFS, -1, u, rhs, dia, res);

  gts_object_destroy (GTS_OBJECT (dp));
}

typedef struct {
  GfsSourceDiffusion * d;
  gdouble lambda2[FTT_DIMENSION];
  gdouble dt;
  GfsVariable * dia;
  GfsFunction * alpha;
} DiffusionCoeff;

static void diffusion_coef (FttCellFace * face, DiffusionCoeff * c)
{
  GfsStateVector * s = GFS_STATE (face->cell);
  gdouble v = c->lambda2[face->d/2]*c->dt*gfs_source_diffusion_face (c->d, face);

  if (GFS_IS_MIXED (face->cell))
    v *= s->solid->s[face->d];
  s->f[face->d].v = v;

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
    GFS_STATE (cell)->solid->v = c->dt*gfs_source_diffusion_cell (c->d, cell);
  GFS_VARIABLE (cell, c->dia->i) = c->alpha ? 1./gfs_function_value (c->alpha, cell) : 1.;
}

/**
 * gfs_diffusion_coefficients:
 * @domain: a #GfsDomain.
 * @d: a #GfsSourceDiffusion.
 * @dt: the time-step.
 * @dia: where to store the diagonal weight.
 * @alpha: the inverse of density or %NULL.
 * @beta: the implicitness parameter (0.5 Crank-Nicholson, 1. backward Euler).
 *
 * Initializes the face coefficients for the diffusion equation.
 */
void gfs_diffusion_coefficients (GfsDomain * domain,
				 GfsSourceDiffusion * d,
				 gdouble dt,
				 GfsVariable * dia,
				 GfsFunction * alpha,
				 gdouble beta)
{
  DiffusionCoeff coef;
  FttComponent i;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (d != NULL);
  g_return_if_fail (dia != NULL);
  g_return_if_fail (beta >= 0.5 && beta <= 1.);

  for (i = 0; i < FTT_DIMENSION; i++) {
    gdouble lambda = (&domain->lambda.x)[i];

    coef.lambda2[i] = lambda*lambda;
  }
  coef.d = d;
  coef.dt = beta*dt;
  coef.dia = dia;
  coef.alpha = alpha;
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
  gdouble a, f, h, val;
  FttCellNeighbors neighbor;
  FttCellFace face;
  
  if (GFS_IS_MIXED (cell)) {
    a = GFS_STATE (cell)->solid->a*GFS_VARIABLE (cell, p->dia);
    if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
      f = gfs_cell_dirichlet_gradient_flux (cell, p->u, -1, GFS_STATE (cell)->solid->fv);
    else
      f = GFS_STATE (cell)->solid->fv;
  }
  else {
    a = GFS_VARIABLE (cell, p->dia);
    f = 0.; /* Neumann condition by default */
  }
  h = ftt_cell_size (cell);
  val = GFS_VARIABLE (cell, p->u);
  face.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++) {
    GfsGradient g;

    face.neighbor = neighbor.c[face.d];
    gfs_face_gradient_flux (&face, &g, p->u, -1);
    f += g.b - g.a*val;
  }
  GFS_VARIABLE (cell, p->rhs) += val + p->beta*f/(h*h*a);
}

/**
 * gfs_diffusion_rhs:
 * @domain: a #GfsDomain.
 * @v: a #GfsVariable.
 * @rhs: a #GfsVariable.
 * @dia: the diagonal weight.
 * @beta: the implicitness parameter (0.5 Crank-Nicholson, 1. backward Euler).
 *
 * Adds to the @rhs variable of @cell the right-hand side of the
 * diffusion equation for variable @v.
 *
 * The diffusion coefficients must have been already set using
 * gfs_diffusion_coefficients().
 */
void gfs_diffusion_rhs (GfsDomain * domain, 
			GfsVariable * v, GfsVariable * rhs, GfsVariable * dia,
			gdouble beta)
{
  RelaxParams p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (dia != NULL);
  g_return_if_fail (beta >= 0.5 && beta <= 1.);

  p.u = v->i;
  p.rhs = rhs->i;
  p.dia = dia->i;
  p.beta = (1. - beta)/beta;
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

  if (GFS_IS_MIXED (cell)) {
    a = GFS_STATE (cell)->solid->a*GFS_VARIABLE (cell, p->dia);
    if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
      g.b = gfs_cell_dirichlet_gradient_flux (cell, p->u, p->maxlevel, 0.);
  }
  else
    a = GFS_VARIABLE (cell, p->dia);

  face.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++) {
    GfsGradient ng;

    face.neighbor = neighbor.c[face.d];
    gfs_face_gradient_flux (&face, &ng, p->u, p->maxlevel);
    g.a += ng.a;
    g.b += ng.b;
  }
  a *= h*h;
  g_assert (a > 0.);
  g.a = 1. + g.a/a;
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
  if (GFS_IS_MIXED (cell)) {
    a = GFS_STATE (cell)->solid->a*GFS_VARIABLE (cell, p->dia);
    if (((cell)->flags & GFS_FLAG_DIRICHLET) != 0)
      g.b = gfs_cell_dirichlet_gradient_flux (cell, p->u, -1, GFS_STATE (cell)->solid->fv);
    else
      g.b = GFS_STATE (cell)->solid->fv;
  }
  else
    a = GFS_VARIABLE (cell, p->dia);

  face.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (face.d = 0; face.d < FTT_NEIGHBORS; face.d++) {
    GfsGradient ng;

    face.neighbor = neighbor.c[face.d];
    gfs_face_gradient_flux (&face, &ng, p->u, -1);
    g.a += ng.a;
    g.b += ng.b;
  }
  a *= h*h;
  g_assert (a > 0.);
  g.a = 1. + g.a/a;
  g.b = GFS_VARIABLE (cell, p->rhs) + g.b/a;
  GFS_VARIABLE (cell, p->res) = g.b - g.a*GFS_VARIABLE (cell, p->u);
}

/**
 * gfs_diffusion_residual:
 * @domain: a #GfsDomain.
 * @u: the variable to use as left-hand side.
 * @rhs: the right-hand side.
 * @dia: the diagonal weight.
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
			     GfsVariable * dia,
			     GfsVariable * res)
{
  RelaxParams p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (u != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (dia != NULL);
  g_return_if_fail (res != NULL);

  p.u = u->i;
  p.rhs = rhs->i;
  p.dia = dia->i;
  p.res = res->i;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) diffusion_residual, &p);
}

/**
 * gfs_diffusion_cycle:
 * @domain: the domain on which to solve the diffusion equation.
 * @levelmin: the top level of the multigrid hierarchy.
 * @depth: the total depth of the domain.
 * @nrelax: the number of relaxations to apply at each level.
 * @u: the variable to use as left-hand side.
 * @rhs: the right-hand side.
 * @dia: the diagonal weight.
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
			  GfsVariable * dia,
			  GfsVariable * res)
{
  guint n;
  GfsVariable * dp;
  RelaxParams p;
  gpointer data[2];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (u != NULL);
  g_return_if_fail (rhs != NULL);
  g_return_if_fail (dia != NULL);
  g_return_if_fail (res != NULL);

  dp = gfs_temporary_variable (domain);

  /* compute residual on non-leafs cells */
  gfs_domain_cell_traverse (domain, 
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_get_from_below_intensive, res);

  /* relax top level */
  gfs_domain_cell_traverse (domain, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, levelmin,
			    (FttCellTraverseFunc) gfs_cell_reset, dp);
  p.maxlevel = levelmin;
  p.u = dp->i;
  p.res = res->i;
  p.dia = dia->i;
  for (n = 0; n < 10*nrelax; n++) {
    gfs_domain_homogeneous_bc (domain, 
			       FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS,
			       levelmin, dp, u);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, 
			      FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, 
			      levelmin,
			      (FttCellTraverseFunc) diffusion_relax, &p);
  }

  /* relax from top to bottom */
  for (p.maxlevel = levelmin + 1; p.maxlevel <= depth; p.maxlevel++) {
    /* get initial guess from coarser grid */ 
    gfs_domain_cell_traverse (domain,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_NON_LEAFS,
			      p.maxlevel - 1,
			      (FttCellTraverseFunc) get_from_above, dp);
    for (n = 0; n < nrelax; n++) {
      gfs_domain_homogeneous_bc (domain, 
				 FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS,
				 p.maxlevel, dp, u);
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, 
				FTT_TRAVERSE_LEVEL | FTT_TRAVERSE_LEAFS, p.maxlevel,
				(FttCellTraverseFunc) diffusion_relax, &p);
    }
  }
  /* correct on leaf cells */
  data[0] = u;
  data[1] = dp;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) correct, data);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, u);
  /* compute new residual on leaf cells */
  gfs_diffusion_residual (domain, u, rhs, dia, res);

  gts_object_destroy (GTS_OBJECT (dp));
}

