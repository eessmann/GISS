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
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "domain.h"

#include "config.h"
#include "advection.h"
#include "source.h"
#include "solid.h"
#ifdef HAVE_MPI
#  include "mpi_boundary.h"
#  include "init.h"
#endif /* HAVE_MPI */

/* GfsDomain: Object */

static void domain_write (GtsObject * o, FILE * fp)
{
  GfsDomain * domain = GFS_DOMAIN (o);

  if (GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class->write) (o, fp);

  fputs (" { ", fp);
  if (domain->rootlevel != 0)
    fprintf (fp, "rootlevel = %u ", domain->rootlevel);
  if (domain->refpos.x != 0.)
    fprintf (fp, "x = %g ", domain->refpos.x);
  if (domain->refpos.y != 0.)
    fprintf (fp, "y = %g ", domain->refpos.y);
  if (domain->refpos.z != 0.)
    fprintf (fp, "z = %g ", domain->refpos.z);
  if (domain->lambda.x != 1.)
    fprintf (fp, "lx = %g ", domain->lambda.x);
  if (domain->lambda.y != 1.)
    fprintf (fp, "ly = %g ", domain->lambda.y);
  if (domain->lambda.z != 1.)
    fprintf (fp, "lz = %g ", domain->lambda.z);
  if (domain->max_depth_write > -2) {
    GfsVariable * v = domain->variables_io;

    if (v != NULL) {
      fprintf (fp, "variables = %s", v->name);
      v = v->next;
      while (v) {
	if (v->name)
	  fprintf (fp, ",%s", v->name);
	v = v->next;
      }
      fputc (' ', fp);
    }
  }
  if (domain->binary != FALSE)
    fprintf (fp, "binary = 1 ");
  fputc ('}', fp);
}

static void domain_read (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain = GFS_DOMAIN (*o);
  GtsFileVariable var[] = {
    {GTS_UINT,   "rootlevel", TRUE},
    {GTS_DOUBLE, "x",         TRUE},
    {GTS_DOUBLE, "y",         TRUE},
    {GTS_DOUBLE, "z",         TRUE},
    {GTS_DOUBLE, "lx",        TRUE},
    {GTS_DOUBLE, "ly",        TRUE},
    {GTS_DOUBLE, "lz",        TRUE},
    {GTS_STRING, "variables", TRUE},
    {GTS_INT,    "binary",    TRUE},
    {GTS_NONE}
  };
  gchar * variables = NULL;

  if (GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  var[0].data = &domain->rootlevel;
  var[1].data = &domain->refpos.x;
  var[2].data = &domain->refpos.y;
  var[3].data = &domain->refpos.z;
  var[4].data = &domain->lambda.x;
  var[5].data = &domain->lambda.y;
  var[6].data = &domain->lambda.z;
  var[7].data = &variables;
  var[8].data = &domain->binary;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR) {
    g_free (variables);
    return;
  }

  if (var[4].set && domain->lambda.x <= 0.) {
    gts_file_error (fp, "lx must be strictly positive");
    return;
  }
  if (var[5].set && domain->lambda.y <= 0.) {
    gts_file_error (fp, "ly must be strictly positive");
    return;
  }
  if (var[6].set && domain->lambda.z <= 0.) {
    gts_file_error (fp, "lz must be strictly positive");
    return;
  }

  if (variables != NULL) {
    gchar * variables1, * s;
    gboolean empty = TRUE;

    variables1 = g_strdup (variables);
    s = strtok (variables1, ",");
    while (s) {
      gfs_domain_add_variable (domain, s);
      empty = FALSE;
      s = strtok (NULL, ",");
    }
    g_free (variables1);

    if (!empty) {
      gchar * error;

      if (domain->variables_io != domain->variables)
	gfs_variable_list_destroy (domain->variables_io);
      domain->variables_io = gfs_variables_from_list (domain->variables, 
						      variables, &error);
      g_assert (domain->variables_io);
    }
    g_free (variables);
  }
}

static void set_ref_pos (GfsBox * box, FttVector * pos)
{
  if (box->id == 1)
    gfs_box_set_pos (box, pos);
}

#ifdef HAVE_MPI
static void removed_list (GfsBox * box, gpointer * data)
{
  GfsDomain * domain = data[0];
  GSList ** removed = data[1];
  guint * np = data[2];

  if (box->pid != domain->pid)
    *removed = g_slist_prepend (*removed, box);
  if (box->pid > *np)
    *np = box->pid;
}

static void mpi_links (GfsBox * box, GfsDomain * domain)
{
  FttDirection d;
  GtsObject * neighbor[FTT_NEIGHBORS];
  gint pid = box->pid;
  gint id = box->id;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOX (box->neighbor[d]) && 
	GFS_BOX (box->neighbor[d])->pid == domain->pid)
      neighbor[d] = box->neighbor[d];
    else
      neighbor[d] = NULL;
  gts_object_destroy (GTS_OBJECT (box));

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (neighbor[d])
#ifndef DUMMY_MPI
      gfs_boundary_mpi_new (gfs_boundary_mpi_class (),
			    GFS_BOX (neighbor[d]), 
			    FTT_OPPOSITE_DIRECTION (d), 
			    pid, id);
#else  /* DUMMY_MPI */
      gfs_boundary_new (GFS_BOUNDARY_CLASS (gfs_boundary_outflow_class ()),
			GFS_BOX (neighbor[d]), 
			FTT_OPPOSITE_DIRECTION (d));
#endif /* DUMMY_MPI */
}
#endif /* HAVE_MPI */

static void domain_post_read (GfsDomain * domain, GtsFile * fp)
{
  gts_graph_foreach_edge (GTS_GRAPH (domain), (GtsFunc) gfs_gedge_link_boxes, NULL);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) set_ref_pos, &domain->refpos);

#ifdef HAVE_MPI
  if (domain->pid >= 0) {
    GSList * removed = NULL;
    guint np = 0;
    gpointer data[3];
    int comm_size;
    
    data[0] = domain;
    data[1] = &removed;
    data[2] = &np;
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) removed_list, data);
    MPI_Comm_size (MPI_COMM_WORLD, &comm_size);
    if (np + 1 != comm_size) {
      g_slist_free (removed);
      gts_file_error (fp, "it would be valid if one or %d PE were used", np + 1);
      return;
    }
    g_slist_foreach (removed, (GFunc) mpi_links, domain);
    g_slist_free (removed);
  }
#endif /* HAVE_MPI */

  gfs_domain_match (domain);
}

static void free_pair (gpointer key, gpointer value)
{
  g_free (key);
  g_free (value);
}

static void domain_destroy (GtsObject * o)
{
  GfsDomain * domain = GFS_DOMAIN (o);

  g_timer_destroy (domain->timer);
  gfs_variable_list_destroy (domain->variables);
  if (domain->variables_io != domain->variables)
    gfs_variable_list_destroy (domain->variables_io);
  g_hash_table_foreach (domain->timers, (GHFunc) free_pair, NULL);
  g_hash_table_destroy (domain->timers);

  (* GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class->destroy) (o);
}

static void domain_class_init (GfsDomainClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = domain_read;
  GTS_OBJECT_CLASS (klass)->write = domain_write;
  GTS_OBJECT_CLASS (klass)->destroy = domain_destroy;

  klass->post_read = domain_post_read;
}

static void domain_init (GfsDomain * domain)
{
#ifdef HAVE_MPI
  int size;

  MPI_Comm_size (MPI_COMM_WORLD, &size);
  if (size > 1)
    MPI_Comm_rank (MPI_COMM_WORLD, &domain->pid);
  else
    domain->pid = -1;
#else /* not HAVE_MPI */
  domain->pid = -1;
#endif /* not HAVE_MPI */

  domain->timer = g_timer_new ();
  domain->timers = g_hash_table_new (g_str_hash, g_str_equal);

  gts_range_init (&domain->size);

  domain->profile_bc = FALSE;

  gts_range_init (&domain->mpi_messages);
  gts_range_init (&domain->mpi_wait);

  domain->rootlevel = 0;
  domain->refpos.x = domain->refpos.y = domain->refpos.z = 0.;
  domain->lambda.x = domain->lambda.y = domain->lambda.z = 1.;
  domain->variables = gfs_variable_list_copy (gfs_centered_variables, 
					      GTS_OBJECT (domain));
  domain->variables_size = sizeof (GfsStateVector);
  domain->variables_io = domain->variables;
  domain->max_depth_write = -1;
}

GfsDomainClass * gfs_domain_class (void)
{
  static GfsDomainClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_domain_info = {
      "GfsDomain",
      sizeof (GfsDomain),
      sizeof (GfsDomainClass),
      (GtsObjectClassInitFunc) domain_class_init,
      (GtsObjectInitFunc) domain_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_wgraph_class ()),
				  &gfs_domain_info);
  }

  return klass;
}

static void box_bc (GfsBox * box, gpointer * datum)
{
  FttTraverseFlags * flags = datum[0];
  gint * max_depth = datum[1];
  GfsVariable * v = datum[2];
  GfsVariable * v1 = datum[4];
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++) 
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, v1);

      if (bc) {
	b->v = v;
	bc->v = v;
	b->type = GFS_BOUNDARY_CENTER_VARIABLE;
	ftt_face_traverse_boundary (b->root, b->d,
				    FTT_PRE_ORDER, *flags, *max_depth,
				    bc->bc, bc);
	bc->v = v1;
	gfs_boundary_send (b);
      }
    }
}

static void direction_face_bc (GtsObject * neighbor,
			       GfsVariable * v)
{
  if (GFS_IS_BOUNDARY (neighbor)) {
    GfsBoundary * b = GFS_BOUNDARY (neighbor);
    GfsBc * bc = gfs_boundary_lookup_bc (b, v);

    if (bc) {
      b->v = v;
      b->type = GFS_BOUNDARY_CENTER_VARIABLE;
      ftt_face_traverse_boundary (b->root, b->d,
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				  bc->face_bc, bc);
      b->type = GFS_BOUNDARY_FACE_VARIABLE;
      gfs_boundary_send (b);
    }
  }
}

static void box_face_bc (GfsBox * box, gpointer * datum)
{
  GfsVariable * v = datum[2];
  FttComponent * c = datum[3];

  if (*c == FTT_XYZ) {
    FttDirection d;
    
    for (d = 0; d < FTT_NEIGHBORS; d++)
      direction_face_bc (box->neighbor[d], v);
  }
  else {
    direction_face_bc (box->neighbor[2*(*c)], v);
    direction_face_bc (box->neighbor[2*(*c) + 1], v);
  }
}

static void box_receive_bc (GfsBox * box, gpointer * datum)
{
  FttTraverseFlags * flags = datum[0];
  gint * max_depth = datum[1];
  FttComponent * c = datum[3];

  if (*c == FTT_XYZ) {
    FttDirection d;
    
    for (d = 0; d < FTT_NEIGHBORS; d++) {
      FttDirection od = FTT_OPPOSITE_DIRECTION (d);

      if (GFS_IS_BOUNDARY (box->neighbor[od]))
	gfs_boundary_receive (GFS_BOUNDARY (box->neighbor[od]), *flags, *max_depth);
    }
  }
  else {
    if (GFS_IS_BOUNDARY (box->neighbor[2*(*c) + 1]))
      gfs_boundary_receive (GFS_BOUNDARY (box->neighbor[2*(*c) + 1]), *flags, *max_depth);
    if (GFS_IS_BOUNDARY (box->neighbor[2*(*c)]))
      gfs_boundary_receive (GFS_BOUNDARY (box->neighbor[2*(*c)]), *flags, *max_depth);
  }
}

static void box_match (GfsBox * box)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * boundary = GFS_BOUNDARY (box->neighbor[d]);

      g_assert (GFS_BOUNDARY_CLASS (box->neighbor[d]->klass)->match);
      boundary->type = GFS_BOUNDARY_MATCH_VARIABLE;
      (* GFS_BOUNDARY_CLASS (box->neighbor[d]->klass)->match) (boundary);
      gfs_boundary_send (boundary);
    }
}

static void box_synchronize (GfsBox * box, FttComponent * c)
{
  if (*c == FTT_XYZ) {
    FttDirection d;
    
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (GFS_IS_BOUNDARY (box->neighbor[d]))
	gfs_boundary_synchronize (GFS_BOUNDARY (box->neighbor[d]));
  }
  else {
    if (GFS_IS_BOUNDARY (box->neighbor[2*(*c)]))
      gfs_boundary_synchronize (GFS_BOUNDARY (box->neighbor[2*(*c)]));
    if (GFS_IS_BOUNDARY (box->neighbor[2*(*c) + 1]))
      gfs_boundary_synchronize (GFS_BOUNDARY (box->neighbor[2*(*c) + 1]));
  }
}

/**
 * gfs_domain_copy_bc:
 * @domain: a #GfsDomain.
 * @flags: the traversal flags.
 * @max_depth: the maximum depth of the traversal.
 * @v: a #GfsVariable.
 * @v1: another #GfsVariable.
 *
 * Apply the boundary conditions of variable @v in @domain to variable @v1.
 */
void gfs_domain_copy_bc (GfsDomain * domain,
			 FttTraverseFlags flags,
			 gint max_depth,
			 GfsVariable * v,
			 GfsVariable * v1)
{
  FttComponent c = FTT_XYZ;
  gpointer datum[5];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);
  g_return_if_fail (v1 != NULL);

  if (domain->profile_bc)
    gfs_domain_timer_start (domain, "bc");

  datum[0] = &flags;
  datum[1] = &max_depth;
  datum[2] = v1;
  datum[3] = &c;
  datum[4] = v;

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_bc, datum);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_receive_bc, datum);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_synchronize, &c);

  if (domain->profile_bc)
    gfs_domain_timer_stop (domain, "bc");
}

/**
 * gfs_domain_bc:
 * @domain: a #GfsDomain.
 * @flags: the traversal flags.
 * @max_depth: the maximum depth of the traversal.
 * @v: a #GfsVariable.
 *
 * Apply the boundary conditions in @domain for variable @v.
 */
void gfs_domain_bc (GfsDomain * domain,
		    FttTraverseFlags flags,
		    gint max_depth,
		    GfsVariable * v)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);

  gfs_domain_copy_bc (domain, flags, max_depth, v, v);
}

static void box_homogeneous_bc (GfsBox * box, gpointer * datum)
{
  FttTraverseFlags * flags = datum[0];
  gint * max_depth = datum[1];
  GfsVariable * ov = datum[2];
  GfsVariable * v = datum[4];
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++) 
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, v);

      if (bc) {
	b->v = ov;
	bc->v = ov;
	ftt_face_traverse_boundary (b->root, b->d,
				    FTT_PRE_ORDER, *flags, *max_depth,
				    bc->homogeneous_bc, bc);
	bc->v = v;
	gfs_boundary_send (b);
      }
    }
}

/**
 * gfs_domain_homogeneous_bc:
 * @domain: a #GfsDomain.
 * @flags: the traversal flags.
 * @max_depth: the maximum depth of the traversal.
 * @ov: a #GfsVariable.
 * @v: a #GfsVariable of which @ov is an homogeneous version.
 *
 * Apply the boundary conditions in @domain for variable @ov using the
 * homogeneous version of the boundary condititons for @v.
 */
void gfs_domain_homogeneous_bc (GfsDomain * domain,
				FttTraverseFlags flags,
				gint max_depth,
				GfsVariable * ov,
				GfsVariable * v)
{
  FttComponent c = FTT_XYZ;
  gpointer datum[5];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (ov != NULL);
  g_return_if_fail (v != NULL);

  if (domain->profile_bc)
    gfs_domain_timer_start (domain, "bc");

  datum[0] = &flags;
  datum[1] = &max_depth;
  datum[2] = ov;
  datum[3] = &c;
  datum[4] = v;
  gts_container_foreach (GTS_CONTAINER (domain),
			 (GtsFunc) box_homogeneous_bc, datum);
  gts_container_foreach (GTS_CONTAINER (domain),
			 (GtsFunc) box_receive_bc, datum);
  gts_container_foreach (GTS_CONTAINER (domain),
			 (GtsFunc) box_synchronize, &c);

  if (domain->profile_bc)
    gfs_domain_timer_stop (domain, "bc");
}

/**
 * gfs_domain_face_bc:
 * @domain: a #GfsDomain.
 * @c: a component.
 * @v: a #GfsVariable.
 *
 * Apply the boundary conditions on the faces of @domain for variable @v.
 */
void gfs_domain_face_bc (GfsDomain * domain,
			 FttComponent c,
			 GfsVariable * v)
{
  FttTraverseFlags flags = FTT_TRAVERSE_LEAFS;
  gint max_depth = -1;
  gpointer datum[4];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (c == FTT_XYZ || (c >= 0 && c < FTT_DIMENSION));
  g_return_if_fail (v != NULL);

  if (domain->profile_bc)
    gfs_domain_timer_start (domain, "face_bc");

  datum[0] = &flags;
  datum[1] = &max_depth;
  datum[2] = v;
  datum[3] = &c;

  gts_container_foreach (GTS_CONTAINER (domain), 
			 (GtsFunc) box_face_bc, datum);
  gts_container_foreach (GTS_CONTAINER (domain), 
			 (GtsFunc) box_receive_bc, datum);
  gts_container_foreach (GTS_CONTAINER (domain),
			 (GtsFunc) box_synchronize, &c);
  
  if (domain->profile_bc)
    gfs_domain_timer_stop (domain, "face_bc");
}

static void box_changed (GfsBox * box, gboolean * changed)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d]))
      *changed |= GFS_BOUNDARY (box->neighbor[d])->changed;
}

static gboolean domain_match (GfsDomain * domain)
{
  FttComponent c = FTT_XYZ;
  FttTraverseFlags flags = FTT_TRAVERSE_LEAFS;
  gint max_depth = -1;
  gboolean changed = FALSE;
  gpointer datum[4];

  datum[0] = &flags;
  datum[1] = &max_depth;
  datum[2] = NULL;
  datum[3] = &c;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_match, NULL);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_receive_bc, datum);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_synchronize, &c);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_changed, &changed);
#ifdef HAVE_MPI
  if (domain->pid >= 0) {
    guint global_changed;

    MPI_Allreduce (&changed, &global_changed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    changed = global_changed;
  }
#endif /* HAVE_MPI */
  return changed;
}

/**
 * gfs_domain_match:
 * @domain: a #GfsDomain.
 *
 * Match the boundaries of @domain.
 */
void gfs_domain_match (GfsDomain * domain)
{
  g_return_if_fail (domain != NULL);

  if (domain->profile_bc)
    gfs_domain_timer_start (domain, "match");

  while (domain_match (domain));

  if (domain->profile_bc)
    gfs_domain_timer_stop (domain, "match");
}

static void dirichlet_bc (FttCell * cell)
{
  cell->flags |= GFS_FLAG_DIRICHLET;
  GFS_STATE (cell)->solid->fv = 0.;
}

static void neumann_bc (FttCell * cell)
{
  cell->flags &= ~GFS_FLAG_DIRICHLET;
  GFS_STATE (cell)->solid->fv = 0.;
}

/**
 * gfs_domain_surface_bc:
 * @domain: a #GfsDomain.
 * @v: a #GfsVariable.
 *
 * Apply boundary conditions for variable @v on embedded surfaces. 
 */
void gfs_domain_surface_bc (GfsDomain * domain,
			    GfsVariable * v)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);

  if (v->surface_bc)
    gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL,
      (FttCellTraverseFunc) GFS_SURFACE_GENERIC_BC_CLASS (GTS_OBJECT (v->surface_bc)->klass)->bc, 
			       v->surface_bc);
  else {
    if (GFS_VELOCITY_COMPONENT (v->i) < FTT_DIMENSION)
      gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL,
				 (FttCellTraverseFunc) dirichlet_bc, NULL);
    else
      gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL,
				 (FttCellTraverseFunc) neumann_bc, NULL);
  }
}

static void box_traverse (GfsBox * box, gpointer * datum)
{
  FttTraverseType * order = datum[0];
  FttTraverseFlags * flags = datum[1];
  gint * max_depth = datum[2];
  FttCellTraverseFunc func = (FttCellTraverseFunc) datum[3];
  gpointer data = datum[4];

  ftt_cell_traverse (box->root, *order, *flags, *max_depth, func, data);
}

/**
 * gfs_domain_cell_traverse:
 * @domain: a #GfsDomain.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Traverses the cell trees of @domain. Calls the given function for
 * each cell visited.  
 */
void gfs_domain_cell_traverse (GfsDomain * domain,
			       FttTraverseType order,
			       FttTraverseFlags flags,
			       gint max_depth,
			       FttCellTraverseFunc func,
			       gpointer data)
{
  gpointer datum[5];

  datum[0] = &order;
  datum[1] = &flags;
  datum[2] = &max_depth;
  datum[3] = func;
  datum[4] = data;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (func != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), 
			 (GtsFunc) box_traverse, datum);
}

static void box_traverse_box (GfsBox * box, gpointer * datum)
{
  FttTraverseType * order = datum[0];
  FttTraverseFlags * flags = datum[1];
  gint * max_depth = datum[2];
  FttCellTraverseFunc func = (FttCellTraverseFunc) datum[3];
  gpointer data = datum[4];
  GtsBBox * bb = datum[5];

  ftt_cell_traverse_box (box->root, bb, 
			 *order, *flags, *max_depth, func, data);
}

/**
 * gfs_domain_cell_traverse_box:
 * @domain: a #GfsDomain.
 * @box: a #GtsBBox.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Traverses the cell trees of @domain. Calls the given function for
 * each cell visited. Only the cells overlapping with @box are visited.
 */
void gfs_domain_cell_traverse_box (GfsDomain * domain,
				   GtsBBox * box,
				   FttTraverseType order,
				   FttTraverseFlags flags,
				   gint max_depth,
				   FttCellTraverseFunc func,
				   gpointer data)
{
  gpointer datum[6];

  datum[0] = &order;
  datum[1] = &flags;
  datum[2] = &max_depth;
  datum[3] = func;
  datum[4] = data;
  datum[5] = box;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (box != NULL);
  g_return_if_fail (func != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), 
			 (GtsFunc) box_traverse_box, datum);
}

static void traverse_mixed (GfsBox * box, gpointer * datum)
{
  FttCellTraverseFunc func = (FttCellTraverseFunc) datum[0];
  gpointer data = datum[1];
  FttTraverseType * order = datum[2];
  FttTraverseFlags * flags = datum[3];

  gfs_cell_traverse_mixed (box->root, *order, *flags, func, data);
}

/**
 * gfs_domain_traverse_mixed:
 * @domain: a #GfsDomain.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Calls @func for each mixed cell of @domain.
 */
void gfs_domain_traverse_mixed (GfsDomain * domain,
				FttTraverseType order,
				FttTraverseFlags flags,
				FttCellTraverseFunc func,
				gpointer data)
{
  gpointer datum[4];

  datum[0] = func;
  datum[1] = data;
  datum[2] = &order;
  datum[3] = &flags;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (func != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) traverse_mixed, datum);
}

static void traverse_cut (GfsBox * box, gpointer * datum)
{
  FttCellTraverseCutFunc func = (FttCellTraverseCutFunc) datum[0];
  gpointer data = datum[1];
  FttTraverseType * order = datum[2];
  FttTraverseFlags * flags = datum[3];
  GtsSurface * s = datum[4];

  gfs_cell_traverse_cut (box->root, s, *order, *flags, func, data);
}

/**
 * gfs_domain_traverse_cut:
 * @domain: a #GfsDomain.
 * @s: a #GtsSurface.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Calls @func for each cell of @domain cut by @s.
 */
void gfs_domain_traverse_cut (GfsDomain * domain,
			      GtsSurface * s,
			      FttTraverseType order,
			      FttTraverseFlags flags,
			      FttCellTraverseCutFunc func,
			      gpointer data)
{
  gpointer datum[5];

  datum[0] = func;
  datum[1] = data;
  datum[2] = &order;
  datum[3] = &flags;
  datum[4] = s;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (func != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) traverse_cut, datum);
}

static void box_depth (GfsBox * box, guint * depth)
{
  guint d = ftt_cell_depth (box->root);

  if (d > *depth)
    *depth = d;
}

/**
 * gfs_domain_depth:
 * @domain: a #GfsDomain.
 *
 * Returns: the maximum depth of the cell trees of @domain. This
 * function is global i.e. it returns the maximum depth over all the
 * processes (for parallel execution).
 */
guint gfs_domain_depth (GfsDomain * domain)
{
  guint depth = 0;

  g_return_val_if_fail (domain != NULL, 0);

  gts_container_foreach (GTS_CONTAINER (domain),
			 (GtsFunc) box_depth, &depth);
#ifdef HAVE_MPI
  if (domain->pid >= 0) {
    guint global_depth;

    MPI_Allreduce (&depth, &global_depth, 1, MPI_UNSIGNED, MPI_MAX, 
		   MPI_COMM_WORLD);
    depth = global_depth;
  }
#endif /* HAVE_MPI */
  return depth;
}

#include "ftt_internal.c"

/**
 * gfs_domain_face_traverse:
 * @domain: a #GfsDomain.
 * @c: only the faces orthogonal to this component will be traversed - one of
 * %FTT_X, %FTT_Y, (%FTT_Z), %FTT_XYZ.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children and faces are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCellFace.
 * @data: user data to pass to @func.
 *
 * Traverses a @domain. Calls the given function for each face
 * of the cells of the domain.
 *
 * If %FTT_TRAVERSE_BOUNDARY_FACES is not set in @flags, only
 * "double-sided" faces are traversed i.e. the @neighbor field of the
 * face is never %NULL.  
 */
void gfs_domain_face_traverse (GfsDomain * domain,
			       FttComponent c,
			       FttTraverseType order,
			       FttTraverseFlags flags,
			       gint max_depth,
			       FttFaceTraverseFunc func,
			       gpointer data)
{
  FttDirection d;
  gpointer datum[6];
  gboolean check = FALSE;
  gboolean boundary_faces;
  
  g_return_if_fail (domain != NULL);
  g_return_if_fail (c >= FTT_X && c <= FTT_XYZ);
  g_return_if_fail (func != NULL);

  boundary_faces = ((flags & FTT_TRAVERSE_BOUNDARY_FACES) != 0);
  datum[1] = &max_depth;
  datum[2] = func;
  datum[3] = data;
  datum[4] = &check;
  datum[5] = &boundary_faces;
  if (c == FTT_XYZ) {
    if (boundary_faces) {
      check = TRUE;
      gfs_domain_cell_traverse (domain, order, flags, max_depth, 
	  (FttCellTraverseFunc) traverse_all_faces, 
				datum);
    }
    else {
      gfs_domain_cell_traverse (domain, order, flags, max_depth, 
	  (FttCellTraverseFunc) traverse_all_direct_faces, 
				datum);
      datum[0] = &d;
      check = TRUE;
      for (d = 1; d < FTT_NEIGHBORS; d += 2)
	gfs_domain_cell_traverse_boundary (domain, 
					   d, order, flags, max_depth, 
					   (FttCellTraverseFunc) traverse_face, datum);
    }
  }
  else if (c == FTT_XY) {
    gfs_domain_face_traverse (domain, FTT_X, order, flags, max_depth, func, data);
    gfs_domain_face_traverse (domain, FTT_Y, order, flags, max_depth, func, data);
  }
  else {
    if (boundary_faces) {
      check = TRUE;
      datum[0] = &c;
      gfs_domain_cell_traverse (domain, order, flags, max_depth, 
				(FttCellTraverseFunc) traverse_face_component,
				datum);
    }
    else {
      d = 2*c;
      datum[0] = &d;
      gfs_domain_cell_traverse (domain, order, flags, max_depth, 
				(FttCellTraverseFunc) traverse_face_direction, 
				datum);
      d = 2*c + 1;
      check = TRUE;
      gfs_domain_cell_traverse_boundary (domain, d, order, flags, max_depth, 
					 (FttCellTraverseFunc) traverse_face, datum);
    }
  }
  gfs_domain_cell_traverse (domain, order, flags, max_depth, 
			    (FttCellTraverseFunc) reset_flag, NULL);
}

static void cell_traverse_boundary (GfsBox * box, gpointer * datum)
{
  FttDirection * d = datum[0];

  if (!GFS_IS_BOX (box->neighbor[*d])) {
    FttTraverseType * order = datum[1];
    FttTraverseFlags * flags = datum[2];
    gint * max_depth = datum[3];
    FttCellTraverseFunc func = (FttCellTraverseFunc) datum[4];
    gpointer data = datum[5];

    ftt_cell_traverse_boundary (box->root, 
				*d, *order, *flags, *max_depth, func, data);
  }
}

/**
 * gfs_domain_cell_traverse_boundary:
 * @domain: a #GfsDomain.
 * @d: the direction of the boundary to traverse.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Traverses the boundary of a domain in direction @d. Calls the given
 * function for each cell visited.  
 */
void gfs_domain_cell_traverse_boundary (GfsDomain * domain,
					FttDirection d,
					FttTraverseType order,
					FttTraverseFlags flags,
					gint max_depth,
					FttCellTraverseFunc func,
					gpointer data)
{
  gpointer datum[6];
  
  datum[0] = &d;
  datum[1] = &order;
  datum[2] = &flags;
  datum[3] = &max_depth;
  datum[4] = func;
  datum[5] = data;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (d < FTT_NEIGHBORS);
  g_return_if_fail (func != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) cell_traverse_boundary, datum);
}

static void add_stats (const FttCell * cell, gpointer * data)
{
  GtsRange * s = data[0];
  GfsVariable * v = data[1];

  gts_range_add_value (s, GFS_VARIABLE (cell, v->i));
}

#ifdef HAVE_MPI
static void range_reduce (void * i, void * o, 
			  int * len,
			  MPI_Datatype * type)
{
  gdouble * in = (gdouble *) i;
  gdouble * inout = (gdouble *) o;
  g_assert (*len == 5);
  
  if (in[0] < inout[0]) /* min */
    inout[0] = in[0];
  if (in[1] > inout[1]) /* max */
    inout[1] = in[1];
  inout[2] += in[2];    /* sum */
  inout[3] += in[3];    /* sum2 */
  inout[4] += in[4];    /* n */
}

static void domain_range_reduce (GfsDomain * domain, GtsRange * s)
{
  if (domain->pid >= 0) {
    double in[5];
    double out[5] = { G_MAXDOUBLE, - G_MAXDOUBLE, 0., 0., 0. };
    MPI_Op op;
    
    MPI_Op_create (range_reduce, TRUE, &op);
    in[0] = s->min; in[1] = s->max; in[2] = s->sum; in[3] = s->sum2;
    in[4] = s->n;
    MPI_Allreduce (in, out, 5, MPI_DOUBLE, op, MPI_COMM_WORLD);
    MPI_Op_free (&op);
    s->min = out[0]; s->max = out[1]; s->sum = out[2]; s->sum2 = out[3];
    s->n = out[4];
  }
}
#endif /* HAVE_MPI */

/**
 * gfs_domain_stats_variable:
 * @domain: the domain to obtain statistics from.
 * @v: a #GfsVariable.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 *
 * Traverses the domain defined by @domain using gfs_domain_cell_traverse()
 * and gathers statistics about variable @v.
 *
 * Returns: a #GtsRange containing the statistics about @v.
 */
GtsRange gfs_domain_stats_variable (GfsDomain * domain,
				    GfsVariable * v,
				    FttTraverseFlags flags,
				    gint max_depth)
{
  GtsRange s;
  gpointer data[2];

  g_return_val_if_fail (domain != NULL, s);
  g_return_val_if_fail (v != NULL, s);

  gts_range_init (&s);
  data[0] = &s;
  data[1] = v;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			   (FttCellTraverseFunc) add_stats, data);
#ifdef HAVE_MPI
  domain_range_reduce (domain, &s);
#endif /* HAVE_MPI */
  gts_range_update (&s);

  return s;
}

static void add_stats_solid (FttCell * cell, GtsRange * s)
{
  gts_range_add_value (s, GFS_STATE (cell)->solid->a);
}

/**
 * gfs_domain_stats_solid:
 * @domain: the domain to obtain statistics from.
 *
 * Traverses the domain defined by @domain using gfs_domain_traverse_mixed()
 * and gathers statistics about the solid volume fraction in mixed cells.
 *
 * Returns: statistics about the solid volume fraction @a in mixed cells.
 */
GtsRange gfs_domain_stats_solid (GfsDomain * domain)
{
  GtsRange s;

  g_return_val_if_fail (domain != NULL, s);

  gts_range_init (&s);
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			    (FttCellTraverseFunc) add_stats_solid, &s);
#ifdef HAVE_MPI
  domain_range_reduce (domain, &s);
#endif /* HAVE_MPI */
  gts_range_update (&s);

  return s;
}

static void add_stats_merged (GSList * m, gpointer * data)
{
  GtsRange * solid =  data[0];
  GtsRange * number = data[1];
  gdouble a = 0.;
  guint n = 0;

  while (m) {
    FttCell * c = m->data;

    a += GFS_IS_MIXED (c) ? GFS_STATE (c)->solid->a : 1.;
    n++;
    m = m->next;
  }
  if (n > 1 || a < 1.)
    gts_range_add_value (solid, a);
  if (n > 1)
    gts_range_add_value (number, n);
}

/**
 * gfs_domain_stats_merged:
 * @domain: the domain to obtain statistics from.
 * @solid: #GtsRange in which to return stats for the total solid
 * volume fraction of merged cells. 
 * @number: #GtsRange in which to return stats for the number of cells
 * used per merged cell.
 *
 * Traverses the domain defined by @domain using
 * gfs_domain_traverse_merged() and gathers statistics about the total
 * solid volume fraction of merged cells and the number of cells used
 * per merged cell.
 */
void gfs_domain_stats_merged (GfsDomain * domain,
			     GtsRange * solid,
			     GtsRange * number)
{
  gpointer data[2];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (solid != NULL);
  g_return_if_fail (number != NULL);

  gts_range_init (solid);
  gts_range_init (number);
  data[0] = solid;
  data[1] = number;
  gfs_domain_traverse_merged (domain,
			     (GfsMergedTraverseFunc) add_stats_merged, data);
#ifdef HAVE_MPI
  domain_range_reduce (domain, solid);
  domain_range_reduce (domain, number);
#endif /* HAVE_MPI */
  gts_range_update (solid);
  gts_range_update (number);
}

static void cell_count (FttCell * cell, guint * count)
{
  (*count)++;
}

#ifdef HAVE_MPI
static void boundary_size (GfsBox * box, guint * count)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY_MPI (box->neighbor[d]))
      ftt_cell_traverse (GFS_BOUNDARY (box->neighbor[d])->root,
			 FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			 (FttCellTraverseFunc) cell_count, count);
}
#endif /* HAVE_MPI */

/**
 * gfs_domain_stats_balance:
 * @domain: the domain to obtain statistics from.
 * @size: #GtsRange in which to return stats for the total size of the domain.
 * @boundary: #GtsRange in which to return stats for the size of the parallel 
 * boundaries of the domain.
 * @mpiwait:  #GtsRange in which to return stats for the average time spend
 * waiting for MPI calls in each PE.
 *
 * Gathers statistics about the sizes of the domains, their parallel
 * boundaries and the execution time on each PE.  
 */
void gfs_domain_stats_balance (GfsDomain * domain,
			       GtsRange * size,
			       GtsRange * boundary,
			       GtsRange * mpiwait)
{
  guint count;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (size != NULL);
  g_return_if_fail (boundary != NULL);
  g_return_if_fail (mpiwait != NULL);

  gts_range_init (size);
  gts_range_init (boundary);
  gts_range_init (mpiwait);
  count = 0;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			   (FttCellTraverseFunc) cell_count, &count);
  gts_range_add_value (size, count);
  if (domain->timestep.n > 0)
    gts_range_add_value (mpiwait, domain->mpi_wait.sum/domain->timestep.n);
#ifdef HAVE_MPI
  count = 0;
  gts_container_foreach (GTS_CONTAINER (domain),
			 (GtsFunc) boundary_size, &count);
  gts_range_add_value (boundary, count);
  domain_range_reduce (domain, size);
  domain_range_reduce (domain, boundary);
  domain_range_reduce (domain, mpiwait);
#endif /* HAVE_MPI */
  gts_range_update (size);
  gts_range_update (boundary);
  gts_range_update (mpiwait);
}

static void add_norm (const FttCell * cell, gpointer * data)
{
  GfsNorm * n = data[0];
  GfsVariable * v = data[1];

  gfs_norm_add (n, GFS_VARIABLE (cell, v->i),
		ftt_cell_volume (cell)*(GFS_IS_MIXED (cell) ?
					GFS_STATE (cell)->solid->a : 1.));
}

#ifdef HAVE_MPI
static void norm_reduce (void * i, void * o, 
			 int * len,
			 MPI_Datatype * type)
{
  gdouble * in = (gdouble *) i;
  gdouble * inout = (gdouble *) o;
  g_assert (*len == 5);
  
  inout[0] += in[0];    /* bias */
  inout[1] += in[1];    /* first */
  inout[2] += in[2];    /* second */
  if (in[3] > inout[3]) /* infty */
    inout[3] = in[3];    
  inout[4] += in[4];    /* w */
}

static void domain_norm_reduce (GfsDomain * domain, GfsNorm * n)
{
  if (domain->pid >= 0) {
    double in[5];
    double out[5] = { 0., 0., 0., - G_MAXDOUBLE, 0. };
    MPI_Op op;

    MPI_Op_create (norm_reduce, TRUE, &op);
    in[0] = n->bias; in[1] = n->first; in[2] = n->second; in[3] = n->infty;
    in[4] = n->w;
    MPI_Allreduce (in, out, 5, MPI_DOUBLE, op, MPI_COMM_WORLD);
    MPI_Op_free (&op);
    n->bias = out[0]; n->first = out[1]; n->second = out[2]; n->infty = out[3];
    n->w = out[4];
  }
}
#endif /* HAVE_MPI */

/**
 * gfs_domain_norm_variable:
 * @domain: the domain to obtain norm from.
 * @v: a #GfsVariable.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 *
 * Traverses the domain defined by @domain using gfs_domain_cell_traverse()
 * and gathers norm statistics about variable @v.
 *
 * The norm is weighted by the volume of each cell.
 *
 * Returns: a #GfsNorm containing the norm statistics about @v.
 */
GfsNorm gfs_domain_norm_variable (GfsDomain * domain,
				  GfsVariable * v,
				  FttTraverseFlags flags,
				  gint max_depth)
{
  GfsNorm n;
  gpointer data[2];

  g_return_val_if_fail (domain != NULL, n);
  g_return_val_if_fail (v != NULL, n);
  
  gfs_norm_init (&n);
  data[0] = &n;
  data[1] = v;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			   (FttCellTraverseFunc) add_norm, data);
#ifdef HAVE_MPI
  domain_norm_reduce (domain, &n);
#endif /* HAVE_MPI */
  gfs_norm_update (&n);

  return n;
}

static void add_norm_residual (const FttCell * cell, GfsNorm * n)
{
  gdouble size = ftt_cell_size (cell);

  gfs_norm_add (n, GFS_STATE (cell)->res/(size*size), 1.);
}

/**
 * gfs_domain_norm_residual:
 * @domain: the domain to obtain the norm from.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 * @dt: the time step.
 *
 * Traverses the domain defined by @domain using gfs_domain_cell_traverse()
 * and gathers norm statistics about the volume weighted relative residual
 * (i.e. the sum of the residual over the volume defined by each cell
 * divided by the total volume of the cell).
 *
 * Returns: a #GfsNorm containing the norm statistics about the volume
 * weighted relative residual.  
 */
GfsNorm gfs_domain_norm_residual (GfsDomain * domain,
				  FttTraverseFlags flags,
				  gint max_depth,
				  gdouble dt)
{
  GfsNorm n;

  g_return_val_if_fail (domain != NULL, n);
  
  gfs_norm_init (&n);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			   (FttCellTraverseFunc) add_norm_residual, &n);
#ifdef HAVE_MPI
  domain_norm_reduce (domain, &n);
#endif /* HAVE_MPI */
  gfs_norm_update (&n);
  
  n.bias *= dt;
  n.first *= dt;
  n.second *= dt;
  n.infty *= dt;
  return n;
}

static void add_norm_velocity (const FttCell * cell, GfsNorm * n)
{
  FttComponent c;
  gdouble unorm = 0.;
  
  for (c = 0; c < FTT_DIMENSION; c++) {
    gdouble uc = GFS_VARIABLE (cell, GFS_VELOCITY_INDEX (c));

    unorm += uc*uc;
  }
  gfs_norm_add (n, sqrt (unorm), 
		ftt_cell_volume (cell)*(GFS_IS_MIXED (cell) ? 
					GFS_STATE (cell)->solid->a : 1.));
}

/**
 * gfs_domain_norm_velocity:
 * @domain: the domain to obtain the norm from.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 *
 * Traverses the domain defined by @domain using gfs_domain_cell_traverse()
 * and gathers norm statistics about velocity.
 *
 * Returns: a #GfsNorm containing the norm statistics about the velocity.
 */
GfsNorm gfs_domain_norm_velocity (GfsDomain * domain,
				  FttTraverseFlags flags,
				  gint max_depth)
{
  GfsNorm n;

  g_return_val_if_fail (domain != NULL, n);
  
  gfs_norm_init (&n);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			   (FttCellTraverseFunc) add_norm_velocity, &n);
#ifdef HAVE_MPI
  domain_norm_reduce (domain, &n);
#endif /* HAVE_MPI */
  gfs_norm_update (&n);

  return n;
}

/**
 * gfs_domain_read:
 * @fp: a #GtsFile.
 *
 * Reads the graph nodes (#GfsBox) and edges and the
 * corresponding boundaries (#GfsBoundaryMpi if necessary) defined in
 * @fp.
 *
 * Returns: the #GfsDomain or %NULL if an error occured, in which case
 * the corresponding @fp fields (@pos and @error) are set.
 */
GfsDomain * gfs_domain_read (GtsFile * fp)
{
  GfsDomain * domain;

  g_return_val_if_fail (fp != NULL, NULL);
						 
  if (!(domain = GFS_DOMAIN (gts_graph_read (fp))))
    return NULL;

  (* GFS_DOMAIN_CLASS (GTS_OBJECT (domain)->klass)->post_read) (domain, fp);

  return domain;
}

static void box_split (GfsBox * box, gpointer * data)
{
  GSList ** boxlist = data[0];
  guint * bid = data[1];
  gboolean * one_box_per_pe = data[2];
  gint * pid = data[3];
  guint refid = FTT_DIMENSION == 2 ? 2 : 6;
  FttCellChildren child;
  FttDirection d;
  guint i;
  GfsDomain * domain = gfs_box_domain (box);

  *boxlist = g_slist_prepend (*boxlist, box);

  if (FTT_CELL_IS_LEAF (box->root))
    ftt_cell_refine_single (box->root, (FttCellInitFunc) gfs_cell_init, domain);

  ftt_cell_children (box->root, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
      GfsBox * newbox = GFS_BOX (gts_object_new (GTS_OBJECT (box)->klass));

      GTS_OBJECT (newbox)->reserved = domain;
      if (*one_box_per_pe)
	newbox->pid = (*pid)++;
      else
	newbox->pid = box->pid;
      if (box->id == 1 && i == refid)
	newbox->id = 1;
      else
	newbox->id = (*bid)++;

      GFS_DOUBLE_TO_POINTER (GFS_STATE (child.c[i])->div) = newbox;

      if (FTT_CELL_IS_LEAF (child.c[i]))
	ftt_cell_refine_single (child.c[i], (FttCellInitFunc) gfs_cell_init, domain);
    }

#if FTT_2D3
  g_assert_not_implemented ();
#endif
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * boundary = GFS_BOUNDARY (box->neighbor[d]);

      ftt_cell_children (boundary->root, &child);
      for (i = 0; i < FTT_CELLS; i++)
	if (child.c[i] && FTT_CELL_IS_LEAF (child.c[i]))
	  ftt_cell_refine_single (child.c[i], (FttCellInitFunc) gfs_cell_init, domain);
      ftt_cell_destroy_root (boundary->root, &child, (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
      boundary->root = NULL;

      ftt_cell_children_direction (box->root, d, &child);
      for (i = 0; i < FTT_CELLS/2; i++)
	if (child.c[i]) {
	  FttCell * neighbor = ftt_cell_neighbor (child.c[i], d);
	  GfsBox * newbox = GFS_DOUBLE_TO_POINTER (GFS_STATE (child.c[i])->div);
	  GfsBoundaryClass * klass = GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass);
	  GfsBoundary * newboundary = gfs_boundary_new (klass, newbox, d);
	  gchar fname[] = "/tmp/XXXXXX";
	  gint fd = mkstemp (fname);
	  FILE * fp = fdopen (fd, "w");
	  GtsFile * gfp;

	  (* GTS_OBJECT_CLASS (klass)->write) (GTS_OBJECT (boundary), fp);
	  fclose (fp);
	  close (fd);
	  fp = fopen (fname, "r");
	  unlink (fname);
	  gfp = gts_file_new (fp);
	  (* GTS_OBJECT_CLASS (klass)->read) ((GtsObject **) &newboundary, gfp);
	  g_assert (gfp->type != GTS_ERROR);
	  gts_file_destroy (gfp);
	  fclose (fp);

	  g_assert (neighbor);
	  newboundary->root = neighbor;
	}
      gts_object_destroy (GTS_OBJECT (boundary));
    }
}

static void box_link (GfsBox * box, GfsDomain * domain)
{
  FttCellChildren child;
  guint i;

  ftt_cell_children (box->root, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
       GfsBox * newbox = GFS_DOUBLE_TO_POINTER (GFS_STATE (child.c[i])->div);
       FttDirection d;
       
       g_assert (newbox);
       gts_container_add (GTS_CONTAINER (domain), GTS_CONTAINEE (newbox));
       for (d = 0; d < FTT_NEIGHBORS; d++)
	 if (newbox->neighbor[d] == NULL) {
	   FttCell * neighbor = ftt_cell_neighbor (child.c[i], d);

	   if (neighbor) {
	     GfsBox * newbox1 = GFS_DOUBLE_TO_POINTER (GFS_STATE (neighbor)->div);
	     FttDirection od = FTT_OPPOSITE_DIRECTION (d);
	     GfsGEdge * edge;

	     g_assert (newbox1);
	     newbox->neighbor[d] = GTS_OBJECT (newbox1);
	     g_assert (newbox1->neighbor[od] == NULL);
	     newbox1->neighbor[od] = GTS_OBJECT (newbox);
	     edge = GFS_GEDGE (gts_gedge_new (GTS_GRAPH (domain)->edge_class,
					      GTS_GNODE (newbox), 
					      GTS_GNODE (newbox1)));
	     edge->d = d;
	   }
	 }
    }
}

static void box_destroy (GfsBox * box)
{
  GfsBox * newbox[FTT_CELLS];
  FttCellChildren child;
  guint i;

  ftt_cell_children (box->root, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i])
      newbox[i] = GFS_DOUBLE_TO_POINTER (GFS_STATE (child.c[i])->div);
    else
      newbox[i] = NULL;

  ftt_cell_destroy_root (box->root, &child, (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
  box->root = NULL;
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i])
      newbox[i]->root = child.c[i];

  gts_object_destroy (GTS_OBJECT (box));
}

static void get_ref_pos (GfsBox * box, FttVector * pos)
{
  if (box->id == 1)
    ftt_cell_pos (box->root, pos);
}

/**
 * gfs_domain_split:
 * @domain: a #GfsDomain.
 * @one_box_per_pe: if %TRUE each new box created is assigned to a
 * different process, otherwise the newly created box inherits the pid
 * of its parent.
 *
 * Splits each box of @domain into its (4 in 2D, 8 in 3D)
 * children. The corresponding newly created boxes are added to the
 * graph and the parent boxes are destroyed.
 */
void gfs_domain_split (GfsDomain * domain, gboolean one_box_per_pe)
{
  GSList * list = NULL;
  guint bid = 2;
  gint pid = 0;
  gpointer data[4];

  g_return_if_fail (domain != NULL);

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, 1,
  			   (FttCellTraverseFunc) gfs_cell_reset, gfs_div);
  data[0] = &list;
  data[1] = &bid;
  data[2] = &one_box_per_pe;
  data[3] = &pid;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_split, data);
  g_slist_foreach (list, (GFunc) box_link, domain);
  g_slist_foreach (list, (GFunc) box_destroy, NULL);
  g_slist_free (list);

  gfs_domain_match (domain);
  domain->rootlevel++;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) get_ref_pos, &domain->refpos);
}

static void box_locate (GfsBox * box, gpointer * data)
{
  FttVector * target = data[0];
  gint * max_depth = data[1];
  FttCell ** cell = data[2];

  if (*cell == NULL)
    *cell = ftt_cell_locate (box->root, *target, *max_depth);
}

/**
 * gfs_domain_locate:
 * @domain: a #GfsDomain.
 * @target: position of the point to look for.
 * @max_depth: maximum depth to consider (-1 means no restriction).
 *
 * Locates the cell of @domain containing @target. This is done
 * efficiently in log(n) operations by using the topology of the cell
 * trees.
 *
 * Returns: a #FttCell of @domain containing (boundary included) the
 * point defined by @target or %NULL if @target is not contained in
 * any cell of @domain.  
 */
FttCell * gfs_domain_locate (GfsDomain * domain,
			     FttVector target,
			     gint max_depth)
{
  FttCell * cell = NULL;
  gpointer data[3];

  g_return_val_if_fail (domain != NULL, NULL);

  data[0] = &target;
  data[1] = &max_depth;
  data[2] = &cell;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_locate, data);

  return cell;
}

/**
 * gfs_domain_advect_point:
 * @domain: a #GfsDomain.
 * @p: a #GtsPoint.
 * @dt: the time step.
 *
 * Updates the coordinates of point @p at time t + @dt using the
 * velocity field defined by @domain.
 *
 * If @p is not contained within @domain, the coordinates are unchanged.
 */
void gfs_domain_advect_point (GfsDomain * domain, 
			      GtsPoint * p,
			      gdouble dt)
{
  FttCell * cell;
  FttVector p0, p1;
  FttComponent c;
  GfsVariable * v1, * v;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (p != NULL);

  p0.x = p1.x = p->x; 
  p0.y = p1.y = p->y;
  p0.z = p1.z = p->z;
  cell = gfs_domain_locate (domain, p0, -1);
  if (cell == NULL)
    return;
  v1 = v = gfs_variable_from_name (domain->variables, "U");
  for (c = 0; c < FTT_DIMENSION; c++, v = v->next)
    (&p1.x)[c] += dt*gfs_interpolate (cell, p0, v)/2.;
  cell = gfs_domain_locate (domain, p1, -1);
  if (cell == NULL)
    return;
  v = v1;
  for (c = 0; c < FTT_DIMENSION; c++, v = v->next)
    (&p->x)[c] += dt*gfs_interpolate (cell, p1, v);
}

static void count (FttCell * cell, guint * n)
{
  (*n)++;
}

/**
 * gfs_domain_size:
 * @domain: a #GfsDomain.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 *
 * Returns: the number of cells of @domain traversed using @flags and
 * @max_depth.
 */
guint gfs_domain_size (GfsDomain * domain,
		       FttTraverseFlags flags,
		       gint max_depth)
{
  guint n = 0;

  g_return_val_if_fail (domain != NULL, 0);
  
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			   (FttCellTraverseFunc) count, &n);
#ifdef HAVE_MPI
  if (domain->pid >= 0) {
    guint sn;

    MPI_Allreduce (&n, &sn, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    n = sn;
  }
#endif /* HAVE_MPI */
  return n;
}

static void minimum_cfl (FttCell * cell, gpointer * data)
{
  gdouble * cfl = data[0];
  GfsVariable * v = data[1];
  gdouble size = ftt_cell_size (cell);
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++, v = v->next) {
    if (GFS_VARIABLE (cell, v->i) != 0.) {
      gdouble cflu = size/fabs (GFS_VARIABLE (cell, v->i));

      if (cflu*cflu < *cfl)
	*cfl = cflu*cflu;
    }
    if (v->sources) {
      gdouble g = gfs_variable_mac_source (v, cell);

      if (g != 0.) {
	gdouble cflg = 2.*size/fabs (g);

	if (cflg < *cfl)
	  *cfl = cflg;
      }
    }
  }
}

/**
 * gfs_domain_cfl:
 * @domain: a #GfsDomain.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 *
 * Returns: the minimum over the cells of @domain (traversed using
 * @flags and @max_depth) of the time scale defined by the size of the
 * cell and the norm of either the local velocity or the local
 * acceleration.
 */
gdouble gfs_domain_cfl (GfsDomain * domain,
			FttTraverseFlags flags,
			gint max_depth)
{
  gdouble cfl = 1.;
  gpointer data[2];

  g_return_val_if_fail (domain != NULL, 0.);

  data[0] = &cfl;
  data[1] = gfs_variable_from_name (domain->variables, "U");
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			    (FttCellTraverseFunc) minimum_cfl, data);
#ifdef HAVE_MPI
  if (domain->pid >= 0) {
    gdouble gcfl;

    MPI_Allreduce (&cfl, &gcfl, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    cfl = gcfl;
  }
#endif /* HAVE_MPI */
  return sqrt (cfl);
}

/**
 * gfs_cell_init:
 * @cell: a #FttCell.
 * @domain: a #GfsDomain containing @cell.
 *
 * Allocates the memory for fluid state data associated to @cell.
 */
void gfs_cell_init (FttCell * cell, GfsDomain * domain)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (cell->data == NULL);
  g_return_if_fail (domain != NULL);

  cell->data = g_malloc0 (domain->variables_size);
}

/**
 * gfs_cell_copy:
 * @from: a #FttCell to copy attributes from.
 * @to: a #FttCell to copy attributes to.
 * @domain: the #GfsDomain containing @from.
 *
 * Copies the attributes of the fluid cell @from to the fluid cell @to.
 */
void gfs_cell_copy (const FttCell * from, 
		    FttCell * to,
		    GfsDomain * domain)
{
  GfsSolidVector * solid;
  GfsStateVector * froms, * tos;

  g_return_if_fail (from != NULL);
  g_return_if_fail (to != NULL);
  g_return_if_fail (from != to);  
  g_return_if_fail (domain != NULL);

  froms = GFS_STATE (from);
  tos = GFS_STATE (to);
  if (froms != NULL) {
    if (tos == NULL) {
      gfs_cell_init (to, domain);
      tos = GFS_STATE (to);
    }
    solid = tos->solid;
    memcpy (to->data, from->data, domain->variables_size);
    if (froms->solid == NULL) {
      if (solid)
	g_free (solid);
    }
    else {
      tos->solid = solid;
      *solid = *(froms->solid);
    }
  }
  else if (tos != NULL)
    gfs_cell_cleanup (to);
}

/**
 * gfs_cell_write:
 * @cell: a #FttCell.
 * @fp: a file pointer.
 * @variables: the #GfsVariable to be written.
 *
 * Writes in @fp the fluid data associated with @cell and described by
 * @variables. This function is generally used in association with
 * ftt_cell_write().  
 */
void gfs_cell_write (const FttCell * cell, FILE * fp,
		     GfsVariable * variables)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (fp != NULL);

  if (GFS_IS_MIXED (cell)) {
    GfsStateVector * s = GFS_STATE (cell);
    guint i;

    for (i = 0; i < FTT_NEIGHBORS; i++)
      fprintf (fp, " %g", s->solid->s[i]);
    fprintf (fp, " %g", s->solid->a);
    for (i = 0; i < FTT_DIMENSION; i++)
      fprintf (fp, " %g", (&s->solid->cm.x)[i]);
  }
  else
    fputs (" -1", fp);
  
  while (variables) {
    if (variables->name)
      fprintf (fp, " %g", GFS_VARIABLE (cell, variables->i));
    variables = variables->next;
  }
}

/**
 * gfs_cell_read:
 * @cell: a #FttCell.
 * @fp: a #GtsFile.
 * @domain: the #GfsDomain containing @cell.
 *
 * Reads from @fp the fluid data associated with @cell and described
 * by @domain->variables_io. This function is generally used in
 * association with ftt_cell_read().  
 */
void gfs_cell_read (FttCell * cell, GtsFile * fp, GfsDomain * domain)
{
  gdouble s0;
  GfsStateVector * s;
  GfsVariable * v;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (fp != NULL);
  g_return_if_fail (domain != NULL);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (solid->s[0])");
    return;
  }
  s0 = atof (fp->token->str);
  gts_file_next_token (fp);

  gfs_cell_init (cell, domain);
  s = cell->data;
  if (s0 >= 0.) {
    guint i;

    s->solid = g_malloc0 (sizeof (GfsSolidVector));
    s->solid->s[0] = s0;

    for (i = 1; i < FTT_NEIGHBORS; i++) {
      if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
	gts_file_error (fp, "expecting a number (solid->s[%d])", i);
	return;
      }
      s->solid->s[i] = atof (fp->token->str);
      gts_file_next_token (fp);
    }
    if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
      gts_file_error (fp, "expecting a number (solid->a)");
      return;
    }
    s->solid->a = atof (fp->token->str);
    gts_file_next_token (fp);
    for (i = 0; i < FTT_DIMENSION; i++) {
      if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
	gts_file_error (fp, "expecting a number (solid->cm[%d])", i);
	return;
      }
      (&s->solid->cm.x)[i] = atof (fp->token->str);
      gts_file_next_token (fp);
    }
  }

  v = domain->variables_io;
  while (v) {
    if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
      gts_file_error (fp, "expecting a number (%s)", v->name);
      return;
    }
    GFS_VARIABLE (cell, v->i) = atof (fp->token->str);
    gts_file_next_token (fp);
    v = v->next;
  }
}

/**
 * gfs_cell_write_binary:
 * @cell: a #FttCell.
 * @fp: a file pointer.
 * @variables: the #GfsVariable to be written.
 *
 * Writes in @fp the fluid data associated with @cell and described by
 * @variables. This function is generally used in association with
 * ftt_cell_write_binary().
 */
void gfs_cell_write_binary (const FttCell * cell, FILE * fp,
			    GfsVariable * variables)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (fp != NULL);

  if (GFS_IS_MIXED (cell)) {
    GfsStateVector * s = GFS_STATE (cell);

    fwrite (s->solid->s, sizeof (gdouble), FTT_NEIGHBORS, fp);
    fwrite (&s->solid->a, sizeof (gdouble), 1, fp);
    fwrite (&s->solid->cm.x, sizeof (gdouble), FTT_DIMENSION, fp);
  }
  else {
    gdouble a = -1.;
    fwrite (&a, sizeof (gdouble), 1, fp);
  }
  
  while (variables) {
    if (variables->name) {
      gdouble a = GFS_VARIABLE (cell, variables->i);
      fwrite (&a, sizeof (gdouble), 1, fp);
    }
    variables = variables->next;
  }
}

/**
 * gfs_cell_read_binary:
 * @cell: a #FttCell.
 * @fp: a #GtsFile.
 * @domain: the #GfsDomain containing @cell.
 *
 * Reads from @fp the fluid data associated with @cell and described
 * by @domain->variables_io. This function is generally used in
 * association with ftt_cell_read_binary().
 */
void gfs_cell_read_binary (FttCell * cell, GtsFile * fp, GfsDomain * domain)
{
  gdouble s0;
  GfsStateVector * s;
  GfsVariable * v;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (fp != NULL);
  g_return_if_fail (domain != NULL);

  if (gts_file_read (fp, &s0, sizeof (gdouble), 1) != 1) {
    gts_file_error (fp, "expecting a number (solid->s[0])");
    return;
  }

  gfs_cell_init (cell, domain);
  s = cell->data;
  if (s0 >= 0.) {
    s->solid = g_malloc0 (sizeof (GfsSolidVector));
    s->solid->s[0] = s0;
    
    if (gts_file_read (fp, &s->solid->s[1], sizeof (gdouble), FTT_NEIGHBORS - 1) != FTT_NEIGHBORS - 1) {
      gts_file_error (fp, "expecting numbers (solid->s[1..%d])", FTT_NEIGHBORS - 1);
      return;
    }
    if (gts_file_read (fp, &s->solid->a, sizeof (gdouble), 1) != 1) {
      gts_file_error (fp, "expecting a number (solid->a)");
      return;
    }
    if (gts_file_read (fp, &s->solid->cm.x, sizeof (gdouble), FTT_DIMENSION) != FTT_DIMENSION) {
      gts_file_error (fp, "expecting numbers (solid->cm[0..%d])", FTT_DIMENSION - 1);
      return;
    }
  }

  v = domain->variables_io;
  while (v) {
    gdouble a;

    if (gts_file_read (fp, &a, sizeof (gdouble), 1) != 1) {
      gts_file_error (fp, "expecting a number (%s)", v->name);
      return;
    }
    GFS_VARIABLE (cell, v->i) = a;
    v = v->next;
  }
}

/**
 * gfs_domain_replace_variable:
 * @domain: a #GfsDomain.
 * @v: the #GfsVariable to replace.
 * @with: the new #GfsVariable.
 *
 * Replaces existing variable @v with new variable @with.
 */
void gfs_domain_replace_variable (GfsDomain * domain,
				  GfsVariable * v,
				  GfsVariable * with)
{
  GfsVariable * v1, * prev = NULL;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (gts_container_size (GTS_CONTAINER (domain)) == 0);
  g_return_if_fail (v != NULL);
  g_return_if_fail (with != NULL);

  v1 = domain->variables;
  while (v1 && v1 != v) {
    prev = v1;
    v1 = v1->next;
  }
  g_return_if_fail (v1 == v);
  with->i = v->i;
  v->i = -1;
  with->p = GTS_OBJECT (domain);
  v->p = NULL;
  with->next = v->next;
  v->next = NULL;
  if (prev)
    prev->next = with;
  else
    domain->variables = with;
}

/**
 * gfs_domain_add_new_variable:
 * @domain: a #GfsDomain.
 * @v: the #GfsVariable to add.
 *
 * Adds a new variable @v to @domain.
 */
void gfs_domain_add_new_variable (GfsDomain * domain,
				  GfsVariable * v)
{
  GfsVariable * v1, * last;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (gts_container_size (GTS_CONTAINER (domain)) == 0);
  g_return_if_fail (v != NULL);
  g_return_if_fail (v->name == NULL || 
		    gfs_variable_from_name (domain->variables, v->name) == NULL);

  v1 = last = domain->variables;
  while (v1) {
    last = v1;
    v1 = v1->next;
  }
  g_assert (last);

  last->next = v;
  v->i = last->i + 1;
  v->p = GTS_OBJECT (domain);
  domain->variables_size += sizeof (gdouble);
}

/**
 * gfs_domain_add_variable:
 * @domain: a #GfsDomain.
 * @name: the name of the variable to add or %NULL.
 *
 * Adds a new variable @name to @domain.
 *
 * Returns: the new variable or %NULL if a variable with the same name
 * already exists.  
 */
GfsVariable * gfs_domain_add_variable (GfsDomain * domain, 
				       const gchar * name)
{
  GfsVariable * v;

  g_return_val_if_fail (domain != NULL, NULL);
  g_return_val_if_fail (gts_container_size (GTS_CONTAINER (domain)) == 0, NULL);

  if (name && gfs_variable_from_name (domain->variables, name))
    return NULL;

  v = GFS_VARIABLE1 (gts_object_new (GTS_OBJECT_CLASS (gfs_variable_class ())));
  v->name = g_strdup (name);
  gfs_domain_add_new_variable (domain, v);

  return v;
}

static void add_pressure_force (FttCell * cell, gdouble * f)
{
  GfsSolidVector * s = GFS_STATE (cell)->solid;
  gdouble p = gfs_cell_dirichlet_value (cell, gfs_p, -1);
  FttComponent c;
  gdouble size = ftt_cell_size (cell);

#if (!FTT_2D)
  size *= size;
#endif /* 3D */

  for (c = 0; c < FTT_DIMENSION; c++)
    f[c] += p*(s->s[2*c + 1] - s->s[2*c])*size;
}

static GfsSourceDiffusion * source_diffusion (GfsVariable * v)
{
  if (v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
    
    while (i) {
      GtsObject * o = i->data;
      
      if (GFS_IS_SOURCE_DIFFUSION (o))
	return GFS_SOURCE_DIFFUSION (o);
      i = i->next;
    }
  }
  return NULL;
}

static void add_viscous_force (FttCell * cell, gpointer * data)
{
  FttVector * f = data[0];
  GfsVariable * v = data[1];
  GfsSourceDiffusion * d = data[2];
  gdouble D;
  GfsSolidVector * s = GFS_STATE (cell)->solid;
  FttVector n, g;

  g_assert (((cell)->flags & GFS_FLAG_DIRICHLET) != 0);
  gfs_cell_dirichlet_gradient (cell, v->i, -1, s->fv, &g);

  D = gfs_source_diffusion_cell (d, cell);
  n.x = s->s[1] - s->s[0];
  n.y = s->s[3] - s->s[2];
#if FTT_2D
  switch (GFS_VELOCITY_COMPONENT (v->i)) {
  case FTT_X:
    f->x -= D*(2.*g.x*n.x + g.y*n.y);
    f->y -= D*g.y*n.x;
    break;
  case FTT_Y:
    f->x -= D*g.x*n.y;
    f->y -= D*(2.*g.y*n.y + g.x*n.x);
    break;
  default:
    g_assert_not_reached ();
  }
#else /* 3D */
  n.z = s->s[5] - s->s[4];
  D *= ftt_cell_size (cell);
  switch (GFS_VELOCITY_COMPONENT (v->i)) {
  case FTT_X:
    f->x -= D*(2.*g.x*n.x + g.y*n.y + g.z*n.z);
    f->y -= D*g.y*n.x;
    f->z -= D*g.z*n.x;
    break;
  case FTT_Y:
    f->y -= D*(2.*g.y*n.y + g.x*n.x + g.z*n.z);
    f->x -= D*g.x*n.y;
    f->z -= D*g.z*n.y;
    break;
  case FTT_Z:
    f->z -= D*(2.*g.z*n.z + g.x*n.x + g.y*n.y);
    f->x -= D*g.x*n.z;
    f->y -= D*g.y*n.z;
    break;
  default:
    g_assert_not_reached ();
  }
#endif /* 3D */
}

/**
 * gfs_domain_solid_force:
 * @domain: a #GfsDomain.
 * @pf: a #FttVector.
 * @vf: a #FttVector.
 *
 * Fills @pf and @vf with the components of the net pressure and
 * viscous forces applied by the fluid on the solid surface embbeded
 * in @domain.
 */
void gfs_domain_solid_force (GfsDomain * domain, 
			     FttVector * pf,
			     FttVector * vf)
{
  FttComponent c;
  GfsVariable * v;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (pf != NULL);
  g_return_if_fail (vf != NULL);

  pf->x = pf->y = pf->z = 0.;
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseFunc) add_pressure_force, pf);

  vf->x = vf->y = vf->z = 0.;
  v = gfs_variable_from_name (domain->variables, "U");
  for (c = 0; c < FTT_DIMENSION; c++, v = v->next) {
    GfsSourceDiffusion * D = source_diffusion (v);

    if (D) {
      gpointer data[3];

      gfs_domain_surface_bc (domain, v);
      data[0] = vf;
      data[1] = v;
      data[2] = D;
      gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
				 (FttCellTraverseFunc) add_viscous_force, 
				 data);
    }
  }
}

static void tag_cell_fraction (FttCell * cell, GfsVariable * c, guint tag, guint * size)
{
  FttDirection d;
  FttCellNeighbors n;

  g_assert (FTT_CELL_IS_LEAF (cell));
  GFS_STATE (cell)->div = tag;
  (*size)++;
  ftt_cell_neighbors (cell, &n);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (n.c[d] && GFS_STATE (n.c[d])->div == 0. &&
	GFS_VARIABLE (n.c[d], c->i) > 1e-4) {
      if (FTT_CELL_IS_LEAF (n.c[d]))
	tag_cell_fraction (n.c[d], c, tag, size);
      else {
	FttCellChildren child;
	FttDirection od = FTT_OPPOSITE_DIRECTION (d);
	guint i;

#if FTT_2D3
	g_assert_not_implemented ();
#endif	
	ftt_cell_children_direction (n.c[d], od, &child);
	for (i = 0; i < FTT_CELLS/2; i++)
	  if (child.c[i] && GFS_STATE (child.c[i])->div == 0. &&
	      GFS_VARIABLE (child.c[i], c->i) > 1e-4)
	    tag_cell_fraction (child.c[i], c, tag, size);
      }
    }
}

static void tag_new_fraction_region (FttCell * cell, gpointer * data)
{
  if (GFS_STATE (cell)->div == 0.) {
    GfsVariable * c = data[0];
    GArray * sizes = data[1];
    guint size = 0;
    
    tag_cell_fraction (cell, c, sizes->len + 1, &size);
    g_array_append_val (sizes, size);
  }
}

static void reset_small_fraction (FttCell * cell, gpointer * data)
{
  GfsVariable * c = data[0];
  GArray * sizes = data[1];
  guint * min = data[2], i = GFS_STATE (cell)->div - 1.;
  
  g_assert (GFS_STATE (cell)->div > 0.);
  if (g_array_index (sizes, guint, i) < *min)
    GFS_VARIABLE (cell, c->i) = 0.;
}

static int greater (const void * a, const void * b)
{
  return *((guint *)a) > *((guint *)b) ? -1 : 1;
}

/**
 * gfs_domain_remove_droplets:
 * @domain: a #GfsDomain.
 * @c: a #GfsVariable.
 * @min: the minimum size (in cells) of the droplets.
 *
 * Removes all the droplets of @domain smaller than @min cells
 * if @min is positive, or all the droplets but the -$min largest ones
 * if @min is negative.
 */
void gfs_domain_remove_droplets (GfsDomain * domain,
				 GfsVariable * c,
				 gint min)
{
  GArray * sizes;
  gpointer data[3];
  guint minsize;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (c != NULL);

  sizes = g_array_new (FALSE, FALSE, sizeof (guint));
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, gfs_div);
  data[0] = c;
  data[1] = sizes;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) tag_new_fraction_region, data);
  g_assert (sizes->len > 0);
  if (min >= 0)
    minsize = min;
  else if (-min >= sizes->len)
    minsize = 0;
  else {
    guint * tmp = g_malloc (sizes->len*sizeof (guint));
    memcpy (tmp, sizes->data, sizes->len*sizeof (guint));
    qsort (tmp, sizes->len, sizeof (guint), greater);
    minsize = tmp[-1 - min];
    g_free (tmp);
  }
  data[2] = &minsize;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) reset_small_fraction, data);
  g_array_free (sizes, TRUE);
}

static void tag_cell (FttCell * cell, guint tag, guint * size)
{
  FttDirection d;
  FttCellNeighbors n;
  GfsSolidVector * solid = GFS_STATE (cell)->solid;

  g_assert (FTT_CELL_IS_LEAF (cell));
  GFS_STATE (cell)->div = tag;
  (*size)++;
  ftt_cell_neighbors (cell, &n);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (n.c[d] && GFS_STATE (n.c[d])->div == 0. &&
	!GFS_CELL_IS_BOUNDARY (n.c[d]) &&
	(!solid || solid->s[d] > 0.)) {
      if (FTT_CELL_IS_LEAF (n.c[d]))
	tag_cell (n.c[d], tag, size);
      else {
	FttCellChildren child;
	FttDirection od = FTT_OPPOSITE_DIRECTION (d);
	guint i, j;
	
	j = ftt_cell_children_direction (n.c[d], od, &child);
	for (i = 0; i < j; i++)
	  if (child.c[i] && GFS_STATE (child.c[i])->div == 0. &&
	      (!GFS_IS_MIXED (child.c[i]) || GFS_STATE (child.c[i])->solid->s[od] > 0.))
	    tag_cell (child.c[i], tag, size);
      }
    }
}

static void tag_new_region (FttCell * cell, GArray * sizes)
{
  if (GFS_STATE (cell)->div == 0.) {
    guint size = 0;

    tag_cell (cell, sizes->len + 1, &size);
    g_array_append_val (sizes, size);
  }
}

static gboolean remove_small (FttCell * cell, gpointer * data)
{
  if (FTT_CELL_IS_LEAF (cell)) {
    GArray * sizes = data[0];
    guint * min = data[1], i = GFS_STATE (cell)->div - 1.;

    g_assert (GFS_STATE (cell)->div > 0.);
    if (g_array_index (sizes, guint, i) < *min) {
      if (FTT_CELL_IS_ROOT (cell))
	g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, "root cell belongs to a pond");
      else
	ftt_cell_destroy (cell, data[2], data[3]);
      return TRUE;
    }
    return FALSE;
  }
  else {
    FttCellChildren child;
    guint i;
    gboolean changed = FALSE;

    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i] && remove_small (child.c[i], data))
	changed = TRUE;
    if (FTT_CELL_IS_LEAF (cell)) {
      /* all the children have been destroyed i.e. the cell belongs to a small pond */
      if (FTT_CELL_IS_ROOT (cell))
	g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, "root cell belongs to a pond");
      else
	ftt_cell_destroy (cell, data[2], data[3]);
    }
    else if (changed)
      gfs_cell_init_solid_fractions_from_children (cell);
    return changed;
  }
}

static void remove_small_box (GfsBox * box, gpointer * data)
{
  gboolean * changed = data[4];

  if (remove_small (box->root, data))
    *changed = TRUE;
}

/**
 * gfs_domain_remove_ponds:
 * @domain: a #GfsDomain.
 * @min: the minimum size (in cells) of the ponds.
 * @cleanup: a #FttCellCleanupFunc or %NULL.
 * @data: user data to pass to @cleanup.
 *
 * Removes all the fluid "ponds" of @domain smaller than @min cells
 * if @min is positive, or all the ponds but the -@min largest ones
 * if @min is negative.
 *
 * If the domain is modified its boundaries are re"matched" using
 * gfs_domain_match().
 */
void gfs_domain_remove_ponds (GfsDomain * domain, 
			      gint min,
			      FttCellCleanupFunc cleanup,
			      gpointer data)
{
  GArray * sizes;
  gpointer dat[5];
  guint minsize;
  gboolean changed = FALSE;

  g_return_if_fail (domain != NULL);

  sizes = g_array_new (FALSE, FALSE, sizeof (guint));
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, gfs_div);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) tag_new_region, sizes);
  g_assert (sizes->len > 0);
  if (min >= 0)
    minsize = min;
  else if (-min >= sizes->len)
    minsize = 0;
  else {
    guint * tmp = g_malloc (sizes->len*sizeof (guint));
    memcpy (tmp, sizes->data, sizes->len*sizeof (guint));
    qsort (tmp, sizes->len, sizeof (guint), greater);
    minsize = tmp[-1 - min];
    g_free (tmp);
  }
  dat[0] = sizes;
  dat[1] = &minsize;
  dat[2] = cleanup;
  dat[3] = data;
  dat[4] = &changed;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) remove_small_box, dat);
  g_array_free (sizes, TRUE);
  if (changed)
    gfs_domain_match (domain);
}

static gboolean tag_speck (FttCell * cell)
{
  if (GFS_STATE (cell)->div == 0.) {
    FttDirection d;
    FttCellNeighbors n;
    GfsSolidVector * solid = GFS_STATE (cell)->solid;
    
    g_assert (FTT_CELL_IS_LEAF (cell));
    ftt_cell_neighbors (cell, &n);
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (!n.c[d])
	return FALSE;
    GFS_STATE (cell)->div = 1.;
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (GFS_STATE (n.c[d])->div == 0. && 
	  !GFS_CELL_IS_BOUNDARY (n.c[d]) &&
	  solid->s[d] > 0. && solid->s[d] < 1.) {
	g_assert (GFS_IS_MIXED (n.c[d]));
	if (FTT_CELL_IS_LEAF (n.c[d])) {
	  if (!tag_speck (n.c[d])) {
	    GFS_STATE (cell)->div = 0.;
	    return FALSE;
	  }
	}
	else {
	  FttCellChildren child;
	  FttDirection od = FTT_OPPOSITE_DIRECTION (d);
	  guint i;
	  
#if FTT_2D3
	  g_assert_not_implemented ();
#endif	
	  ftt_cell_children_direction (n.c[d], od, &child);
	  for (i = 0; i < FTT_CELLS/2; i++)
	    if (!child.c[i] || (GFS_STATE (child.c[i])->div == 0 && 
				GFS_IS_MIXED (child.c[i]) &&
				!tag_speck (child.c[i]))) {
	      GFS_STATE (cell)->div = 0.;
	      return FALSE;
	    }
	}
      }
  }
  return TRUE;
}

static void fill_speck (FttCell * cell, gboolean * changed)
{
  if (GFS_STATE (cell)->div == 1.) {
    g_free (GFS_STATE (cell)->solid);
    GFS_STATE (cell)->solid = NULL;
    *changed = TRUE;
  }
}

/**
 * gfs_domain_remove_specks:
 * @domain: a #GfsDomain.
 *
 * Removes all the solid "specks" of @domain. Solid specks are islands
 * which do not contain any empty cell.
 *
 * Note that the domain's boundaries are not "matched" automatically.
 */
void gfs_domain_remove_specks (GfsDomain * domain)
{
  gboolean changed = FALSE;

  g_return_if_fail (domain != NULL);

  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL,
			     (FttCellTraverseFunc) gfs_cell_reset, gfs_div);
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseFunc) tag_speck, NULL);
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseFunc) fill_speck, &changed);
  if (changed)
    gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_init_solid_fractions_from_children, 
			      NULL);
    
}

/**
 * gfs_domain_timer_start:
 * @domain: a #GfsDomain.
 * @name: the name of the timer.
 *
 * Starts timer @name of @domain. If @name does not exist it is
 * created first.
 */
void gfs_domain_timer_start (GfsDomain * domain, const gchar * name)
{
  GfsTimer * t;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (name != NULL);

  t = g_hash_table_lookup (domain->timers, name);
  if (t == NULL) {
    t = g_malloc (sizeof (GfsTimer));
    gts_range_init (&t->r);
    g_hash_table_insert (domain->timers, g_strdup (name), t);
  }
  else
    g_return_if_fail (t->start < 0.);
  t->start = g_timer_elapsed (domain->timer, NULL);
}

/**
 * gfs_domain_timer_stop:
 * @domain: a #GfsDomain.
 * @name: the name of the timer.
 *
 * Stops timer @name of @domain. This function fails if @name is not a
 * timer of @domain.
 */
void gfs_domain_timer_stop (GfsDomain * domain, const gchar * name)
{
  GfsTimer * t;
  gdouble end;

  g_return_if_fail (domain != NULL);
  end = g_timer_elapsed (domain->timer, NULL);
  g_return_if_fail (name != NULL);

  t = g_hash_table_lookup (domain->timers, name);
  g_return_if_fail (t != NULL);
  g_return_if_fail (t->start >= 0.);

  gts_range_add_value (&t->r, end - t->start);
  gts_range_update (&t->r);
  t->start = -1.;
}

static void cell_combine_traverse (FttCell * cell,
				   FttCell * parent,
				   FttCellCombineTraverseFunc inside,
				   gpointer idata,
				   FttCellTraverseFunc outside,
				   gpointer odata)
{
  FttCell * locate;
  FttVector p;

  ftt_cell_pos (cell, &p);
  locate = ftt_cell_locate (parent, p, ftt_cell_level (cell));
  if (locate == NULL) {
    if (outside)
      ftt_cell_traverse (cell, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, outside, odata);
  }
  else {
    if (FTT_CELL_IS_LEAF (cell))
      (* inside) (cell, locate, idata);
    else {
      FttCellChildren child;
      guint i;

      ftt_cell_children (cell, &child);
      for (i = 0; i < FTT_CELLS; i++)
	if (child.c[i])
	  cell_combine_traverse (child.c[i], locate, inside, idata, outside, odata);
    }
  }  
}

static void box_combine_traverse (GfsBox * box, gpointer * data)
{
  FttVector p;
  FttCell * locate;

  ftt_cell_pos (box->root, &p);
  locate = gfs_domain_locate (data[0], p, ftt_cell_level (box->root));
  if (locate == NULL) {
    if (data[3])
      ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, data[3], data[4]);
  }
  else
    cell_combine_traverse (box->root, locate, data[1], data[2], data[3], data[4]);
}

/**
 * gfs_domain_combine_traverse:
 * @domain1: a #GfsDomain.
 * @domain2: another #GfsDomain.
 * @inside: function to call for each pair of cells.
 * @idata: user data to pass to @inside.
 * @outside: function to call for cells falling outside of @domain2 or
 * %NULL.
 * @odata: user data to pass to @outside.
 *
 * Calls @inside for each leaf cell of @domain1 contained in
 * @domain2. The second cell argument to @inside is set to the cell of
 * @domain2 containing the first cell argument.
 *
 * If @outside is not %NULL it is called for each leaf cell of
 * @domain1 which is outside of @domain2.
 */
void gfs_domain_combine_traverse (GfsDomain * domain1,
				  GfsDomain * domain2,
				  FttCellCombineTraverseFunc inside,
				  gpointer idata,
				  FttCellTraverseFunc outside,
				  gpointer odata)				  
{
  gpointer data[5];

  g_return_if_fail (domain1 != NULL);
  g_return_if_fail (domain2 != NULL);
  g_return_if_fail (inside != NULL);

  data[0] = domain2;
  data[1] = inside;
  data[2] = idata;
  data[3] = outside;
  data[4] = odata;

  gts_container_foreach (GTS_CONTAINER (domain1), (GtsFunc) box_combine_traverse, data);
}
