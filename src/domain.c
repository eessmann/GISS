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
    GSList * i = domain->variables_io;

    if (i != NULL) {
      fprintf (fp, "variables = %s", GFS_VARIABLE1 (i->data)->name);
      i = i->next;
      while (i) {
	fprintf (fp, ",%s", GFS_VARIABLE1 (i->data)->name);
	i = i->next;
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

#if FTT_2D
  if (var[3].set) {
    gts_file_variable_error (fp, var, "z", "unknown identifier `z'");
    return;
  }
  if (var[6].set) {
    gts_file_variable_error (fp, var, "lz", "unknown identifier `lz'");
    return;
  }
#endif

  if (var[4].set && domain->lambda.x <= 0.) {
    gts_file_variable_error (fp, var, "lx", "lx must be strictly positive");
    return;
  }
  if (var[5].set && domain->lambda.y <= 0.) {
    gts_file_variable_error (fp, var, "ly", "ly must be strictly positive");
    return;
  }
  if (var[6].set && domain->lambda.z <= 0.) {
    gts_file_variable_error (fp, var, "lz", "lz must be strictly positive");
    return;
  }

  if (variables != NULL) {
    gchar * variables1, * s;

    variables1 = g_strdup (variables);
    s = strtok (variables1, ",");
    while (s) {
      gfs_domain_add_variable (domain, s, NULL);
      s = strtok (NULL, ",");
    }
    g_free (variables1);
    domain->variables_io = gfs_variables_from_list (domain->variables, variables, &s);
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
  else {
    FttDirection d;
    GfsBox * matching;

    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (GFS_IS_BOUNDARY_PERIODIC (box->neighbor[d]) &&
	  (matching = GFS_BOUNDARY_PERIODIC (box->neighbor[d])->matching)->pid != domain->pid) {
	gts_object_destroy (GTS_OBJECT (box->neighbor[d]));
	gfs_boundary_mpi_new (gfs_boundary_mpi_class (), box, d, matching->pid, matching->id);
      }
  }
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
    if (GFS_IS_BOX (box->neighbor[d]) && GFS_BOX (box->neighbor[d])->pid == domain->pid)
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
  GSList * i;

  gfs_clock_destroy (domain->timer);

  i = domain->variables;
  while (i) {
    GSList * next = i->next;
    gts_object_destroy (i->data);
    i = next;
  }
  g_assert (domain->variables == NULL);

  g_slist_foreach (domain->derived_variables, (GFunc) gts_object_destroy, NULL);
  g_slist_free (domain->derived_variables);
  domain->derived_variables = NULL;

  g_array_free (domain->allocated, TRUE);

  g_hash_table_foreach (domain->timers, (GHFunc) free_pair, NULL);
  g_hash_table_destroy (domain->timers);

  g_slist_free (domain->variables_io);

  (* GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class->destroy) (o);
}

static void add_item (gpointer item, GPtrArray * a)
{
  g_ptr_array_add (a, item);
}

static int compare_boxes (const void * b1, const void * b2)
{
  return (*(GfsBox **)b1)->id < (*(GfsBox **)b2)->id ? -1 : 1;
}

static void domain_foreach (GtsContainer * c, 
			    GtsFunc func, 
			    gpointer data)
{
  GPtrArray * a = g_ptr_array_new ();
  (* GTS_CONTAINER_CLASS (GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class)->foreach)
    (c, (GtsFunc) add_item, a);
  qsort (a->pdata, a->len, sizeof (gpointer), compare_boxes);
  guint i;
  for (i = 0; i < a->len; i++)
    (* func) (a->pdata[i], data);
  g_ptr_array_free (a, TRUE);
}

static void domain_class_init (GfsDomainClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = domain_read;
  GTS_OBJECT_CLASS (klass)->write = domain_write;
  GTS_OBJECT_CLASS (klass)->destroy = domain_destroy;

  GTS_CONTAINER_CLASS (klass)->foreach = domain_foreach;

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

  domain->timer = gfs_clock_new ();
  domain->timers = g_hash_table_new (g_str_hash, g_str_equal);

  gts_range_init (&domain->size);

  domain->profile_bc = FALSE;

  gts_range_init (&domain->mpi_messages);
  gts_range_init (&domain->mpi_wait);

  domain->rootlevel = 0;
  domain->refpos.x = domain->refpos.y = domain->refpos.z = 0.;
  domain->lambda.x = domain->lambda.y = domain->lambda.z = 1.;

  domain->allocated = g_array_new (FALSE, TRUE, sizeof (gboolean));
  domain->variables = NULL;

  domain->variables_io = NULL;
  domain->max_depth_write = -1;

  domain->cell_init = (FttCellInitFunc) gfs_cell_fine_init;
  domain->cell_init_data = domain;
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
	b->type = GFS_BOUNDARY_CENTER_VARIABLE;
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
  gfs_all_reduce (domain, changed, MPI_INT, MPI_MAX);
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

static gboolean is_velocity (GfsVariable * v, GfsDomain * domain)
{
  FttComponent c;
  GfsVariable ** u = gfs_domain_velocity (domain);

  for (c = 0; c < FTT_DIMENSION; c++)
    if (v == u[c])
      return TRUE;
  return FALSE;
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
  else if (is_velocity (v, domain))
    gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL,
			       (FttCellTraverseFunc) dirichlet_bc, NULL);
  else
    gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL,
			       (FttCellTraverseFunc) neumann_bc, NULL);
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

static void box_traverse_condition (GfsBox * box, gpointer * datum)
{
  FttTraverseType * order = datum[0];
  FttTraverseFlags * flags = datum[1];
  gint * max_depth = datum[2];
  FttCellTraverseFunc func = (FttCellTraverseFunc) datum[3];
  gpointer data = datum[4];
  gboolean (* condition) (FttCell *, gpointer) = datum[5];
  gpointer cdata = datum[6];

  ftt_cell_traverse_condition (box->root, *order, *flags, *max_depth, func, data,
			       condition, cdata);
}

/**
 * gfs_domain_cell_traverse_condition:
 * @domain: a #GfsDomain.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 * @condition: the condition.
 * @cdata: user data to pass to @condition.
 *
 * Traverses the cell trees of @domain. Calls the given function for
 * each cell visited.
 *
 * Traversal of any branch of the tree is stopped whenever @condition
 * is not verified.
 */
void gfs_domain_cell_traverse_condition (GfsDomain * domain,
					 FttTraverseType order,
					 FttTraverseFlags flags,
					 gint max_depth,
					 FttCellTraverseFunc func,
					 gpointer data,
					 gboolean (* condition) (FttCell *, gpointer),
					 gpointer cdata)
{
  gpointer datum[7];

  datum[0] = &order;
  datum[1] = &flags;
  datum[2] = &max_depth;
  datum[3] = func;
  datum[4] = data;
  datum[5] = condition;
  datum[6] = cdata;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (func != NULL);
  g_return_if_fail (condition != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_traverse_condition, datum);
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

typedef struct {
  FttCellTraverseCutFunc func;
  gpointer data;
  FttTraverseType order;
  FttTraverseFlags flags;
  GfsGenericSurface * s;
} TraverseCut;

static void traverse_cut (GfsBox * box, TraverseCut * p)
{
  gfs_cell_traverse_cut (box->root, p->s, p->order, p->flags, p->func, p->data);
}

/**
 * gfs_domain_traverse_cut:
 * @domain: a #GfsDomain.
 * @s: a #GfsGenericSurface.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Calls @func for each cell of @domain cut by @s.
 */
void gfs_domain_traverse_cut (GfsDomain * domain,
			      GfsGenericSurface * s,
			      FttTraverseType order,
			      FttTraverseFlags flags,
			      FttCellTraverseCutFunc func,
			      gpointer data)
{
  TraverseCut p;

  p.func = func;
  p.data= data;
  p.order = order;
  p.flags = flags;
  p.s = s;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (func != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) traverse_cut, &p);
}

static void traverse_cut_2D (GfsBox * box, TraverseCut * p)
{
  gfs_cell_traverse_cut_2D (box->root, p->s, p->order, p->flags, p->func, p->data);
}

/**
 * gfs_domain_traverse_cut_2D:
 * @domain: a #GfsDomain.
 * @s: a #GfsGenericSurface.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Calls @func for each cell of @domain cut by @s.
 *
 * The cells are flattened in the z-direction.
 */
void gfs_domain_traverse_cut_2D (GfsDomain * domain,
				 GfsGenericSurface * s,
				 FttTraverseType order,
				 FttTraverseFlags flags,
				 FttCellTraverseCutFunc func,
				 gpointer data)
{
  TraverseCut p;

  p.func = func;
  p.data = data;
  p.order = order;
  p.flags = flags;
  p.s = s;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (func != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) traverse_cut_2D, &p);
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
  gfs_all_reduce (domain, depth, MPI_UNSIGNED, MPI_MAX);
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
  gdouble v = GFS_VARIABLE (cell, GFS_VARIABLE1 (data[1])->i);

  if (v < G_MAXDOUBLE)
    gts_range_add_value (s, v);
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

#define BPID(b) ((b)->pid > 0 ? (b)->pid : 0)

static void box_count (GfsBox * b, GArray * a)
{
  guint count = 0, pid = BPID(b);
  ftt_cell_traverse (b->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) cell_count, &count);
  if (pid >= a->len)
    g_array_set_size (a, pid + 1);
  g_array_index (a, guint, pid) += count;
}

static void boundary_size (GfsBox * box, GArray * a)
{
  FttDirection d;
  guint count = 0;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (
#ifdef HAVE_MPI
	GFS_IS_BOUNDARY_MPI (box->neighbor[d]) ||
#endif
	(GFS_IS_BOX (box->neighbor[d]) && GFS_BOX (box->neighbor[d])->pid != box->pid)
       )
      ftt_cell_traverse_boundary (box->root, d, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				  (FttCellTraverseFunc) cell_count, &count);
  g_array_index (a, guint, BPID (box)) += count;
}

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
  g_return_if_fail (domain != NULL);
  g_return_if_fail (size != NULL);
  g_return_if_fail (boundary != NULL);
  g_return_if_fail (mpiwait != NULL);

  gts_range_init (size);
  gts_range_init (boundary);
  gts_range_init (mpiwait);

  if (domain->timestep.n > 0)
    gts_range_add_value (mpiwait, domain->mpi_wait.sum/domain->timestep.n);

  GArray * a = g_array_new (FALSE, TRUE, sizeof (guint));
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_count, a);
  guint i;
  for (i = 0; i < a->len; i++) {
    guint v = g_array_index (a, guint, i);
    if (v > 0) {
      gts_range_add_value (size, v);
      g_array_index (a, guint, i) = 0;
    }
  }
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) boundary_size, a);
  for (i = 0; i < a->len; i++) {
    guint v = g_array_index (a, guint, i);
    if (v > 0)
      gts_range_add_value (boundary, v);
  }
#ifdef HAVE_MPI
  domain_range_reduce (domain, size);
  domain_range_reduce (domain, boundary);
  domain_range_reduce (domain, mpiwait);
#endif /* HAVE_MPI */
  g_array_free (a, TRUE);
  gts_range_update (size);
  gts_range_update (boundary);
  gts_range_update (mpiwait);
}

static void add_norm (const FttCell * cell, gpointer * data)
{
  GfsNorm * n = data[0];
  GfsVariable * v = data[1];

  gfs_norm_add (n, GFS_VARIABLE (cell, v->i), gfs_cell_volume (cell));
}

static void add_norm_weighted (FttCell * cell, gpointer * data)
{
  GfsNorm * n = data[0];
  GfsVariable * v = data[1];
  GfsFunction * w = data[2];

  gfs_norm_add (n, GFS_VARIABLE (cell, v->i), gfs_cell_volume (cell)*gfs_function_value (w, cell));
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
 * @w: a #GfsFunction or %NULL.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 *
 * Traverses the domain defined by @domain using gfs_domain_cell_traverse()
 * and gathers norm statistics about variable @v.
 *
 * The norm is weighted by the volume of each cell times the value of
 * function @w (if @w is not %NULL).
 *
 * Returns: a #GfsNorm containing the norm statistics about @v.
 */
GfsNorm gfs_domain_norm_variable (GfsDomain * domain,
				  GfsVariable * v,
				  GfsFunction * w,
				  FttTraverseFlags flags,
				  gint max_depth)
{
  GfsNorm n;
  gpointer data[3];

  g_return_val_if_fail (domain != NULL, n);
  g_return_val_if_fail (v != NULL, n);
  
  gfs_norm_init (&n);
  data[0] = &n;
  data[1] = v;
  data[2] = w;
  if (w)
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			      (FttCellTraverseFunc) add_norm_weighted, data);
  else
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			      (FttCellTraverseFunc) add_norm, data);
#ifdef HAVE_MPI
  domain_norm_reduce (domain, &n);
#endif /* HAVE_MPI */
  gfs_norm_update (&n);

  return n;
}

static void add_norm_residual (const FttCell * cell, gpointer * data)
{
  gdouble size = ftt_cell_size (cell);
  GfsVariable * res = data[0];
  GfsNorm * n = data[1];
  
  gfs_norm_add (n, GFS_VARIABLE (cell, res->i)/(size*size), 1.);
}

/**
 * gfs_domain_norm_residual:
 * @domain: the domain to obtain the norm from.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 * @dt: the time step.
 * @res: the residual.
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
				  gdouble dt,
				  GfsVariable * res)
{
  GfsNorm n;
  gpointer data[2];

  g_return_val_if_fail (domain != NULL, n);
  g_return_val_if_fail (res != NULL, n);
  
  gfs_norm_init (&n);
  data[0] = res;
  data[1] = &n;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			   (FttCellTraverseFunc) add_norm_residual, data);
#ifdef HAVE_MPI
  domain_norm_reduce (domain, &n);
#endif /* HAVE_MPI */
  gfs_norm_update (&n);

  dt *= dt;
  n.bias *= dt;
  n.first *= dt;
  n.second *= dt;
  n.infty *= dt;
  return n;
}

/**
 * gfs_domain_velocity:
 * @domain: a #GfsDomain.
 *
 * Returns: the components of the velocity vector for @domain.
 */
GfsVariable ** gfs_domain_velocity (GfsDomain * domain)
{
  FttComponent c;
  static gchar name[][2] = {"U","V","W"};

  g_return_val_if_fail (domain != NULL, NULL);
  
  for (c = 0; c < FTT_DIMENSION; c++) {
    GfsVariable * v = gfs_variable_from_name (domain->variables, name[c]);
    g_return_val_if_fail (v != NULL, NULL);
    domain->velocity[c] = v;
  }
  return domain->velocity;
}

static void add_norm_velocity (FttCell * cell, gpointer * data)
{
  GfsVariable ** u = data[0];
  GfsNorm * n = data[1];
  
  gfs_norm_add (n, gfs_vector_norm (cell, u), gfs_cell_volume (cell));
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
  gpointer data[2];

  g_return_val_if_fail (domain != NULL, n);
  
  gfs_norm_init (&n);
  data[0] = gfs_domain_velocity (domain);
  data[1] = &n;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			   (FttCellTraverseFunc) add_norm_velocity, data);
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
  if (fp->type == GTS_ERROR) {
    gts_object_destroy (GTS_OBJECT (domain));
    return NULL;
  }

  return domain;
}

typedef struct {
  GSList * boxlist;
  guint bid;
  gboolean one_box_per_pe;
  gint pid;
  GfsVariable * newboxp;
  GfsDomain * domain;
} SplitPar;

static void box_split (GfsBox * box, SplitPar * p)
{
  guint refid = FTT_DIMENSION == 2 ? 2 : 6;
  FttCellChildren child;
  FttDirection d;
  guint i;
  GfsDomain * domain = gfs_box_domain (box);

  p->boxlist = g_slist_prepend (p->boxlist, box);

  if (FTT_CELL_IS_LEAF (box->root))
    ftt_cell_refine_single (box->root, (FttCellInitFunc) gfs_cell_init, domain);

  ftt_cell_children (box->root, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
      GfsBox * newbox = GFS_BOX (gts_object_new (GTS_OBJECT (box)->klass));

      GTS_OBJECT (newbox)->reserved = domain;
      if (p->one_box_per_pe)
	newbox->pid = (p->pid)++;
      else
	newbox->pid = box->pid;
      if (box->id == 1 && i == refid)
	newbox->id = 1;
      else
	newbox->id = (p->bid)++;

      GFS_DOUBLE_TO_POINTER (GFS_VARIABLE (child.c[i], p->newboxp->i)) = newbox;

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
	  GfsBox * newbox = GFS_DOUBLE_TO_POINTER (GFS_VARIABLE (child.c[i], p->newboxp->i));
	  GfsBoundaryClass * klass = GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass);
	  GtsObject * newboundary = GTS_OBJECT (gfs_boundary_new (klass, newbox, d));

	  if (GFS_IS_BOUNDARY_PERIODIC (newboundary))
	    GFS_BOUNDARY_PERIODIC (newboundary)->matching = 
	      GFS_BOUNDARY_PERIODIC (boundary)->matching;
	  else {
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
	    (* GTS_OBJECT_CLASS (klass)->read) (&newboundary, gfp);
	    g_assert (gfp->type != GTS_ERROR);
	    gts_file_destroy (gfp);
	    fclose (fp);
	  }
	  g_assert (neighbor);
	  GFS_BOUNDARY (newboundary)->root = neighbor;
	}
      gts_object_destroy (GTS_OBJECT (boundary));
    }
}

static GtsGEdge * node_is_linked (GtsGNode * n1, GtsGNode * n2, FttDirection d)
{
  GSList * i = GTS_SLIST_CONTAINER (n1)->items;
  while (i) {
    if (GTS_GNODE_NEIGHBOR (n1, i->data) == n2 &&
	GFS_GEDGE (i->data)->d == d)
      return i->data;
    i = i->next;
  }
  return NULL;
}

static void box_link (GfsBox * box, SplitPar * p)
{
  FttCellChildren child;
  guint i;

  ftt_cell_children (box->root, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
       GfsBox * newbox = GFS_DOUBLE_TO_POINTER (GFS_VARIABLE (child.c[i], p->newboxp->i));
       FttDirection d;
       
       g_assert (newbox);
       gts_container_add (GTS_CONTAINER (p->domain), GTS_CONTAINEE (newbox));

       for (d = 0; d < FTT_NEIGHBORS; d++)
	 if (newbox->neighbor[d] != NULL && GFS_IS_BOUNDARY_PERIODIC (newbox->neighbor[d])) {
	   GfsBox * matching =  GFS_BOUNDARY_PERIODIC (newbox->neighbor[d])->matching;
	   static FttDirection match[FTT_CELLS][FTT_DIMENSION] = {
#if FTT_2D
	     {0,2}, {1,2}, {0,3}, {1,3}
#elif FTT_2D3
#else /* 3D */
	     {0,2,4}, {1,2,4}, {0,3,4}, {1,3,4},
	     {0,2,5}, {1,2,5}, {0,3,5}, {1,3,5}
#endif /* 3D */
	   };
	   FttCell * neighbor = ftt_cell_child_corner (matching->root, 
						       match[FTT_CELL_ID (child.c[i])]);
	   g_assert (neighbor);
	   GfsBox * newbox1 = GFS_DOUBLE_TO_POINTER (GFS_VARIABLE (neighbor, p->newboxp->i));
	   g_assert (newbox1);
	   GFS_BOUNDARY_PERIODIC (newbox->neighbor[d])->matching = newbox1;
	   if (!node_is_linked (GTS_GNODE (newbox1), GTS_GNODE (newbox), 
				FTT_OPPOSITE_DIRECTION (d))) {
	     GfsGEdge * edge = GFS_GEDGE (gts_gedge_new (GTS_GRAPH (p->domain)->edge_class,
							 GTS_GNODE (newbox), 
							 GTS_GNODE (newbox1)));
	     edge->d = d;
	   }
	 }
	 else if (newbox->neighbor[d] == NULL) {
	   FttCell * neighbor = ftt_cell_neighbor (child.c[i], d);

	   if (neighbor) {
	     GfsBox * newbox1 = GFS_DOUBLE_TO_POINTER (GFS_VARIABLE (neighbor, p->newboxp->i));
	     FttDirection od = FTT_OPPOSITE_DIRECTION (d);
	     GfsGEdge * edge;

	     g_assert (newbox1);
	     newbox->neighbor[d] = GTS_OBJECT (newbox1);
	     g_assert (newbox1->neighbor[od] == NULL);
	     newbox1->neighbor[od] = GTS_OBJECT (newbox);
	     edge = GFS_GEDGE (gts_gedge_new (GTS_GRAPH (p->domain)->edge_class,
					      GTS_GNODE (newbox), 
					      GTS_GNODE (newbox1)));
	     edge->d = d;
	   }
	 }
    }
}

static void box_destroy (GfsBox * box, GfsVariable * newboxp)
{
  GfsBox * newbox[FTT_CELLS];
  FttCellChildren child;
  guint i;

  ftt_cell_children (box->root, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i])
      newbox[i] = GFS_DOUBLE_TO_POINTER (GFS_VARIABLE (child.c[i], newboxp->i));
    else
      newbox[i] = NULL;

  ftt_cell_destroy_root (box->root, &child, (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
  box->root = NULL;
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
      newbox[i]->root = child.c[i];
      FTT_ROOT_CELL (newbox[i]->root)->parent = newbox[i];
    }

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
  SplitPar p;

  g_return_if_fail (domain != NULL);

  p.newboxp = gfs_temporary_variable (domain);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, 1,
  			   (FttCellTraverseFunc) gfs_cell_reset, p.newboxp);
  p.boxlist = NULL;
  p.bid = 2;
  p.pid = 0;
  p.one_box_per_pe = one_box_per_pe;
  p.domain = domain;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_split, &p);
  g_slist_foreach (p.boxlist, (GFunc) box_link, &p);
  g_slist_foreach (p.boxlist, (GFunc) box_destroy, p.newboxp);
  g_slist_free (p.boxlist);
  gts_object_destroy (GTS_OBJECT (p.newboxp));

  gfs_domain_match (domain);
  domain->rootlevel++;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) get_ref_pos, &domain->refpos);
}

typedef struct {
  FttVector target;
  gint max_depth;
  FttCell * cell;
} LocateArgs;

static void box_locate (GfsBox * box, LocateArgs * a)
{
  if (a->cell == NULL)
    a->cell = ftt_cell_locate (box->root, a->target, a->max_depth);
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
  LocateArgs a;

  g_return_val_if_fail (domain != NULL, NULL);

  a.target = target;
  a.max_depth = max_depth;
  a.cell = NULL;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_locate, &a);

  return a.cell;
}

static void box_boundary_locate (GfsBox * box, LocateArgs * a)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (a->cell == NULL && box->neighbor[d] && GFS_IS_BOUNDARY (box->neighbor[d])) {
      FttCell * cell = ftt_cell_locate (GFS_BOUNDARY (box->neighbor[d])->root, a->target, 
					a->max_depth);
      if (cell && GFS_CELL_IS_BOUNDARY (cell))
	a->cell = cell;
    }
}

/**
 * gfs_domain_boundary_locate:
 * @domain: a #GfsDomain.
 * @target: position of the point to look for.
 * @max_depth: maximum depth to consider (-1 means no restriction).
 *
 * Locates the cell of the boundary of @domain containing @target.
 *
 * Returns: a #FttCell of the boundary of @domain containing the
 * point defined by @target or %NULL if @target is not contained in
 * any cell of the boundary of @domain.  
 */
FttCell * gfs_domain_boundary_locate (GfsDomain * domain,
				      FttVector target,
				      gint max_depth)
{
  LocateArgs a;

  g_return_val_if_fail (domain != NULL, NULL);

  a.target = target;
  a.max_depth = max_depth;
  a.cell = NULL;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_boundary_locate, &a);

  return a.cell;
}

static void box_distance2 (GfsBox * box, GPtrArray * a)
{
  g_ptr_array_add (a, box);
}

static void bubble_sort (GPtrArray * a, gdouble * d)
{
  guint i, j;

  for (i = 0; i < a->len - 1; i++)
    for (j = 0; j < a->len - 1 - i; j++)
      if (d[j+1] < d[j]) {
	gdouble tmp = d[j];
	gpointer data = a->pdata[j];
	d[j] = d[j+1];
	d[j+1] = tmp;
	a->pdata[j] = a->pdata[j+1];
	a->pdata[j+1] = data;
      }
}

/**
 * gfs_domain_cell_point_distance2:
 * @domain: a #GfsDomain.
 * @p: a #GtsPoint.
 * @distance2: the squared distance function.
 * @data: user data to pass to @distance2.
 * @closest: where to return the closest cell or %NULL.
 *
 * For non-leafs cells @distance2 must return a lower-bound for the
 * minimum distance (using for example ftt_cell_point_distance2_min()).
 *
 * Returns: the square of the minimum distance measured according to
 * @distance2 between @p and a leaf cell of @domain.
 */
gdouble gfs_domain_cell_point_distance2 (GfsDomain * domain,
					 GtsPoint * p,
					 gdouble (* distance2) (FttCell *, GtsPoint *, gpointer),
					 gpointer data,
					 FttCell ** closest)
{
  gdouble dmin = G_MAXDOUBLE;
  GPtrArray * a;
  gdouble * d;
  guint i;

  g_return_val_if_fail (domain != NULL, dmin);
  g_return_val_if_fail (p != NULL, dmin);
  g_return_val_if_fail (distance2 != NULL, dmin);

  if (closest)
    *closest = NULL;
  a = g_ptr_array_new ();
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_distance2, a);
  d = g_malloc (sizeof (gdouble)*a->len);
  for (i = 0; i < a->len; i++)
    d[i] = (* distance2) (GFS_BOX (a->pdata[i])->root, p, data);
  bubble_sort (a, d);
  for (i = 0; i < a->len; i++)
    if (d[i] < dmin)
      ftt_cell_point_distance2_internal (GFS_BOX (a->pdata[i])->root, p, d[i],
					 distance2, data, closest, &dmin);
  g_free (d);
  g_ptr_array_free (a, TRUE);
  return dmin;
}

/**
 * gfs_domain_advect_point:
 * @domain: a #GfsDomain.
 * @p: a #FttVector.
 * @dt: the time step.
 *
 * Updates the coordinates of point @p at time t + @dt using the
 * velocity field defined by @domain.
 *
 * If @p is not contained within @domain, the coordinates are unchanged.
 */
void gfs_domain_advect_point (GfsDomain * domain, 
			      FttVector * p,
			      gdouble dt)
{
  FttCell * cell;
  FttVector p0, p1;
  FttComponent c;
  GfsVariable ** u;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (p != NULL);

  p0 = p1 = *p;
  cell = gfs_domain_locate (domain, p0, -1);
  if (cell == NULL)
    return;
  u = gfs_domain_velocity (domain);
  for (c = 0; c < FTT_DIMENSION; c++)
    (&p1.x)[c] += dt*gfs_interpolate (cell, p0, u[c])/2.;
  cell = gfs_domain_locate (domain, p1, -1);
  if (cell == NULL)
    return;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&p->x)[c] += dt*gfs_interpolate (cell, p1, u[c]);
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
  gfs_all_reduce (domain, n, MPI_UNSIGNED, MPI_SUM);
  return n;
}

typedef struct {
  gdouble cfl;
  GfsVariable ** v;
} CflData;

static void minimum_mac_cfl (FttCellFace * face, CflData * p)
{
  gdouble un = GFS_STATE (face->cell)->f[face->d].un;
  if (un != 0.) {
    gdouble cflu = ftt_cell_size (face->cell)/fabs (un);
    if (cflu*cflu < p->cfl)
      p->cfl = cflu*cflu;
  }
  FttComponent c = face->d/2;
  if (p->v[c]->sources) {
    gdouble g = 0.;
    GSList * i = GTS_SLIST_CONTAINER (p->v[c]->sources)->items;
    while (i) {
      GfsSourceGeneric * s = i->data;
      if (s->face_value)
	g += (* s->face_value) (s, face, p->v[c]);
      i = i->next;
    }
    if (g != 0.) {
      gdouble cflg = 2.*ftt_cell_size (face->cell)/fabs (g);
      if (cflg < p->cfl)
	p->cfl = cflg;
    }
  }
}

static void minimum_cfl (FttCell * cell, CflData * p)
{
  gdouble size = ftt_cell_size (cell);
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++) {
    if (GFS_VARIABLE (cell, p->v[c]->i) != 0.) {
      gdouble cflu = size/fabs (GFS_VARIABLE (cell, p->v[c]->i));

      if (cflu*cflu < p->cfl)
	p->cfl = cflu*cflu;
    }
    if (p->v[c]->sources) {
      gdouble g = gfs_variable_mac_source (p->v[c], cell);

      if (g != 0.) {
	gdouble cflg = 2.*size/fabs (g);

	if (cflg < p->cfl)
	  p->cfl = cflg;
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
  CflData p;

  g_return_val_if_fail (domain != NULL, 0.);

  p.cfl = G_MAXDOUBLE;
  p.v = gfs_domain_velocity (domain);
  gfs_domain_face_traverse (domain, FTT_XYZ, FTT_PRE_ORDER, flags, max_depth, 
			    (FttFaceTraverseFunc) minimum_mac_cfl, &p);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			    (FttCellTraverseFunc) minimum_cfl, &p);
  gfs_all_reduce (domain, p.cfl, MPI_DOUBLE, MPI_MIN);
  return sqrt (p.cfl);
}

/**
 * gfs_cell_init:
 * @cell: a #FttCell.
 * @domain: a #GfsDomain containing @cell.
 *
 * Allocates the memory for fluid state data associated to @cell or its children.
 */
void gfs_cell_init (FttCell * cell, GfsDomain * domain)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (domain != NULL);

  if (FTT_CELL_IS_LEAF (cell)) {
    g_return_if_fail (cell->data == NULL);
    cell->data = g_malloc0 (gfs_domain_variables_size (domain));
  }
  else {
    FttCellChildren child;
    guint n;

    ftt_cell_children (cell, &child);
    for (n = 0; n < FTT_CELLS; n++) {
      g_return_if_fail (child.c[n]->data == NULL);
      child.c[n]->data = g_malloc0 (gfs_domain_variables_size (domain));
    }
  }
}

/**
 * gfs_cell_reinit:
 * @cell: a #FttCell.
 * @domain: a #GfsDomain containing @cell.
 *
 * Re-allocates the memory for fluid state data associated to @cell.
 */
void gfs_cell_reinit (FttCell * cell, GfsDomain * domain)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (cell->data != NULL);
  g_return_if_fail (domain != NULL);

  cell->data = g_realloc (cell->data, gfs_domain_variables_size (domain));
}

/**
 * gfs_cell_fine_init:
 * @parent: a #FttCell.
 * @domain: a #GfsDomain containing @parent.
 *
 * Initialises the children of @parent.
 */
void gfs_cell_fine_init (FttCell * parent, GfsDomain * domain)
{
  GSList * i;

  g_return_if_fail (parent != NULL);
  g_return_if_fail (!FTT_CELL_IS_LEAF (parent));
  g_return_if_fail (domain != NULL);

  gfs_cell_init (parent, domain);

  if (!GFS_CELL_IS_BOUNDARY (parent) && GFS_IS_MIXED (parent))
    gfs_solid_coarse_fine (parent);

  i = domain->variables;
  while (i) {
    GfsVariable * v = i->data;
  
    (* v->coarse_fine) (parent, v);
    i = i->next;
  }
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
    memcpy (to->data, from->data, gfs_domain_variables_size (domain));
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
 * @variables: the list of #GfsVariable to be written.
 *
 * Writes in @fp the fluid data associated with @cell and described by
 * @variables. This function is generally used in association with
 * ftt_cell_write().  
 */
void gfs_cell_write (const FttCell * cell, FILE * fp,
		     GSList * variables)
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
    fprintf (fp, " %g", GFS_VARIABLE (cell, GFS_VARIABLE1 (variables->data)->i));
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
  GSList * i;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (fp != NULL);
  g_return_if_fail (domain != NULL);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (solid->s[0])");
    return;
  }
  s0 = atof (fp->token->str);
  if (s0 < 0. && s0 != -1.) {
    gts_file_error (fp, "solid->s[0] must be positive");
    return;
  }
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

  i = domain->variables_io;
  while (i) {
    GfsVariable * v = i->data;

    if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
      gts_file_error (fp, "expecting a number (%s)", v->name);
      return;
    }
    GFS_VARIABLE (cell, v->i) = atof (fp->token->str);
    gts_file_next_token (fp);
    i = i->next;
  }
}

/**
 * gfs_cell_write_binary:
 * @cell: a #FttCell.
 * @fp: a file pointer.
 * @variables: the list of #GfsVariable to be written.
 *
 * Writes in @fp the fluid data associated with @cell and described by
 * @variables. This function is generally used in association with
 * ftt_cell_write_binary().
 */
void gfs_cell_write_binary (const FttCell * cell, FILE * fp,
			    GSList * variables)
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
    gdouble a = GFS_VARIABLE (cell, GFS_VARIABLE1 (variables->data)->i);
    fwrite (&a, sizeof (gdouble), 1, fp);
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
  GSList * i;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (fp != NULL);
  g_return_if_fail (domain != NULL);

  if (gts_file_read (fp, &s0, sizeof (gdouble), 1) != 1) {
    gts_file_error (fp, "expecting a number (solid->s[0])");
    return;
  }
  if (s0 < 0. && s0 != -1.) {
    gts_file_error (fp, "solid->s[0] must be positive");
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

  i = domain->variables_io;
  while (i) {
    GfsVariable * v = i->data;
    gdouble a;

    if (gts_file_read (fp, &a, sizeof (gdouble), 1) != 1) {
      gts_file_error (fp, "expecting a number (%s)", v->name);
      return;
    }
    GFS_VARIABLE (cell, v->i) = a;
    i = i->next;
  }
}

static void box_realloc (GfsBox * box, GfsDomain * domain)
{
  FttDirection d;

  ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		     (FttCellTraverseFunc) gfs_cell_reinit, domain);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d]))
      ftt_cell_traverse (GFS_BOUNDARY (box->neighbor[d])->root, 
			 FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			 (FttCellTraverseFunc) gfs_cell_reinit, domain);
}

/**
 * gfs_domain_alloc:
 * @domain: a #GfsDomain.
 *
 * Returns: the index of a memory location newly allocated for each
 * cell of @domain.
 */
guint gfs_domain_alloc (GfsDomain * domain)
{
  guint i = 0;

  g_return_val_if_fail (domain != NULL, -1);

  while (i < domain->allocated->len && g_array_index (domain->allocated, gboolean, i))
    i++;
  if (i == domain->allocated->len) {
    g_array_set_size (domain->allocated, domain->allocated->len + 1);
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_realloc, domain);
  }
  g_array_index (domain->allocated, gboolean, i) = TRUE;
  return i;
}

/**
 * gfs_domain_free:
 * @domain: a #GfsDomain.
 * @i: a memory location index previously allocated using gfs_domain_alloc().
 *
 * Frees the memory location of @domain defined by @i.
 */
void gfs_domain_free (GfsDomain * domain, guint i)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (i < domain->allocated->len);
  g_return_if_fail (g_array_index (domain->allocated, gboolean, i));

  g_array_index (domain->allocated, gboolean, i) = FALSE;
}

/**
 * gfs_domain_add_variable:
 * @domain: a #GfsDomain.
 * @name: the name of the variable to add or %NULL.
 * @description: the variable description or %NULL.
 *
 * Adds a new variable @name to @domain.
 *
 * Returns: the new variable or %NULL if a variable with the same name
 * already exists.  
 */
GfsVariable * gfs_domain_add_variable (GfsDomain * domain,
				       const gchar * name,
				       const gchar * description)
{
  GfsVariable * v;

  g_return_val_if_fail (domain != NULL, NULL);

  if ((v = gfs_variable_new (gfs_variable_class (), domain, name, description)) == NULL)
    return NULL;
  domain->variables = g_slist_append (domain->variables, v);
  return v;
}

/**
 * gfs_domain_get_or_add_variable:
 * @domain: a #GfsDomain.
 * @name: the name of the variable to add or get.
 * @description: the variable description or %NULL.
 *
 * Adds a new variable @name to @domain or returns the variable of
 * @domain with the same name. In either case the description of the
 * variable name is set to @description.
 *
 * Returns: the new or already existing variable or %NULL if @name is a
 * reserved variable name.
 */
GfsVariable * gfs_domain_get_or_add_variable (GfsDomain * domain,
					      const gchar * name,
					      const gchar * description)
{
  GfsVariable * v;

  g_return_val_if_fail (domain != NULL, NULL);
  g_return_val_if_fail (name != NULL, NULL);

  v = gfs_variable_from_name (domain->variables, name);
  if (v != NULL) {
    if (v->description)
      g_free (v->description);
    v->description = description ? g_strdup (description) : NULL;
  }
  else
    v = gfs_domain_add_variable (domain, name, description);
  return v;
}

static void add_pressure_force (FttCell * cell, gpointer * data)
{
  gdouble * f = data[0];
  gdouble * m = data[1];
  gdouble * r = &GFS_STATE (cell)->solid->ca.x;
  GfsVariable * p = data[2];
  FttVector ff, mm;
  FttComponent c;

  gfs_pressure_force (cell, p, &ff);
  gts_vector_cross (&mm.x, r, &ff.x);
  for (c = 0; c < 3; c++) {
    f[c] += (&ff.x)[c];
    m[c] += (&mm.x)[c];
  }
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
  gdouble * f = data[0];
  gdouble * m = data[1];
  GfsVariable * v = data[2];
  GfsSourceDiffusion * d = data[3];
  gdouble D;
  GfsSolidVector * s = GFS_STATE (cell)->solid;
  gdouble * r = &s->ca.x;
  FttVector ff, mm, n, g;
  FttComponent c;

  g_assert (((cell)->flags & GFS_FLAG_DIRICHLET) != 0);
  gfs_cell_dirichlet_gradient (cell, v->i, -1, s->fv, &g);

  D = - gfs_source_diffusion_cell (d, cell);
  n.x = s->s[1] - s->s[0];
  n.y = s->s[3] - s->s[2];
#if FTT_2D
  ff.z = 0.;
  switch (v->component) {
  case FTT_X:
    ff.x = D*(2.*g.x*n.x + g.y*n.y);
    ff.y = D*g.y*n.x;
    break;
  case FTT_Y:
    ff.x = D*g.x*n.y;
    ff.y = D*(2.*g.y*n.y + g.x*n.x);
    break;
  default:
    g_assert_not_reached ();
  }
#else /* 3D */
  n.z = s->s[5] - s->s[4];
  D *= ftt_cell_size (cell);
  switch (v->component) {
  case FTT_X:
    ff.x = D*(2.*g.x*n.x + g.y*n.y + g.z*n.z);
    ff.y = D*g.y*n.x;
    ff.z = D*g.z*n.x;
    break;
  case FTT_Y:
    ff.y = D*(2.*g.y*n.y + g.x*n.x + g.z*n.z);
    ff.x = D*g.x*n.y;
    ff.z = D*g.z*n.y;
    break;
  case FTT_Z:
    ff.z = D*(2.*g.z*n.z + g.x*n.x + g.y*n.y);
    ff.x = D*g.x*n.z;
    ff.y = D*g.y*n.z;
    break;
  default:
    g_assert_not_reached ();
  }
#endif /* 3D */
  gts_vector_cross (&mm.x, r, &ff.x);
  for (c = 0; c < 3; c++) {
    f[c] += (&ff.x)[c];
    m[c] += (&mm.x)[c];
  }
}

/**
 * gfs_domain_solid_force:
 * @domain: a #GfsDomain.
 * @pf: a #FttVector.
 * @vf: a #FttVector.
 * @pm: a #FttVector.
 * @vm: a #FttVector.
 *
 * Fills @pf and @vf (resp. @pm and @vm) with the components of the
 * net pressure and viscous forces (resp. pressure and viscous
 * moments) applied by the fluid on the solid surface embbeded in
 * @domain.
 *
 * The reference point for the moments is the origin of the coordinate system.
 */
void gfs_domain_solid_force (GfsDomain * domain, 
			     FttVector * pf,
			     FttVector * vf,
			     FttVector * pm,
			     FttVector * vm)
{
  FttComponent c;
  GfsVariable ** v;
  gpointer data[3];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (pf != NULL);
  g_return_if_fail (vf != NULL);
  g_return_if_fail (pm != NULL);
  g_return_if_fail (vm != NULL);

  if (GFS_IS_AXI (domain))
    g_assert_not_implemented ();

  pf->x = pf->y = pf->z = 0.;
  pm->x = pm->y = pm->z = 0.;
  data[0] = pf;
  data[1] = pm;
  data[2] = gfs_variable_from_name (domain->variables, "P");
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseFunc) add_pressure_force, data);

  vf->x = vf->y = vf->z = 0.;
  vm->x = vm->y = vm->z = 0.;
  v = gfs_domain_velocity (domain);
  for (c = 0; c < FTT_DIMENSION; c++) {
    GfsSourceDiffusion * D = source_diffusion (v[c]);

    if (D) {
      gpointer data[4];

      gfs_domain_surface_bc (domain, v[c]);
      data[0] = vf;
      data[1] = vm;
      data[2] = v[c];
      data[3] = D;
      gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
				 (FttCellTraverseFunc) add_viscous_force, data);
    }
  }
}

#define THRESHOLD 1e-4

static void tag_cell_fraction (GtsFifo * fifo,
			       FttCell * cell,
			       GfsVariable * c, GfsVariable * v,
			       guint tag)
{
  FttDirection d;
  FttCellNeighbors n;

  g_assert (FTT_CELL_IS_LEAF (cell));
  ftt_cell_neighbors (cell, &n);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (n.c[d] && GFS_VALUE (n.c[d], v) == 0. && GFS_VALUE (n.c[d], c) > THRESHOLD) {
      if (FTT_CELL_IS_LEAF (n.c[d])) {
	GFS_VALUE (n.c[d], v) = tag;
	gts_fifo_push (fifo, n.c[d]);
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
	  if (child.c[i] && GFS_VALUE (child.c[i], v) == 0. &&
	      GFS_VALUE (child.c[i], c) > THRESHOLD) {
	    GFS_VALUE (child.c[i], v) = tag;
	    gts_fifo_push (fifo, child.c[i]);
	  }
      }
    }
}

typedef struct {
  GfsVariable * v, * c;
  guint tag;
} TagPar;

static void tag_new_fraction_region (FttCell * cell, TagPar * p)
{
  if (GFS_VALUE (cell, p->v) == 0. && GFS_VALUE (cell, p->c) > THRESHOLD) {
    GtsFifo * fifo = gts_fifo_new ();

    GFS_VALUE (cell, p->v) = ++p->tag;
    gts_fifo_push (fifo, cell);
    while ((cell = gts_fifo_pop (fifo)))
      tag_cell_fraction (fifo, cell, p->c, p->v, p->tag);
    gts_fifo_destroy (fifo);
  }
}

/**
 * gfs_domain_tag_droplets:
 * @domain: a #GfsDomain.
 * @c: the volume fraction.
 * @tag: a #GfsVariable.
 *
 * Fills the @tag variable of the cells of @domain with the (strictly
 * positive) index of the droplet they belong to. The cells belonging
 * to the background phase have an index of zero.
 *
 * Note that the volume fraction @c must be defined on all levels.
 *
 * Returns: the number of droplets.
 */
guint gfs_domain_tag_droplets (GfsDomain * domain,
			       GfsVariable * c,
			       GfsVariable * tag)
{
  /* fixme: this function may not work as expected for parallel domain
     and/or periodic boundaries: droplets sitting on PE boundaries
     will be seen as two independent droplets... */

  g_return_val_if_fail (domain != NULL, 0);
  g_return_val_if_fail (c != NULL, 0);
  g_return_val_if_fail (tag != NULL, 0);

  TagPar p;
  p.c = c;
  p.v = tag;
  p.tag = 0;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, tag);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) tag_new_fraction_region, &p);
  return p.tag;
}

typedef struct {
  GfsVariable * tag, * c;
  guint * sizes;
  guint n, min;
} RemoveDropletsPar;

static void compute_droplet_size (FttCell * cell, RemoveDropletsPar * p)
{
  guint i = GFS_VALUE (cell, p->tag);
  if (i > 0)
    p->sizes[i - 1]++;
}

static void reset_small_fraction (FttCell * cell, RemoveDropletsPar * p)
{
  guint i = GFS_VALUE (cell, p->tag);
  if (i > 0 && p->sizes[i - 1] < p->min)
    GFS_VALUE (cell, p->c) = 0.;
}

static int greater (const void * a, const void * b)
{
  return *((guint *)a) > *((guint *)b) ? -1 : 1;
}

/**
 * gfs_domain_remove_droplets:
 * @domain: a #GfsDomain.
 * @c: a #GfsVariable.
 * @v: a #GfsVariable.
 * @min: the minimum size (in cells) of the droplets.
 *
 * Resets the @v variable of all the droplets (defined by the @c
 * variable) smaller than @min cells if @min is positive, or all the
 * droplets but the -$min largest ones if @min is negative.
 */
void gfs_domain_remove_droplets (GfsDomain * domain,
				 GfsVariable * c,
				 GfsVariable * v,
				 gint min)
{
  RemoveDropletsPar p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (c != NULL);
  g_return_if_fail (v != NULL);

  p.c = c;
  p.tag = gfs_temporary_variable (domain);
  p.n = gfs_domain_tag_droplets (domain, c, p.tag);
  if (p.n > 0 && -min < (gint) p.n) {
    p.sizes = g_malloc0 (p.n*sizeof (guint));
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) compute_droplet_size, &p);
    if (min >= 0)
      p.min = min;
    else {
      guint * tmp = g_malloc (p.n*sizeof (guint));
      memcpy (tmp, p.sizes, p.n*sizeof (guint));
      qsort (tmp, p.n, sizeof (guint), greater);
      /* fixme: this won't work for parallel jobs */
      p.min = tmp[-1 - min];
      g_free (tmp);
    }
    p.c = v;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) reset_small_fraction, &p);
    g_free (p.sizes);
  }
  gts_object_destroy (GTS_OBJECT (p.tag));
}

static void tag_cell (FttCell * cell, GfsVariable * v, guint tag, guint * size)
{
  FttDirection d;
  FttCellNeighbors n;
  GfsSolidVector * solid = GFS_STATE (cell)->solid;

  g_assert (FTT_CELL_IS_LEAF (cell));
  GFS_VARIABLE (cell, v->i) = tag;
  (*size)++;
  ftt_cell_neighbors (cell, &n);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (n.c[d] && GFS_VARIABLE (n.c[d], v->i) == 0. &&
	!GFS_CELL_IS_BOUNDARY (n.c[d]) &&
	(!solid || solid->s[d] > 0.)) {
      if (FTT_CELL_IS_LEAF (n.c[d]))
	tag_cell (n.c[d], v, tag, size);
      else {
	FttCellChildren child;
	FttDirection od = FTT_OPPOSITE_DIRECTION (d);
	guint i, j;
	
	j = ftt_cell_children_direction (n.c[d], od, &child);
	for (i = 0; i < j; i++)
	  if (child.c[i] && GFS_VARIABLE (child.c[i], v->i) == 0. &&
	      (!GFS_IS_MIXED (child.c[i]) || GFS_STATE (child.c[i])->solid->s[od] > 0.))
	    tag_cell (child.c[i], v, tag, size);
      }
    }
}

static void tag_new_region (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];

  if (GFS_VARIABLE (cell, v->i) == 0.) {
    GArray * sizes = data[1];
    guint size = 0;

    tag_cell (cell, v, sizes->len + 1, &size);
    g_array_append_val (sizes, size);
  }
}

static gboolean remove_small (FttCell * cell, gpointer * data)
{
  if (FTT_CELL_IS_LEAF (cell)) {
    GArray * sizes = data[0];
    GfsVariable * v = data[5];
    guint * min = data[1], i = GFS_VARIABLE (cell, v->i) - 1.;

    g_assert (GFS_VARIABLE (cell, v->i) > 0.);
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
  GfsVariable * v;

  g_return_if_fail (domain != NULL);

  v = gfs_temporary_variable (domain);
  sizes = g_array_new (FALSE, FALSE, sizeof (guint));
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, v);
  dat[0] = v;
  dat[1] = sizes;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) tag_new_region, dat);
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
  dat[5] = v;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) remove_small_box, dat);
  g_array_free (sizes, TRUE);
  gts_object_destroy (GTS_OBJECT (v));
  if (changed)
    gfs_domain_match (domain);
}

static gboolean tag_speck (FttCell * cell, GfsVariable * v)
{
  if (GFS_VARIABLE (cell, v->i) == 0.) {
    FttDirection d;
    FttCellNeighbors n;
    GfsSolidVector * solid = GFS_STATE (cell)->solid;
    
    g_assert (FTT_CELL_IS_LEAF (cell));
    ftt_cell_neighbors (cell, &n);
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (!n.c[d])
	return FALSE;
    GFS_VARIABLE (cell, v->i) = 1.;
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (GFS_VARIABLE (n.c[d], v->i) == 0. && 
	  !GFS_CELL_IS_BOUNDARY (n.c[d]) &&
	  solid->s[d] > 0. && solid->s[d] < 1.) {
	g_assert (GFS_IS_MIXED (n.c[d]));
	if (FTT_CELL_IS_LEAF (n.c[d])) {
	  if (!tag_speck (n.c[d], v)) {
	    GFS_VARIABLE (cell, v->i) = 0.;
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
	    if (!child.c[i] || (GFS_VARIABLE (child.c[i], v->i) == 0 && 
				GFS_IS_MIXED (child.c[i]) &&
				!tag_speck (child.c[i], v))) {
	      GFS_VARIABLE (cell, v->i) = 0.;
	      return FALSE;
	    }
	}
      }
  }
  return TRUE;
}

static void fill_speck (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];

  if (GFS_VARIABLE (cell, v->i) == 1.) {
    gboolean * changed = data[1];
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
  GfsVariable * v;
  gpointer data[2];

  g_return_if_fail (domain != NULL);

  v = gfs_temporary_variable (domain);
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL,
			     (FttCellTraverseFunc) gfs_cell_reset, v);
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseFunc) tag_speck, v);
  data[0] = v;
  data[1] = &changed;
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseFunc) fill_speck, data);
  gts_object_destroy (GTS_OBJECT (v));
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
  t->start = gfs_clock_elapsed (domain->timer);
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
  end = gfs_clock_elapsed (domain->timer);
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

/**
 * gfs_domain_add_derived_variable:
 * @domain: a #GfsDomain.
 * @info: the #GfsDerivedVariableInfo.
 *
 * Adds a derived variable described by @info to @domain.
 *
 * Returns: the #GfsDerivedVariable if the variable was successfully
 * added to @domain or %NULL if a variable with the same name already
 * exists.
 */
GfsDerivedVariable * gfs_domain_add_derived_variable (GfsDomain * domain, 
						      GfsDerivedVariableInfo info)
{
  GfsDerivedVariable * v;

  g_return_val_if_fail (domain != NULL, NULL);

  if (gfs_variable_from_name (domain->variables, info.name) ||
      gfs_derived_variable_from_name (domain->derived_variables, info.name))
    return NULL;
  v = GFS_DERIVED_VARIABLE (gts_object_new (GTS_OBJECT_CLASS (gfs_derived_variable_class ())));
  v->name = g_strdup (info.name);
  v->description = info.description ? g_strdup (info.description) : NULL;
  v->func = info.func;
  v->data = info.data;
  domain->derived_variables = g_slist_prepend (domain->derived_variables, v);
  GTS_OBJECT (v)->reserved = domain;
  return v;
}

/**
 * gfs_domain_remove_derived_variable:
 * @domain: a #GfsDomain.
 * @name: the name of a #GfsDerivedVariable.
 *
 * Removes derived variable @name from @domain.
 *
 * Returns: %TRUE if the variable was successfully removed from @domain or
 * %FALSE if a derived variable with the this name does not exist.
 */
gboolean gfs_domain_remove_derived_variable (GfsDomain * domain, const gchar * name)
{
  GSList * i;
  
  g_return_val_if_fail (domain != NULL, FALSE);
  g_return_val_if_fail (name != NULL, FALSE);

  i = domain->derived_variables;
  while (i) {
    GfsDerivedVariable * u = i->data;

    if (!strcmp (u->name, name)) {
      gts_object_destroy (GTS_OBJECT (u));
      domain->derived_variables = g_slist_remove_link (domain->derived_variables, i);
      g_slist_free (i);
      return TRUE;
    }
    i = i->next;
  }
  return FALSE;
}

typedef struct {
  FttDirection d;
  GfsFunction * f;
  GfsVariable * v;
} SumData;

static gdouble product (FttCell * cell, GfsFunction * f)
{
  GfsSolidVector * solid = GFS_STATE (cell)->solid;
  return ftt_cell_volume (cell)*(solid ? solid->a : 1.)*gfs_function_value (f, cell);
}

static void sum (FttCell * cell, SumData * data)
{
  FttCell * n = ftt_cell_neighbor (cell, data->d);
  GfsSolidVector * solid = GFS_STATE (cell)->solid;

  if (!n || GFS_CELL_IS_BOUNDARY (n) || (solid && solid->s[data->d] == 0.)) {
    gdouble s = 0.;

    n = cell;
    do {
      /* fixme: does not work if the resolution varies along data->d */
      g_assert (ftt_cell_level (n) == ftt_cell_level (cell));
      s += product (n, data->f);
      GFS_VARIABLE (n, data->v->i) = s;
      n = ftt_cell_neighbor (n, FTT_OPPOSITE_DIRECTION (data->d));
    } while (n && !GFS_CELL_IS_BOUNDARY (n) && 
	     (!GFS_IS_MIXED (n) || GFS_STATE (n)->solid->s[data->d] > 0.));
  }
}

/**
 * gfs_domain_sum:
 * @domain: a #GfsDomain.
 * @d: the #FttDirection.
 * @f: a #GfsFunction.
 * @v: a #GfsVariable.
 *
 * Fills variable @v of each cell of @domain with the sum in direction
 * @d of the volume-weighted function @f.
 */
void gfs_domain_sum (GfsDomain * domain, FttDirection d, GfsFunction * f, GfsVariable * v)
{
  SumData data;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (d >= 0 && d < FTT_NEIGHBORS);
  g_return_if_fail (f != NULL);
  g_return_if_fail (v != NULL);

  data.d = d;
  data.f = f;
  data.v = v;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) sum, &data);
}

static void filter (FttCell * cell, gpointer * data)
{
  FttDirection d[4*(FTT_DIMENSION - 1)][FTT_DIMENSION] = {
#if FTT_2D
    {FTT_RIGHT, FTT_TOP}, {FTT_RIGHT, FTT_BOTTOM}, {FTT_LEFT, FTT_TOP}, {FTT_LEFT, FTT_BOTTOM}
#else
    {FTT_RIGHT, FTT_TOP, FTT_FRONT}, {FTT_RIGHT, FTT_BOTTOM, FTT_FRONT}, 
    {FTT_LEFT, FTT_TOP, FTT_FRONT}, {FTT_LEFT, FTT_BOTTOM, FTT_FRONT},
    {FTT_RIGHT, FTT_TOP, FTT_BACK}, {FTT_RIGHT, FTT_BOTTOM, FTT_BACK}, 
    {FTT_LEFT, FTT_TOP, FTT_BACK}, {FTT_LEFT, FTT_BOTTOM, FTT_BACK}
#endif
  };
  guint i;
  gdouble val = 0.;
  GfsVariable * a = data[0];
  GfsVariable * b = data[1];

  for (i = 0; i < 4*(FTT_DIMENSION - 1); i++)
    val += gfs_cell_corner_value (cell, d[i], a, -1);
  GFS_VARIABLE (cell, b->i) = val/(4*(FTT_DIMENSION - 1));
}

/**
 * gfs_domain_filter:
 * @domain: a #GfsDomain.
 * @v: a #GfsVariable.
 * @fv: the filtered variable or %NULL.
 *
 * Apply a "corner-averaging" filter to variable @v on all leaf cells
 * of @domain.
 *
 * If @fv is %NULL, @v is replaced by its filtered value.
 */
void gfs_domain_filter (GfsDomain * domain, GfsVariable * v, GfsVariable * fv)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);

  gpointer data[2];
  data[0] = v;
  data[1] = fv ? fv : gfs_temporary_variable (domain);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) filter, data);
  if (fv == NULL) {
    gfs_variables_swap (data[0], data[1]);
    gts_object_destroy (data[1]);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, v);
  }
  else
    gfs_domain_copy_bc (domain, FTT_TRAVERSE_LEAFS, -1, v, fv);
}
