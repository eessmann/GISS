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

#include "domain.h"
#include "mpi_boundary.h"
#include "adaptive.h"

/* #define DEBUG */

static guint tag_shift = 32767/FTT_NEIGHBORS;

#define TAG(boundary)           (tag_shift*(boundary)->d + (boundary)->box->id)
#define MATCHING_TAG(boundary)  (tag_shift*FTT_OPPOSITE_DIRECTION ((boundary)->d) +\
                                 GFS_BOUNDARY_MPI (boundary)->id)

static void boundary_mpi_destroy (GtsObject * object)
{
  GfsBoundaryMpi * boundary = GFS_BOUNDARY_MPI (object);

  g_array_free (boundary->sndbuf, TRUE);
  g_array_free (boundary->rcvbuf, TRUE);
  
  (* GTS_OBJECT_CLASS (gfs_boundary_mpi_class ())->parent_class->destroy) 
    (object);
}

static void boundary_mpi_read (GtsObject ** object, GtsFile * fp)
{
  boundary_mpi_destroy (*object);
}

static void center_mpi (FttCellFace * face, GfsBc * b)
{
  GfsBoundaryMpi * boundary_mpi = GFS_BOUNDARY_MPI (b->b);

  g_assert (boundary_mpi->sndcount < boundary_mpi->sndbuf->len);
  g_assert (ftt_face_type (face) == FTT_FINE_FINE);
  g_assert (!FTT_CELL_IS_LEAF (face->cell) || FTT_CELL_IS_LEAF (face->neighbor));
  g_array_index (boundary_mpi->sndbuf, gdouble, boundary_mpi->sndcount++) =
    GFS_VARIABLE (face->neighbor, b->v->i);
}

static void face_mpi (FttCellFace * face, GfsBc * b)
{
  GfsBoundaryMpi * boundary_mpi = GFS_BOUNDARY_MPI (b->b);

  g_assert (boundary_mpi->sndcount < boundary_mpi->sndbuf->len);
  g_array_index (boundary_mpi->sndbuf, gdouble, boundary_mpi->sndcount++) =
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v;
}

static void boundary_size (FttCell * cell, guint * count)
{
  (*count)++;
}

static void set_buffers_size (GfsBoundaryMpi * boundary)
{
  guint count = 0;

  ftt_cell_traverse (GFS_BOUNDARY (boundary)->root, 
		     FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) boundary_size, &count);
  g_array_set_size (boundary->rcvbuf, count);
  g_array_set_size (boundary->sndbuf, count);
}

static void boundary_tree (FttCell * cell, GfsBoundaryMpi * boundary)
{
  gdouble is_leaf = FTT_CELL_IS_LEAF (cell);

  if (boundary->sndcount == boundary->sndbuf->len)
    g_array_append_val (boundary->sndbuf, is_leaf);
  else
    g_array_index (boundary->sndbuf, gdouble, boundary->sndcount) = is_leaf;
  boundary->sndcount++;

  if (!is_leaf) {
    FttCellChildren child;
    guint i, n;

    n = ftt_cell_children_direction (cell, GFS_BOUNDARY (boundary)->d, &child);
    for (i = 0; i < n; i++) {
      gdouble is_destroyed = (child.c[i] == NULL);
      
      if (boundary->sndcount == boundary->sndbuf->len)
	g_array_append_val (boundary->sndbuf, is_destroyed);
      else
	g_array_index (boundary->sndbuf, gdouble, boundary->sndcount) = is_destroyed;
      boundary->sndcount++;
    }
  }
}

static void match (GfsBoundary * boundary)
{
  (* gfs_boundary_class ()->match) (boundary);

  g_assert (GFS_BOUNDARY_MPI (boundary)->sndcount == 0);
  ftt_cell_traverse (boundary->root,
		     FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		     (FttCellTraverseFunc) boundary_tree, boundary);
}

static void send (GfsBoundary * bb)
{
  GfsBoundaryMpi * boundary = GFS_BOUNDARY_MPI (bb);
  GfsDomain * domain = gfs_box_domain (bb->box);

  g_assert (boundary->sndcount <= boundary->sndbuf->len);
  if (GFS_BOUNDARY (boundary)->type == GFS_BOUNDARY_MATCH_VARIABLE) {
#ifdef DEBUG
fprintf (stderr, "%d send %d tag: %d\n",
	 domain->pid, 
	 boundary->process,
	 TAG (GFS_BOUNDARY (boundary)));
#endif
    MPI_Isend (&boundary->sndcount, 1, MPI_UNSIGNED,
	       boundary->process,
	       TAG (GFS_BOUNDARY (boundary)),
	       boundary->comm,
	       &(boundary->request[boundary->nrequest++]));
    gts_range_add_value (&domain->mpi_messages, sizeof (guint));
  }
#ifdef DEBUG
fprintf (stderr, "%d send %d tag: %d size: %d\n",
	 domain->pid, 
	 boundary->process,
	 TAG (GFS_BOUNDARY (boundary)),
	 boundary->sndcount);
#endif
  MPI_Isend (boundary->sndbuf->data, boundary->sndcount, MPI_DOUBLE,
	     boundary->process,
	     TAG (GFS_BOUNDARY (boundary)),
	     boundary->comm,
	     &(boundary->request[boundary->nrequest++]));
  gts_range_add_value (&domain->mpi_messages, 
                       sizeof (gdouble)*boundary->sndcount);
}

static void center_update (FttCell * cell,
			   GfsBoundaryMpi * boundary)
{
  g_assert (boundary->rcvcount < boundary->rcvbuf->len);
  GFS_VARIABLE (cell, GFS_BOUNDARY (boundary)->v->i) =
    g_array_index (boundary->rcvbuf, gdouble, boundary->rcvcount++);
}

static void face_update (FttCellFace * face,
			 GfsBoundaryMpi * boundary)
{
  g_assert (boundary->rcvcount < boundary->rcvbuf->len);
  GFS_STATE (face->cell)->f[face->d].v = 
    g_array_index (boundary->rcvbuf, gdouble, boundary->rcvcount++);
}

static void match_ignore (GfsBoundaryMpi * boundary)
{
  gboolean is_leaf;

  g_assert (boundary->rcvcount < boundary->rcvbuf->len);
  is_leaf = g_array_index (boundary->rcvbuf, gdouble, boundary->rcvcount++);

  if (!is_leaf) {
    gboolean is_destroyed[FTT_CELLS/2];
    guint i;

    for (i = 0; i < FTT_CELLS/2; i++) {
      g_assert (boundary->rcvcount < boundary->rcvbuf->len);
      is_destroyed[i] = g_array_index (boundary->rcvbuf, gdouble, boundary->rcvcount++);
    }
    for (i = 0; i < FTT_CELLS/2; i++)
      if (!is_destroyed[i])
	match_ignore (boundary);
  }
}

static void match_update (FttCell * cell,
			  GfsBoundaryMpi * boundary)
{
  gboolean is_leaf;

  g_assert (boundary->rcvcount < boundary->rcvbuf->len);
  is_leaf = g_array_index (boundary->rcvbuf, gdouble, boundary->rcvcount++);

  if (!is_leaf) {
    GfsDomain * domain = gfs_box_domain (GFS_BOUNDARY (boundary)->box);
    FttCellChildren child;
    gboolean is_destroyed[FTT_CELLS/2];
    guint i, n;

    if (FTT_CELL_IS_LEAF (cell)) {
      FttCell * neighbor = ftt_cell_neighbor (cell, GFS_BOUNDARY (boundary)->d);

      g_assert (neighbor);
      ftt_cell_refine_single (cell, (FttCellInitFunc) gfs_cell_fine_init, domain);
      if (FTT_CELL_IS_LEAF (neighbor))
	ftt_cell_refine_single (neighbor, (FttCellInitFunc) gfs_cell_fine_init, domain);
      /* what about solid fractions? */
      GFS_BOUNDARY (boundary)->changed = TRUE;
    }
    n = ftt_cell_children_direction (cell, GFS_BOUNDARY (boundary)->d, &child);
    for (i = 0; i < n; i++) {
      g_assert (boundary->rcvcount < boundary->rcvbuf->len);
      is_destroyed[i] = g_array_index (boundary->rcvbuf, gdouble, boundary->rcvcount++);
      if (is_destroyed[i] && child.c[i]) {
	ftt_cell_destroy (child.c[i], (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
	child.c[i] = NULL;
	GFS_BOUNDARY (boundary)->changed = TRUE;
      }
    }
    for (i = 0; i < n; i++)
      if (!is_destroyed[i]) {
	if (child.c[i])
	  match_update (child.c[i], boundary);
	else
	  match_ignore (boundary);
      }
  }
}

static void receive (GfsBoundary * bb,
		     FttTraverseFlags flags,
		     gint max_depth)
{
  GfsBoundaryMpi * boundary = GFS_BOUNDARY_MPI (bb);
  MPI_Status status;
  gint count;
#ifdef PROFILE_MPI
  GfsDomain * domain = gfs_box_domain (bb->box);
  gdouble start, end;

  start = MPI_Wtime ();
#endif /* PROFILE_MPI */

  if (GFS_BOUNDARY (boundary)->type == GFS_BOUNDARY_MATCH_VARIABLE) {
#ifdef DEBUG
fprintf (stderr, "%d wait %d %d match variable\n",
	 gfs_box_domain (bb->box)->pid,
	 boundary->process,
	 MATCHING_TAG (GFS_BOUNDARY (boundary)));
#endif
    MPI_Recv (&boundary->rcvcount, 1, MPI_UNSIGNED,
	      boundary->process,
	      MATCHING_TAG (GFS_BOUNDARY (boundary)),
	      boundary->comm,
	      &status);
#ifdef PROFILE_MPI
    end = MPI_Wtime ();
    gts_range_add_value (&domain->mpi_wait, end - start);
    start = MPI_Wtime ();
#endif /* PROFILE_MPI */
    if (boundary->rcvcount > boundary->rcvbuf->len)
      g_array_set_size (boundary->rcvbuf, boundary->rcvcount);
  }
  else
    boundary->rcvcount = boundary->sndcount;
#ifdef DEBUG
fprintf (stderr, "%d wait %d %d\n",
	 gfs_box_domain (bb->box)->pid,
	 boundary->process,
	 MATCHING_TAG (GFS_BOUNDARY (boundary)));
#endif
  g_assert (boundary->rcvcount <= boundary->rcvbuf->len);
  MPI_Recv (boundary->rcvbuf->data,
	    boundary->rcvcount,
	    MPI_DOUBLE,
	    boundary->process,
	    MATCHING_TAG (GFS_BOUNDARY (boundary)),
	    boundary->comm,
	    &status);
  MPI_Get_count (&status, MPI_DOUBLE, &count);
  g_assert (count == boundary->rcvcount);

#ifdef PROFILE_MPI
  end = MPI_Wtime ();
  gts_range_add_value (&domain->mpi_wait, end - start);
#endif /* PROFILE_MPI */ 
 
  boundary->rcvcount = 0;
  switch (GFS_BOUNDARY (boundary)->type) {
  case GFS_BOUNDARY_FACE_VARIABLE:
    ftt_face_traverse_boundary (GFS_BOUNDARY (boundary)->root,
				GFS_BOUNDARY (boundary)->d,
				FTT_PRE_ORDER, flags, max_depth,
				(FttFaceTraverseFunc) face_update, boundary);
    break;

  case GFS_BOUNDARY_MATCH_VARIABLE:
    match_update (GFS_BOUNDARY (boundary)->root, boundary);
    ftt_cell_flatten (GFS_BOUNDARY (boundary)->root, 
		      GFS_BOUNDARY (boundary)->d,
		      (FttCellCleanupFunc) gfs_cell_cleanup, NULL);
    break;

  default:
    ftt_cell_traverse (GFS_BOUNDARY (boundary)->root,
		       FTT_PRE_ORDER, flags, max_depth,
		       (FttCellTraverseFunc) center_update, boundary);
  }
  g_assert (boundary->rcvcount == count);
}

static void synchronize (GfsBoundary * bb)
{
  GfsBoundaryMpi * boundary = GFS_BOUNDARY_MPI (bb);
  MPI_Status status;
  guint i;
#ifdef PROFILE_MPI
  GfsDomain * domain = gfs_box_domain (bb->box);
  gdouble start, end;

  start = MPI_Wtime ();
#endif /* PROFILE_MPI */

  /* wait for completion of non-blocking send(s) */
  for (i = 0; i < boundary->nrequest; i++)
    MPI_Wait (&(boundary->request[i]), &status);
#ifdef PROFILE_MPI
  end = MPI_Wtime ();
  gts_range_add_value (&domain->mpi_wait, end - start);
#endif /* PROFILE_MPI */
  boundary->nrequest = 0;
  boundary->sndcount = 0;

  if (bb->type == GFS_BOUNDARY_MATCH_VARIABLE)
    set_buffers_size (boundary);
}

static GtsColor mpi_color (GtsObject * o)
{
  GtsColor c = { 1., 0., 0. }; /* red */

  return c;
}

static void gfs_boundary_mpi_class_init (GfsBoundaryMpiClass * klass)
{
  GfsBoundaryClass * parent_class = GFS_BOUNDARY_CLASS (klass);

  parent_class->match             = match;
  parent_class->send              = send;
  parent_class->receive           = receive;
  parent_class->synchronize       = synchronize;

  GTS_OBJECT_CLASS (klass)->color = mpi_color;
  GTS_OBJECT_CLASS (klass)->destroy = boundary_mpi_destroy;
  GTS_OBJECT_CLASS (klass)->read    = boundary_mpi_read;
}

static void gfs_boundary_mpi_init (GfsBoundaryMpi * boundary)
{
  GfsBc * b = GFS_BOUNDARY (boundary)->default_bc;

  b->bc                = (FttFaceTraverseFunc) center_mpi;
  b->homogeneous_bc    = (FttFaceTraverseFunc) center_mpi;
  b->face_bc           = (FttFaceTraverseFunc) face_mpi;

  boundary->comm = MPI_COMM_WORLD;
  boundary->process = -1; 
  boundary->id = -1;

  boundary->nrequest = 0;
  boundary->sndbuf = g_array_new (FALSE, FALSE, sizeof (gdouble));
  boundary->rcvbuf = g_array_new (FALSE, FALSE, sizeof (gdouble));
  boundary->sndcount = boundary->rcvcount = 0;
}

GfsBoundaryMpiClass * gfs_boundary_mpi_class (void)
{
  static GfsBoundaryMpiClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_boundary_mpi_info = {
      "GfsBoundaryMpi",
      sizeof (GfsBoundaryMpi),
      sizeof (GfsBoundaryMpiClass),
      (GtsObjectClassInitFunc) gfs_boundary_mpi_class_init,
      (GtsObjectInitFunc) gfs_boundary_mpi_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    int * tagub, flag, maxtag;

    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_boundary_class ()),
				  &gfs_boundary_mpi_info);
    MPI_Attr_get (MPI_COMM_WORLD, MPI_TAG_UB, &tagub, &flag);
    if (flag)
      maxtag = *tagub;
    else
      maxtag = 32767; /* minimum value from MPI standard specification */
    tag_shift = maxtag/FTT_NEIGHBORS;
  }

  return klass;
}

GfsBoundaryMpi * gfs_boundary_mpi_new (GfsBoundaryMpiClass * klass,
				       GfsBox * box,
				       FttDirection d,
				       gint process,
				       gint id)
{
  GfsBoundaryMpi * boundary;
  int comm_size;

  MPI_Comm_size (MPI_COMM_WORLD, &comm_size);

  g_return_val_if_fail (process >= 0 && process < comm_size, NULL);

  boundary = GFS_BOUNDARY_MPI (gfs_boundary_new (GFS_BOUNDARY_CLASS (klass),
					       box,
					       d));
  boundary->process = process;
  boundary->id = id;

  if (id >= tag_shift)
    g_warning ("GfsBoundaryMpi id (%d) is larger than the maximum MPI tag value\n"
	       "allowed on this system (%d)", id, tag_shift);

  set_buffers_size (boundary);

  return boundary;
}
