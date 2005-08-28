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
#include <string.h>
#include "solid.h"
#include "vof.h"
#include "variable.h"

/**
 * gfs_cell_fluid:
 * @cell: a #FttCell.
 * 
 * Sets @cell and all its descendants as full fluid cells.
 */
void gfs_cell_fluid (FttCell * cell)
{
  g_return_if_fail (cell != NULL);

  if (GFS_STATE (cell)->solid) {
    g_free (GFS_STATE (cell)->solid);
    GFS_STATE (cell)->solid = NULL;
  }

  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren child;
    guint i;
 
    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i])
	gfs_cell_fluid (child.c[i]);
  }
}

static gdouble segment_triangle_intersection (GtsPoint * E, GtsPoint * D,
					      GtsTriangle * t,
					      gboolean * inside)
{
  GtsVertex * vA, * vB, * vC;
  GtsPoint * A, * B, * C;
  gint ABCE, ABCD, ADCE, ABDE, BCDE;
  GtsEdge * AB, * BC, * CA;
  gdouble a, b;
  gboolean reversed = FALSE;

  gts_triangle_vertices_edges (t, NULL, &vA, &vB, &vC, &AB, &BC, &CA);
  A = GTS_POINT (vA);
  B = GTS_POINT (vB);
  C = GTS_POINT (vC);
  ABCE = gts_point_orientation_3d_sos (A, B, C, E);
  ABCD = gts_point_orientation_3d_sos (A, B, C, D);
  if (ABCE < 0 || ABCD > 0) {
    GtsPoint * tmpp;
    gint tmp;

    tmpp = E; E = D; D = tmpp;
    tmp = ABCE; ABCE = ABCD; ABCD = tmp;
    reversed = TRUE;
  }
  if (ABCE < 0 || ABCD > 0)
    return -1.;
  ADCE = gts_point_orientation_3d_sos (A, D, C, E);
  if (ADCE < 0)
    return -1.;
  ABDE = gts_point_orientation_3d_sos (A, B, D, E);
  if (ABDE < 0)
    return -1.;
  BCDE = gts_point_orientation_3d_sos (B, C, D, E);
  if (BCDE < 0)
    return -1.;
  *inside = reversed ? (ABCD < 0) : (ABCE < 0);
  a = gts_point_orientation_3d (A, B, C, E);
  b = gts_point_orientation_3d (A, B, C, D);
  if (a != b)
    return reversed ? 1. - a/(a - b) : a/(a - b);
  /* D and E are contained within ABC */
  g_assert (a == 0.);
  return 0.5;
}

typedef struct {
  GtsPoint p[4];
  gdouble x[4];
  guint n[4];
  gint inside[4];
} CellFace;

static void face_fractions (CellFace * f, GfsSolidVector * solid, FttVector * h)
{
  static guint etod[] = { 3, 0, 2, 1 };
  guint k, m;
  gboolean ins;
  guint o = 0;
  GtsPoint r[2];
  gdouble a;

  solid->a = 0.;
  solid->cm.x = solid->cm.y = solid->cm.z = 0.;
  solid->ca.z = 0.;
      
  for (m = 0; m < 4 && f->n[m] == 0; m++);
  ins = f->inside[m] < 0;
  for (k = m; k < m + 4; k++) {
    guint i = k % 4, i1 = (i + 1) % 4;
    gdouble x1 = f->p[i].x, y1 = f->p[i].y, x2 = f->p[i1].x, y2 = f->p[i1].y;
    if (f->n[i] > 0) {
      g_assert (ins == (f->inside[i] < 0));
      solid->s[etod[i]] = ins ? f->x[i] : 1. - f->x[i];
      r[o].x = (1. - f->x[i])*x1 + f->x[i]*x2;
      r[o].y = (1. - f->x[i])*y1 + f->x[i]*y2;
      if (ins) {
	x2 = r[o].x; y2 = r[o].y;
      }
      else {
	x1 = r[o].x; y1 = r[o].y;
      }
      solid->a += (x1 + x2)*(y2 - y1);
      solid->cm.x += (x1 - x2)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
      solid->cm.y += (y2 - y1)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
      o++;
      if (o == 2) {
	o = 0;
	if (ins) {
	  x1 = r[1].x; y1 = r[1].y;
	  x2 = r[0].x; y2 = r[0].y;	    
	}
	else {
	  x1 = r[0].x; y1 = r[0].y;
	  x2 = r[1].x; y2 = r[1].y;	    
	}
	solid->a += (x1 + x2)*(y2 - y1);
	solid->cm.x += (x1 - x2)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
	solid->cm.y += (y2 - y1)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
	solid->ca.x = (x1 + x2)/2.;
	solid->ca.y = (y1 + y2)/2.;
      }
      ins = !ins;
    }
    else if (ins) {
      solid->s[etod[i]] = 1.;
      solid->a += (x1 + x2)*(y2 - y1);
      solid->cm.x += (x1 - x2)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
      solid->cm.y += (y2 - y1)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
    }
    else
      solid->s[etod[i]] = 0.;
  }

  a = solid->a/(2.*h->x*h->y);
  if (a > 1e-4) {
    solid->cm.x /= 3.*solid->a;
    solid->cm.y /= 3.*solid->a;
  }
  else {
    guint n = 0;

    solid->cm.x = solid->cm.y = 0.;
    for (m = 0; m < 4 && f->n[m] == 0; m++);
    ins = f->inside[m] < 0;
    for (k = m; k < m + 4; k++) {
      guint i = k % 4, i1 = (i + 1) % 4;
      gdouble x1 = f->p[i].x, y1 = f->p[i].y, x2 = f->p[i1].x, y2 = f->p[i1].y;
      if (f->n[i] > 0) {
	gdouble x = (1. - f->x[i])*x1 + f->x[i]*x2;
	gdouble y = (1. - f->x[i])*y1 + f->x[i]*y2;

	g_assert (ins == (f->inside[i] < 0));
	solid->cm.x += x;
	solid->cm.y += y;
	n++;
	if (ins) {
	  solid->cm.x += x1;
	  solid->cm.y += y1;
	  n++;
	}
	ins = !ins;
      }
      else if (ins) {
	solid->cm.x += x1;
	solid->cm.y += y1;
	n++;
      }
    }
    g_assert (n > 0);
    solid->cm.x /= n;
    solid->cm.y /= n;
  }
  solid->a = a;
}

#if FTT_2D
static void triangle_face_intersection (GtsTriangle * t, CellFace * f)
{
  guint i;

  for (i = 0; i < 4; i++) {
    gboolean ins;
    gdouble x = segment_triangle_intersection (&(f->p[i]), &(f->p[(i + 1) % 4]), t, &ins);

    if (x != -1.) {
      f->x[i] += x; f->n[i]++;
      f->inside[i] += ins ? 1 : -1;
    }
  }
}

static void set_solid_fractions_from_surface (FttCell * cell,
					      GtsSurface * s)
{
  GfsSolidVector * solid = GFS_STATE (cell)->solid;
  FttVector h;
  FttVector p;
  CellFace f;
  guint i, n1 = 0;

  ftt_cell_pos (cell, &p);
  h.x = h.y = ftt_cell_size (cell);
  f.p[0].x = p.x - h.x/2.; f.p[0].y = p.y - h.y/2.; f.p[0].z = 0.;
  f.p[1].x = p.x + h.x/2.; f.p[1].y = p.y - h.y/2.; f.p[1].z = 0.;
  f.p[2].x = p.x + h.x/2.; f.p[2].y = p.y + h.y/2.; f.p[2].z = 0.;
  f.p[3].x = p.x - h.x/2.; f.p[3].y = p.y + h.y/2.; f.p[3].z = 0.;
  f.x[0] = f.x[1] = f.x[2] = f.x[3] = 0.;
  f.n[0] = f.n[1] = f.n[2] = f.n[3] = 0;
  f.inside[0] = f.inside[1] = f.inside[2] = f.inside[3] = 0;

  gts_surface_foreach_face (s, (GtsFunc) triangle_face_intersection, &f);

  for (i = 0; i < 4; i++)
    if (f.n[i] % 2 != 0) {
      f.x[i] /= f.n[i];
      n1++;
    }
    else
      f.n[i] = 0;

  switch (n1) {
  case 0:
    if (solid) {
      g_free (solid);
      GFS_STATE (cell)->solid = NULL;
    }
    break;
  case 2: case 4: {
    if (!solid)
      GFS_STATE (cell)->solid = solid = g_malloc (sizeof (GfsSolidVector));
    face_fractions (&f, solid, &h);
    break;
  }
  default:
    g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	   "the surface is probably not closed (n1 = %d)", n1);
  }
}
#else /* 2D3 or 3D */
#include "isocube.h"

typedef struct {
  GtsPoint p[8];
  gdouble x[12];
  guint n[12];
  gint inside[12];
} CellCube;

static void triangle_cube_intersection (GtsTriangle * t, CellCube * cube)
{
  guint i;

  for (i = 0; i < 12; i++) {
    gboolean ins;
    gdouble x = segment_triangle_intersection (&cube->p[edge1[i][0]], &cube->p[edge1[i][1]], 
					       t, &ins);
    if (x != -1.) {
      cube->x[i] += x; cube->n[i]++;
      cube->inside[i] += ins ? 1 : -1;
    }
  }
}

static void rotate (CellFace * f, FttVector * h, FttComponent c)
{
  guint i;

  switch (c) {
  case FTT_X: 
    for (i = 0; i < 4; i++) {
      f->p[i].x = f->p[i].y; f->p[i].y = f->p[i].z;
    }
    h->x = h->y; h->y = h->z;
    break;
  case FTT_Y:
    for (i = 0; i < 4; i++)
      f->p[i].y = f->p[i].z;
    h->y = h->z;
    break;
  case FTT_Z:
    break;
  default:
    g_assert_not_reached ();
  }
}

static void cell_size (FttCell * cell, FttVector * h)
{
  h->x = h->y = ftt_cell_size (cell);
#if FTT_2D3
  h->z = 1.;
#else  /* 3D */
  h->z = h->x;
#endif /* 3D */
}

static void set_solid_fractions_from_surface (FttCell * cell, GtsSurface * s)
{
  GfsSolidVector * solid = GFS_STATE (cell)->solid;
  CellCube cube;
  FttVector o, ca = {0., 0., 0.}, h;
  guint i, n1 = 0;
  gint inside[8] = {0,0,0,0,0,0,0,0};

  ftt_cell_pos (cell, &o);
  cell_size (cell, &h);
  for (i = 0; i < FTT_DIMENSION; i++)
    (&o.x)[i] -= (&h.x)[i]/2.;
  for (i = 0; i < 12; i++) {
    cube.x[i] = 0.; cube.n[i] = 0; cube.inside[i] = 0;
  }
  for (i = 0; i < 8; i++) { /* for each vertex of the cube */
    cube.p[i].x = o.x + h.x*vertex[i].x;
    cube.p[i].y = o.y + h.y*vertex[i].y;
    cube.p[i].z = o.z + h.z*vertex[i].z;
  }

  gts_surface_foreach_face (s, (GtsFunc) triangle_cube_intersection, &cube);

  for (i = 0; i < 12; i++) /* for each edge of the cube */
    if (cube.n[i] % 2 != 0) { /* only for odd number of intersections */
      guint j = edge1[i][0], k = edge1[i][1];

      /* intersection vertex position is the average of all the n[i] intersections */
      cube.x[i] /= cube.n[i];

      /* average of all intersections */
      ca.x += (1. - cube.x[i])*cube.p[j].x + cube.x[i]*cube.p[k].x;
      ca.y += (1. - cube.x[i])*cube.p[j].y + cube.x[i]*cube.p[k].y;
      ca.z += (1. - cube.x[i])*cube.p[j].z + cube.x[i]*cube.p[k].z;

      g_assert (inside[j] == 0 || inside[j] == cube.inside[i]);
      g_assert (inside[k] == 0 || inside[k] == - cube.inside[i]);
      inside[j] = cube.inside[i];
      inside[k] = - cube.inside[i];
      n1++;
    }
    else
      cube.n[i] = 0;

  if (n1 == 0) { /* no intersections */
    if (solid) {
      g_free (solid);
      GFS_STATE (cell)->solid = NULL;
    }
    return;
  }

  if (!solid)
    GFS_STATE (cell)->solid = solid = g_malloc0 (sizeof (GfsSolidVector));

  /* compute face fractions */
  for (i = 0; i < FTT_NEIGHBORS; i++) {
    CellFace f;
    guint j, n2;

    n2 = 0;
    for (j = 0; j < 4; j++) { /* initialise face i */
      guint e = face[i][j][0];

      f.p[j] = cube.p[face_v[i][j]];
      f.n[j] = cube.n[e];
      if (f.n[j]) n2++;
      if (face[i][j][1]) {
	f.x[j] = 1. - cube.x[e];
	f.inside[j] = - cube.inside[e];
      }
      else {
	f.x[j] = cube.x[e];
	f.inside[j] = cube.inside[e];
      }
    }

    switch (n2) {
    case 0: { /* the face is not cut */
      gint ins = 0;

      /* checks whether the face vertices are inside or outside */
      for (j = 0; j < 4; j++) {
	gint k = inside[face_v[i][j]];
	if (k) {
	  g_assert (ins == 0 || ins == k);
	  ins = k;
	}
      }
      g_assert (ins != 0);
      solid->s[i] = ins > 0 ? 0. : 1.;
      break;
    }
    case 2: case 4: { /* the face is cut 2 or 4 times */
      GfsSolidVector sol;
      FttVector h1;

      h1 = h;
      rotate (&f, &h1, i/2);
      face_fractions (&f, &sol, &h1);
      solid->s[i] = sol.a;
      break;
    }
    default:
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	     "the surface is probably not closed (n2 = %d)", n2);
    }
  }

  /* now compute cell fraction, center of area, center of mass */

  /* fixme: the calculation of ca below is not strictly equivalent to
     the true center of area of the piece of surface contained in the
     cell if there are more than 4 intersection points (n1 > 4) */
  ca.x /= n1; ca.y /= n1; ca.z /= n1; 
  solid->ca = ca;
  {
    FttVector m;
    gdouble alpha, n = 0.;
    gboolean sym[FTT_DIMENSION];
    FttComponent c;

    for (c = 0; c < FTT_DIMENSION; c++) {
      (&ca.x)[c] = ((&ca.x)[c] - (&o.x)[c])/(&h.x)[c];
      (&m.x)[c] = solid->s[2*c + 1] - solid->s[2*c];
      if ((&m.x)[c] < 0.) {
	(&m.x)[c] = - (&m.x)[c];
	(&ca.x)[c] = 1. - (&ca.x)[c];
	sym[c] = TRUE;
      }
      else
	sym[c] = FALSE;
      n += (&m.x)[c];
    }
    if (n > 0.) {
      m.x /= n; m.y /= n; m.z /= n;
      alpha = m.x*ca.x + m.y*ca.y + m.z*ca.z;
      solid->a = gfs_plane_volume (&m, alpha, 1.);
      gfs_plane_center (&m, alpha, solid->a, &solid->cm);
      for (c = 0; c < FTT_DIMENSION; c++)
	(&solid->cm.x)[c] = (&o.x)[c] + 
	  (sym[c] ? 1. - (&solid->cm.x)[c] : (&solid->cm.x)[c])*(&h.x)[c];
    }
    else { /* degenerate intersections */
      solid->a = 0.;
      for (i = 0; i < FTT_NEIGHBORS; i++)
	solid->a += solid->s[i];
      solid->a /= FTT_NEIGHBORS;
      if (solid->a == 0. || solid->a == 1.) {
	g_free (solid);
	GFS_STATE (cell)->solid = NULL;
      }
    }
  }
}
#endif /* 2D3 or 3D */

static gdouble solid_sa (GfsSolidVector * s)
{
  gdouble sa2 = 0.;
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++) {
    gdouble n = s->s[2*c] - s->s[2*c + 1];

    sa2 += n*n;
  }
  return sqrt (sa2);
}

/**
 * gfs_cell_init_solid_fractions_from_children:
 * @cell: a #FttCell.
 *
 * Uses the values of the solid fractions of the children of @cell to
 * compute the values of its solid fractions.
 *
 * This function fails if @cell is a leaf of the cell tree.  
 */
void gfs_cell_init_solid_fractions_from_children (FttCell * cell)
{
  FttCellChildren child;
  guint i, j;
  gdouble w = 0., wa = 0.;
  gboolean cell_is_solid = TRUE;
  gboolean cell_is_mixed = FALSE;
  FttVector cm = { 0., 0., 0.};
  FttVector ca = { 0., 0., 0.};

  g_return_if_fail (cell != NULL);
  g_return_if_fail (!FTT_CELL_IS_LEAF (cell));

  ftt_cell_children (cell, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
      if (GFS_IS_FLUID (child.c[i])) {
	FttVector p;

	w += 1.;
	ftt_cell_pos (child.c[i], &p);
	cm.x += p.x; cm.y += p.y; cm.z += p.z;
	cell_is_solid = FALSE;
      }
      else {
	GfsSolidVector * solid = GFS_STATE (child.c[i])->solid;
	gdouble sa = solid_sa (solid) + 1e-9;

	w += solid->a; wa += sa;
	cm.x += solid->cm.x*solid->a;
	cm.y += solid->cm.y*solid->a;
	cm.z += solid->cm.z*solid->a;
	ca.x += solid->ca.x*sa;
	ca.y += solid->ca.y*sa;
	ca.z += solid->ca.z*sa;
	cell_is_mixed = TRUE;
      }
    }

  if (cell_is_mixed) {
    GfsSolidVector * solid = GFS_STATE (cell)->solid;

    if (solid == NULL)
      GFS_STATE (cell)->solid = solid = g_malloc0 (sizeof (GfsSolidVector));

    solid->a = w/FTT_CELLS; 
    solid->cm.x = cm.x/w;
    solid->cm.y = cm.y/w;
    solid->cm.z = cm.z/w;
    solid->ca.x = ca.x/wa;
    solid->ca.y = ca.y/wa;
    solid->ca.z = ca.z/wa;
    for (i = 0; i < FTT_NEIGHBORS; i++) {
      guint n = ftt_cell_children_direction (cell, i, &child);

      w = 0.;
      for (j = 0; j < n; j++) {
	FttCell * c = child.c[j];

	if (c)
	  w += GFS_IS_FLUID (c) ? 1. : GFS_STATE (c)->solid->s[i];
      }
      solid->s[i] = w/n;
    }
  }
  else { /* !cell_is_mixed */
    if (GFS_STATE (cell)->solid) {
      g_free (GFS_STATE (cell)->solid);
      GFS_STATE (cell)->solid = NULL;
    }
    g_assert (!cell_is_solid);
  }
}

static void push_leaf (GtsFifo * fifo, FttCell * cell, FttDirection d, gdouble a,
		       GfsVariable * status)
{
  if (FTT_CELL_IS_LEAF (cell)) {
    if (!GFS_IS_MIXED (cell) && GFS_VARIABLE (cell, status->i) == 0.) {
      GFS_VARIABLE (cell, status->i) = a;
      gts_fifo_push (fifo, cell);
    }
  }
  else {
    FttCellChildren child;
    guint i, n;
    
    n = ftt_cell_children_direction (cell, FTT_OPPOSITE_DIRECTION (d), &child);
    for (i = 0; i < n; i++)
      if (child.c[i] && !GFS_IS_MIXED (child.c[i]) && GFS_VARIABLE (child.c[i], status->i) == 0.) {
	g_assert (FTT_CELL_IS_LEAF (child.c[i]));
	GFS_VARIABLE (child.c[i], status->i) = a;
	gts_fifo_push (fifo, child.c[i]);
      }
  }  
}

static void paint_leaf (GtsFifo * fifo, gdouble a, GfsVariable * status)
{
  FttCell * cell;

  while ((cell = gts_fifo_pop (fifo))) {
    FttDirection i;
    FttCellNeighbors n;
    
    ftt_cell_neighbors (cell, &n);
    for (i = 0; i < FTT_NEIGHBORS; i++)
      if (n.c[i] && !GFS_CELL_IS_BOUNDARY (n.c[i]))
	push_leaf (fifo, n.c[i], i, a, status);
  }
}

static void paint_mixed_leaf (FttCell * cell, GfsVariable * status)
{
  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * solid = GFS_STATE (cell)->solid;
    GtsFifo * fifo;
    FttCell * n;
    FttDirection i;

    fifo = gts_fifo_new ();
    for (i = 0; i < FTT_NEIGHBORS; i++)
      if ((n = ftt_cell_neighbor (cell, i)) && !GFS_CELL_IS_BOUNDARY (n)) {
	if (solid->s[i] == 0. || solid->s[i] == 1.) {
	  push_leaf (fifo, n, i, solid->s[i] + 1., status);
	  paint_leaf (fifo, solid->s[i] + 1., status);
	}
	else if (!FTT_CELL_IS_LEAF (n)) {
	  FttCellChildren child;
	  guint j, k;
	  gdouble w = 0.;

	  k = ftt_cell_children_direction (n, FTT_OPPOSITE_DIRECTION (i), &child);
	  for (j = 0; j < k; j++)
	    if (child.c[j])
	      w += GFS_IS_FLUID (child.c[j]) ? 1. : 
		GFS_STATE (child.c[j])->solid->s[FTT_OPPOSITE_DIRECTION (i)];
	  solid->s[i] = w/k;
	  g_assert (solid->s[i] > 0. && solid->s[i] < 1.);
	}
      }
    gts_fifo_destroy (fifo);
  }
}

typedef struct {
  gboolean destroy_solid;
  FttCellCleanupFunc cleanup;
  gpointer data;
  GfsVariable * status;
} InitSolidParams;

static void solid_fractions_from_children (FttCell * cell, InitSolidParams * p)
{
  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren child;
    guint i;
    
    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i])
	solid_fractions_from_children (child.c[i], p);
    if (FTT_CELL_IS_LEAF (cell))
      /* all the children have been destroyed i.e. the cell is solid */
      GFS_VARIABLE (cell, p->status->i) = 1.;
    else {
      gfs_cell_init_solid_fractions_from_children (cell);
      if (p->destroy_solid)
	GFS_VARIABLE (cell, p->status->i) = 0.;
      else if (!GFS_IS_MIXED (cell)) {
	ftt_cell_children (cell, &child);
	GFS_VARIABLE (cell, p->status->i) = 1.;
	for (i = 0; i < FTT_CELLS; i++)
	  if (child.c[i] && GFS_VARIABLE (child.c[i], p->status->i) == 2.)
	    GFS_VARIABLE (cell, p->status->i) = 2.;
      }
    }
  }
  if (p->destroy_solid && GFS_VARIABLE (cell, p->status->i) == 1.) {
    if (FTT_CELL_IS_ROOT (cell))
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	     "root cell is entirely outside of the fluid domain\n"
	     "the solid surface orientation may be incorrect");
    else
      ftt_cell_destroy (cell, p->cleanup, p->data);
  }
}

static void foreach_box (GfsBox * box, InitSolidParams * p)
{
  solid_fractions_from_children (box->root, p);
}

/**
 * gfs_domain_init_solid_fractions:
 * @domain: a #GfsDomain.
 * @s: an orientable surface defining the solid boundary.
 * @destroy_solid: controls what to do with solid cells.
 * @cleanup: a #FttCellCleanupFunc or %NULL.
 * @data: user data to pass to @cleanup.
 * @status: a temporary variable or %NULL.
 *
 * Initializes the solid fractions of all the cells of @domain.
 *
 * If @destroy_solid is set to %TRUE, the cells entirely contained in
 * the solid are destroyed using @cleanup as cleanup function.  
 */
void gfs_domain_init_solid_fractions (GfsDomain * domain,
				      GtsSurface * s,
				      gboolean destroy_solid,
				      FttCellCleanupFunc cleanup,
				      gpointer data,
				      GfsVariable * status)
{
  InitSolidParams p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (s != NULL);

  p.destroy_solid = destroy_solid;
  p.cleanup = cleanup;
  p.data = data;
  p.status = status ? status : gfs_temporary_variable (domain);
  gfs_domain_traverse_cut (domain, s, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			   (FttCellTraverseCutFunc) set_solid_fractions_from_surface, NULL);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, p.status);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) paint_mixed_leaf, p.status);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) foreach_box, &p);
  if (status == NULL)
    gts_object_destroy (GTS_OBJECT (p.status));
}

static gboolean check_area_fractions (const FttCell * root)
{
  guint i, level;
  FttCellNeighbors neighbor;
  gboolean ret = TRUE;
  GfsSolidVector * solid;

  level = ftt_cell_level (root);
  ftt_cell_neighbors (root, &neighbor);
  solid = GFS_STATE (root)->solid;

  if (solid) {
    GtsBBox bb;

    ftt_cell_bbox (root, &bb);
    if (!gts_bbox_point_is_inside (&bb, &solid->cm)) {
      g_warning ("file %s: line %d (%s): cm is not inside cell",
		 __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION);
      ret = FALSE;
      g_assert_not_reached ();
    }
    if (!gts_bbox_point_is_inside (&bb, &solid->ca)) {
      g_warning ("file %s: line %d (%s): ca is not inside cell",
		 __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION);
      ret = FALSE;
      g_assert_not_reached ();
    }
  }

  for (i = 0; i < FTT_NEIGHBORS; i++)
    if (neighbor.c[i]) {
      GfsSolidVector * nsolid = GFS_STATE (neighbor.c[i])->solid;

      if (ftt_cell_level (neighbor.c[i]) == level) {
	if (GFS_IS_FLUID (root)) {
	  if (!GFS_IS_FLUID (neighbor.c[i]) && 
	      1. - nsolid->s[FTT_OPPOSITE_DIRECTION (i)] >= 1e-10) {
	    FttVector p;
	    ftt_cell_pos (root, &p);
	    g_warning ("file %s: line %d (%s): (%g,%g,%g): s[%d]: %g",		       
		       __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		       p.x, p.y, p.z,
		       FTT_OPPOSITE_DIRECTION (i),
		       nsolid->s[FTT_OPPOSITE_DIRECTION (i)]);
	    ret = FALSE;
	    nsolid->s[FTT_OPPOSITE_DIRECTION (i)] = 1.;
	  }
	}
	else if (GFS_IS_MIXED (neighbor.c[i])) {
	  if (fabs (solid->s[i] - 
		    nsolid->s[FTT_OPPOSITE_DIRECTION (i)]) >= 1e-10) {
	    FttVector p;
	    ftt_cell_pos (root, &p);
	    g_warning ("file %s: line %d (%s): (%g,%g,%g): s[%d]: %g neighbor->s[%d]: %g",
		       __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		       p.x, p.y, p.z,
		       i, solid->s[i],		       
		       FTT_OPPOSITE_DIRECTION (i),
		       nsolid->s[FTT_OPPOSITE_DIRECTION (i)]);
	    ret = FALSE;
	    nsolid->s[FTT_OPPOSITE_DIRECTION (i)] = solid->s[i];
	  }
	}
	else if (1. - solid->s[i] >= 1e-10) {
	  FttVector p;
	  ftt_cell_pos (root, &p);
	  g_warning ("file %s: line %d (%s): (%g,%g,%g): s[%d]: %g",		     
		     __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		     p.x, p.y, p.z,
		     i, solid->s[i]);
	  ret = FALSE;
	  solid->s[i] = 1.;
	}
      }
      else { /* fine/coarse boundary */
	g_assert (ftt_cell_level (neighbor.c[i]) == level - 1);
	if (GFS_IS_FLUID (neighbor.c[i]) && GFS_IS_MIXED (root) && 1. - solid->s[i] >= 1e-10) {
	  FttVector p;
	  ftt_cell_pos (root, &p);
	  g_warning ("file %s: line %d (%s): (%g,%g,%g): s[%d]: %g",
		     __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		     p.x, p.y, p.z,
		     i, solid->s[i]);
	  ret = FALSE;
	  solid->s[i] = 1.;
	}
      }
    }
  
  if (!FTT_CELL_IS_LEAF (root)) {
    FttCellChildren child;

    ftt_cell_children (root, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i] && !check_area_fractions (child.c[i]))
	ret = FALSE;
  }

  return ret;
}

static void check_solid_fractions (FttCell * cell, gboolean * ret)
{
  FttCellChildren children;
  guint n;

  ftt_cell_children (cell, &children);
  if (!GFS_IS_MIXED (cell)) {
    for (n = 0; n < FTT_CELLS; n++)
      if (children.c[n] && GFS_IS_MIXED (children.c[n])) {
	g_warning ("file %s: line %d (%s): children[%d] is mixed (%g)"
		   " parent is not",
                   __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		   n, GFS_STATE (children.c[n])->solid->a);
	*ret = FALSE;
      }
  }
  else {
    gdouble a = 0.;

    for (n = 0; n < FTT_CELLS; n++)
      if (children.c[n]) {
	if (GFS_IS_MIXED (children.c[n]))
	  a += GFS_STATE (children.c[n])->solid->a;
	else
	  a += 1.;
      }
    a /= FTT_CELLS;
    if (fabs (GFS_STATE (cell)->solid->a - a) >= 1e-10) {
      g_warning ("file %s: line %d (%s): children->a: %g parent->a: %g",
		 __FILE__, __LINE__, G_GNUC_PRETTY_FUNCTION,
		 a, GFS_STATE (cell)->solid->a);
	*ret = FALSE;
    }
  }
}

/**
 * gfs_cell_check_solid_fractions:
 * @root: the root #FttCell of the cell tree to check.
 * 
 * Checks the consistency of the solid fractions of each cell of the
 * cell tree relative to the neighboring solid fractions.
 *
 * Returns: %TRUE if the solid fractions are consistent, %FALSE otherwise.
 */
gboolean gfs_cell_check_solid_fractions (FttCell * root)
{
  gboolean ret = TRUE;

  g_return_val_if_fail (root != NULL, FALSE);

  ftt_cell_traverse (root, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
		     (FttCellTraverseFunc) check_solid_fractions, &ret);
  return ret & check_area_fractions (root);
}

static void save_solid (FttCell * cell, GfsVariable * c)
{
  GFS_DOUBLE_TO_POINTER (GFS_VARIABLE (cell, c->i)) = GFS_STATE (cell)->solid;
  GFS_STATE (cell)->solid = NULL;
}

static void restore_solid (FttCell * cell, gpointer * data)
{
  GfsVariable * c = data[0];
  gboolean * not_cut = data[1];
  GfsVariable * status = data[2];
  GfsSolidVector * solid = GFS_STATE (cell)->solid;

  GFS_STATE (cell)->solid = GFS_DOUBLE_TO_POINTER (GFS_VARIABLE (cell, c->i));
  if (solid) {
    GFS_VARIABLE (cell, c->i) = solid->a;
    g_free (solid);
    *not_cut = FALSE;
  }
  else if (GFS_VARIABLE (cell, status->i) == 0.) {
    g_assert (*not_cut);
    GFS_VARIABLE (cell, c->i) = 0.;
  }
  else {
    g_assert (GFS_VARIABLE (cell, status->i) == 1. || GFS_VARIABLE (cell, status->i) == 2.);
    GFS_VARIABLE (cell, c->i) = GFS_VARIABLE (cell, status->i) - 1.;
  }
}

/**
 * gfs_domain_init_fraction:
 * @domain: a #GfsDomain.
 * @s: an orientable surface defining the interface boundary.
 * @c: a #GfsVariable.
 *
 * Initializes the fraction @c of the interface @s contained in all
 * the cells of @domain.
 */
void gfs_domain_init_fraction (GfsDomain * domain,
			       GtsSurface * s,
			       GfsVariable * c)
{
  gboolean not_cut = TRUE;
  gpointer data[3];
  GfsVariable * status;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (c != NULL);

  status = gfs_temporary_variable (domain);

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) save_solid, c);
  gfs_domain_init_solid_fractions (domain, s, FALSE, NULL, NULL, status);
  data[0] = c;
  data[1] = &not_cut;
  data[2] = status;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) restore_solid, data);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, c);

  gts_object_destroy (GTS_OBJECT (status));
}

/**
 * gfs_cell_cm:
 * @cell: a #FttCell.
 * @cm: a #FttVector.
 *
 * Fills @cm with the coordinates of the center of mass of @cell.
 */
void gfs_cell_cm (const FttCell * cell, FttVector * cm)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (cm != NULL);

  if (GFS_IS_MIXED (cell))
    *cm = GFS_STATE (cell)->solid->cm;
  else
    ftt_cell_pos (cell, cm);
}

/**
 * gfs_solid_normal:
 * @cell: a #FttCell.
 * @n: a #FttVector.
 *
 * Fills @n with the components of the average unit normal to the
 * fraction of solid boundary contained in @cell, multiplied by the
 * area of the fraction of solid boundary contained in @cell.
 */
void gfs_solid_normal (const FttCell * cell, FttVector * n)
{
  GfsSolidVector * s;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (n != NULL);

  if ((s = GFS_STATE (cell)->solid)) {
    gdouble size = ftt_cell_size (cell);
    FttComponent c;

#if (!FTT_2D)
    size *= size;
#endif /* 3D */

    for (c = 0; c < FTT_DIMENSION; c++)
      (&n->x)[c] = (s->s[2*c + 1] - s->s[2*c])*size;
  }
  else
    n->x = n->y = n->z = 0.;
}

/**
 * gfs_face_ca:
 * @face: a #FttCellFace.
 * @ca: a #FttVector.
 *
 * Fills @ca with the coordinates of the center of area of @face.
 */
void gfs_face_ca (const FttCellFace * face, FttVector * ca)
{
  gdouble f;

  g_return_if_fail (face != NULL);
  g_return_if_fail (ca != NULL);

  ftt_face_pos (face, ca);
  if ((f = GFS_FACE_FRACTION (face)) < 1.) {
    GfsSolidVector * s = GFS_STATE (face->cell)->solid;
    gdouble h = ftt_cell_size (face->cell);
#if FTT_2D
    FttComponent cp = FTT_ORTHOGONAL_COMPONENT (face->d/2);

    (&ca->x)[cp] += (s->s[2*cp] > s->s[2*cp + 1]) ? (1. - f)/2.*h : (f - 1.)/2.*h;
#else /* 3D */
    static guint perpendicular[FTT_DIMENSION][2] = {
      {FTT_Y, FTT_Z}, {FTT_Z, FTT_X}, {FTT_X, FTT_Y}
    };
    FttComponent c0 = face->d/2;
    FttComponent c1 = perpendicular[c0][0];
    FttComponent c2 = perpendicular[c0][1];
    gboolean s1, s2;
    FttVector m, p;
    gdouble n, alpha;

    m.x = s->s[2*c1 + 1] - s->s[2*c1];
    m.y = s->s[2*c2 + 1] - s->s[2*c2];
    s1 = (m.x < 0.);
    s2 = (m.y < 0.);
    m.x = fabs (m.x);
    m.y = fabs (m.y);
    n = m.x + m.y;
    if (n > 0.) {
      m.x /= n;
      m.y /= n;
      alpha = gfs_line_alpha (&m, f);
      gfs_line_center (&m, alpha, f, &p);
      if (s1) p.x = 1. - p.x;
      if (s2) p.y = 1. - p.y;
      (&ca->x)[c1] += (p.x - 0.5)*h;
      (&ca->x)[c2] += (p.y - 0.5)*h;
    }
#endif /* 3D */
  }
}
