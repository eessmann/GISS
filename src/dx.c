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

#include <string.h>
#include <math.h>
#include <dx/dx.h>

#include "simulation.h"
#include "init.h"
#include "solid.h"

static void ftt_cell_neighbors_not_periodic (FttCell * cell,
					     FttCellNeighbors * neighbor)
{
  FttDirection d;
  gdouble size = 2.*ftt_cell_size (cell);
  FttVector p;

  ftt_cell_neighbors (cell, neighbor);
  ftt_cell_pos (cell, &p);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (neighbor->c[d]) {
      FttVector p1;

      ftt_cell_pos (neighbor->c[d], &p1);
      if (fabs ((&p.x)[d/2] - (&p1.x)[d/2]) >= size)
	neighbor->c[d] = NULL;
    }
}

static FttCell * ftt_cell_neighbor_not_periodic (FttCell * cell,
						 FttDirection d)
{
  FttCell * neighbor = ftt_cell_neighbor (cell, d);

  if (neighbor) {
    gdouble size = 2.*ftt_cell_size (cell);
    FttVector p, p1;

    ftt_cell_pos (cell, &p);
    ftt_cell_pos (neighbor, &p1);
    if (fabs ((&p.x)[d/2] - (&p1.x)[d/2]) >= size)
      neighbor = NULL;
  }
  return neighbor;
}

#if FTT_2D

typedef struct {
  FttCell * c1, * c2;
} Edge;

static guint cell_edges (FttCell * cell, 
			 Edge * edge)
{
  guint len = 0;
  FttCellNeighbors neighbor;
  FttCell * n1, * n2, * d1, * d2;
  guint level = ftt_cell_level (cell);

#if FTT_2D3
  g_assert_not_implemented ();
#endif

  ftt_cell_neighbors_not_periodic (cell, &neighbor);

  d1 = NULL;
  if ((n1 = neighbor.c[2])) {
    if (FTT_CELL_IS_LEAF (n1)) {
      if (ftt_cell_level (n1) == level ||
	  FTT_CELL_ID (cell) == 1) {
	d1 = ftt_cell_neighbor_not_periodic (n1, 0);
	if (d1 && !FTT_CELL_IS_LEAF (d1)) {
	  FttCellChildren child;
	  
	  ftt_cell_children_direction (d1, 1, &child);
	  d1 = child.c[1];
	}
      }
    }
    else {
      FttCellChildren child;

      ftt_cell_children_direction (n1, 3, &child);
      edge[len].c1 = child.c[0];
      edge[len++].c2 = child.c[1];
      n1 = child.c[1];
      d1 = ftt_cell_neighbor_not_periodic (child.c[1], 0);
    }
    if (d1) {
      edge[len].c1 = n1;
      edge[len++].c2 = d1;
    }
  }
  d2 = NULL;
  if ((n2 = neighbor.c[0])) {
    if (FTT_CELL_IS_LEAF (n2)) {
      if (ftt_cell_level (n2) == level ||
	  FTT_CELL_ID (cell) == 1) {
	d2 = ftt_cell_neighbor_not_periodic (n2, 2);
	if (d2 && !FTT_CELL_IS_LEAF (d2)) {
	  FttCellChildren child;
	  
	  ftt_cell_children_direction (d2, 3, &child);
	  d2 = child.c[0];
	}
      }
      if (d2) {
	g_assert (!d1 || d1 == d2);
	edge[len].c1 = d2;
	edge[len++].c2 = n2;
      }
    }
    else {
      FttCellChildren child;

      ftt_cell_children_direction (n2, 1, &child);
      d2 = ftt_cell_neighbor_not_periodic (child.c[0], 2);
      if (d2) {
	g_assert (!d1 || d1 == d2);
	edge[len].c1 = d2;
	edge[len++].c2 = child.c[0];
      }
      edge[len].c1 = child.c[0];
      edge[len++].c2 = child.c[1];
    }
  }
#if 0
  d1 = NULL;
  if (neighbor.c[3]) {
    if (FTT_CELL_IS_LEAF (neighbor.c[3])) {
      loop[len++] = neighbor.c[3];
      if (ftt_cell_level (neighbor.c[3]) == level ||
	  FTT_CELL_ID (cell) == 2) {
	d1 = ftt_cell_neighbor_not_periodic (neighbor.c[3], 1);
	if (d1 && !FTT_CELL_IS_LEAF (d1)) {
	  FttCellChildren child;
	  
	  ftt_cell_children_direction (d1, 0, &child);
	  d1 = child.c[0];
	}
      }
    }
    else {
      FttCellChildren child;

      ftt_cell_children_direction (neighbor.c[3], 2, &child);
      loop[len++] = child.c[1];
      loop[len++] = child.c[0];
      d1 = ftt_cell_neighbor_not_periodic (child.c[0], 1);
    }
    if (d1)
      loop[len++] = d1;
  }
  d2 = NULL;
  if (neighbor.c[1]) {
    if (FTT_CELL_IS_LEAF (neighbor.c[1])) {
      if (ftt_cell_level (neighbor.c[1]) == level ||
	  FTT_CELL_ID (cell) == 2) {
	d2 = ftt_cell_neighbor_not_periodic (neighbor.c[1], 3);
	if (d2 && !FTT_CELL_IS_LEAF (d2)) {
	  FttCellChildren child;
	  
	  ftt_cell_children_direction (d2, 2, &child);
	  d2 = child.c[1];
	}
      }
      if (!d1 && d2)
	loop[len++] = d2;
      loop[len++] = neighbor.c[1];
    }
    else {
      FttCellChildren child;

      ftt_cell_children_direction (neighbor.c[1], 0, &child);
      d2 = ftt_cell_neighbor_not_periodic (child.c[1], 3);
      if (!d1 && d2)
	loop[len++] = d2;
      loop[len++] = child.c[1];
      loop[len++] = child.c[0];
    }
  }
  g_assert (d1 == NULL || d2 == NULL || d1 == d2);  
#endif
  return len;
}

#else /* not FTT_2D */

typedef struct {
  FttCell * c1, * c2, * c3;
} Face;

#define CELL_FACES_ASSERT(x) if (!(x)) { \
  FILE * fp = fopen ("/tmp/assert", "wt"); \
  \
  fprintf (fp, "(geometry \"cell\" = {\n"); \
  ftt_cell_draw (cell, fp); \
  fprintf (fp, "})\n"); \
  if (n1) { \
    fprintf (fp, "(geometry \"n1\" = {\n"); \
    ftt_cell_draw (n1, fp); \
    fprintf (fp, "})\n"); \
  } \
  if (n2) { \
    fprintf (fp, "(geometry \"n2\" = {\n"); \
    ftt_cell_draw (n2, fp); \
    fprintf (fp, "})\n"); \
  } \
  if (n3) { \
    fprintf (fp, "(geometry \"n3\" = {\n"); \
    ftt_cell_draw (n3, fp); \
    fprintf (fp, "})\n"); \
  } \
  if (d1) { \
    fprintf (fp, "(geometry \"d1\" = {\n"); \
    ftt_cell_draw (d1, fp); \
    fprintf (fp, "})\n"); \
  } \
  if (d2) { \
    fprintf (fp, "(geometry \"d2\" = {\n"); \
    ftt_cell_draw (d2, fp); \
    fprintf (fp, "})\n"); \
  } \
  if (d3) { \
    fprintf (fp, "(geometry \"d3\" = {\n"); \
    ftt_cell_draw (d3, fp); \
    fprintf (fp, "})\n"); \
  } \
  if (d4) { \
    fprintf (fp, "(geometry \"d4\" = {\n"); \
    ftt_cell_draw (d4, fp); \
    fprintf (fp, "})\n"); \
  } \
  if (d5) { \
    fprintf (fp, "(geometry \"d5\" = {\n"); \
    ftt_cell_draw (d5, fp); \
    fprintf (fp, "})\n"); \
  } \
  if (d6) { \
    fprintf (fp, "(geometry \"d6\" = {\n"); \
    ftt_cell_draw (d6, fp); \
    fprintf (fp, "})\n"); \
  } \
  if (d7) { \
    fprintf (fp, "(geometry \"d7\" = {\n"); \
    ftt_cell_draw (d7, fp); \
    fprintf (fp, "})\n"); \
  } \
  if (d8) { \
    fprintf (fp, "(geometry \"d8\" = {\n"); \
    ftt_cell_draw (d8, fp); \
    fprintf (fp, "})\n"); \
  } \
  if (d9) { \
    fprintf (fp, "(geometry \"d9\" = {\n"); \
    ftt_cell_draw (d9, fp); \
    fprintf (fp, "})\n"); \
  } \
  fclose (fp); \
  g_assert_not_reached (); \
}

static guint cell_faces (FttCell * cell, 
			 Face * face)
{
  guint len = 0;
  FttCellNeighbors neighbor;
  FttCell * n1 = NULL, * n2 = NULL, * n3 = NULL;
  FttCell * d1 = NULL, * d3 = NULL, * d5 = NULL;
  FttCell * d2 = NULL, * d6 = NULL, * d7 = NULL;
  FttCell * d4 = NULL, * d8 = NULL, * d9 = NULL;
  guint level = ftt_cell_level (cell);

  ftt_cell_neighbors_not_periodic (cell, &neighbor);

  if ((n1 = neighbor.c[2])) {
    if (FTT_CELL_IS_LEAF (n1)) {
      if (ftt_cell_level (n1) == level ||
	  FTT_CELL_ID (cell) == 1 ||
	  FTT_CELL_ID (cell) == 5) {
	d1 = ftt_cell_neighbor_not_periodic (n1, 0);
	if (d1) {
	  if (!FTT_CELL_IS_LEAF (d1)) {
	    FttCellChildren child;
	    
	    ftt_cell_children_direction (d1, 1, &child);
	    if (ftt_cell_level (n1) == level) {
	      face[len].c1 = n1;
	      face[len].c2 = child.c[1];
	      face[len++].c3 = child.c[3];
	    }
	    d1 = child.c[(ftt_cell_level (n1) == level ||
			  FTT_CELL_ID (cell) == 1) ? 1 : 3];
	    d5 = ftt_cell_neighbor_not_periodic (d1, 4);
	    if (d5)
	      CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d5));
	  }
	  else if (ftt_cell_level (d1) == level ||
		   FTT_CELL_ID (cell) == 1 ||
		   FTT_CELL_ID (cell) == 3) {
	    d5 = ftt_cell_neighbor_not_periodic (d1, 4);
	    if (d5 && !FTT_CELL_IS_LEAF (d5)) {
	      FttCellChildren child;
	    
	      ftt_cell_children_direction (d5, 5, &child);
	      d5 = child.c[(ftt_cell_level (d1) == level ||
			    FTT_CELL_ID (cell) == 1) ? 2 : 0];
	    }
	  }
	}
      }
      if (ftt_cell_level (n1) == level ||
	  FTT_CELL_ID (cell) == 0 ||
	  FTT_CELL_ID (cell) == 1) {
	d3 = ftt_cell_neighbor_not_periodic (n1, 4);
	if (d3) {
	  FttCell * d = NULL;

	  if (!FTT_CELL_IS_LEAF (d3)) {
	    FttCellChildren child;
	    
	    ftt_cell_children_direction (d3, 5, &child);
	    if (ftt_cell_level (n1) == level) {
	      face[len].c1 = n1;
	      face[len].c2 = child.c[2];
	      face[len++].c3 = child.c[3];
	    }
	    d3 = child.c[(ftt_cell_level (n1) == level ||
			  FTT_CELL_ID (cell) == 1) ? 3 : 2];
	    d = ftt_cell_neighbor_not_periodic (d3, 0);
	  }
	  else if (ftt_cell_level (d3) == level ||
		   FTT_CELL_ID (cell) == 1 ||
		   FTT_CELL_ID (cell) == 3) {
	    d = ftt_cell_neighbor_not_periodic (d3, 0);
	    if (d && !FTT_CELL_IS_LEAF (d)) {
	      FttCellChildren child;
	      
	      ftt_cell_children_direction (d, 1, &child);
	      d = child.c[(ftt_cell_level (d3) == level ||
			   FTT_CELL_ID (cell) == 1) ? 3 : 2];
	    }
	  }
	  if (d) {
	    CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d));
	    CELL_FACES_ASSERT (!d5 || d == d5);
	    d5 = d;
	  }
	}
      }
    }
    else { /* !FTT_CELL_IS_LEAF (n1) */
      FttCellChildren child;
      FttCell * d;

      ftt_cell_children_direction (n1, 3, &child);
      face[len].c1 = child.c[0];
      face[len].c2 = child.c[1];
      face[len++].c3 = child.c[2];
      face[len].c1 = child.c[1];
      face[len].c2 = child.c[3];
      face[len++].c3 = child.c[2];
      n1 = child.c[1];

      d1 = ftt_cell_neighbor_not_periodic (n1, 0);
      if (d1) {
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d1));
	face[len].c1 = n1;
	face[len].c2 = d1;
	face[len++].c3 = child.c[3];
	d5 = ftt_cell_neighbor_not_periodic (d1, 4);
	if (d5 && !FTT_CELL_IS_LEAF (d5)) {
	  FttCellChildren child;

	  ftt_cell_children_direction (d5, 5, &child);
	  d5 = child.c[2];
	}
      }
      d = ftt_cell_neighbor_not_periodic (child.c[3], 0);
      if (d && d != d1) {
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d));
	face[len].c1 = child.c[3];
	face[len].c2 = d1;
	face[len++].c3 = d;
      }

      d3 = ftt_cell_neighbor_not_periodic (n1, 4);
      if (d3) {
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d3));
	face[len].c1 = n1;
	face[len].c2 = child.c[0];
	face[len++].c3 = d3;
	d = ftt_cell_neighbor_not_periodic (d3, 0);
	if (d) {
	  if (!FTT_CELL_IS_LEAF (d)) {
	    FttCellChildren child;

	    ftt_cell_children_direction (d, 1, &child);
	    d = child.c[3];
	  }
	  CELL_FACES_ASSERT (!d5 || d == d5);
	  d5 = d;
	}
      }
      d = ftt_cell_neighbor_not_periodic (child.c[0], 4);
      if (d && d != d3) {
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d));
	face[len].c1 = child.c[0];
	face[len].c2 = d;
	face[len++].c3 = d3;
      }
    }

    if (d1 && d5 && d1 != d5) {
      face[len].c1 = d1;
      face[len].c2 = n1;
      face[len++].c3 = d5;
    }
    if (d3 && d5 && d3 != d5) {
      face[len].c1 = n1;
      face[len].c2 = d3;
      face[len++].c3 = d5;
    }
  }  

  /* ----------------------------------------------------------------------*/

  if ((n2 = neighbor.c[0])) {
    if (FTT_CELL_IS_LEAF (n2)) {
      if (ftt_cell_level (n2) == level ||
	  FTT_CELL_ID (cell) == 1 ||
	  FTT_CELL_ID (cell) == 5) {
	d2 = ftt_cell_neighbor_not_periodic (n2, 2);
	if (d2) {
	  if (!FTT_CELL_IS_LEAF (d2)) {
	    FttCellChildren child;
	    
	    ftt_cell_children_direction (d2, 3, &child);
	    if (ftt_cell_level (n2) == level) {
	      face[len].c1 = n2;
	      face[len].c2 = child.c[2];
	      face[len++].c3 = child.c[0];
	    }
	    d2 = child.c[(ftt_cell_level (n2) == level ||
			  FTT_CELL_ID (cell) == 1) ? 0 : 2];
	    d6 = ftt_cell_neighbor_not_periodic (d2, 4);
	  }
	  else if (ftt_cell_level (d2) == level ||
		   FTT_CELL_ID (cell) == 1 ||
		   FTT_CELL_ID (cell) == 0) {
	    d6 = ftt_cell_neighbor_not_periodic (d2, 4);
	    if (d6 && !FTT_CELL_IS_LEAF (d6)) {
	      FttCellChildren child;
	    
	      ftt_cell_children_direction (d6, 5, &child);
	      d6 = child.c[(ftt_cell_level (d2) == level ||
			    FTT_CELL_ID (cell) == 1) ? 2 : 3];
	    }
	  }
	  CELL_FACES_ASSERT (!d1 || d2 == d1);
	}
      }
      if (ftt_cell_level (n2) == level ||
	  FTT_CELL_ID (cell) == 1 ||
	  FTT_CELL_ID (cell) == 3) {
	d7 = ftt_cell_neighbor_not_periodic (n2, 4);
	if (d7) {
	  FttCell * d = NULL;

	  if (!FTT_CELL_IS_LEAF (d7)) {
	    FttCellChildren child;
	    
	    ftt_cell_children_direction (d7, 5, &child);
	    if (ftt_cell_level (n2) == level) {
	      face[len].c1 = n2;
	      face[len].c2 = child.c[0];
	      face[len++].c3 = child.c[2];
	    }
	    d7 = child.c[(ftt_cell_level (n2) == level ||
			  FTT_CELL_ID (cell) == 1) ? 0 : 2];
	    d = ftt_cell_neighbor_not_periodic (d7, 2);
	  }
	  else if (ftt_cell_level (d7) == level ||
		   FTT_CELL_ID (cell) == 1 ||
		   FTT_CELL_ID (cell) == 0) {
	    d = ftt_cell_neighbor_not_periodic (d7, 2);
	    if (d && !FTT_CELL_IS_LEAF (d)) {
	      FttCellChildren child;
	      
	      ftt_cell_children_direction (d, 3, &child);
	      d = child.c[(ftt_cell_level (d7) == level ||
			   FTT_CELL_ID (cell) == 1) ? 2 : 3];
	    }
	  }
	  if (d) {
	    CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d));
	    CELL_FACES_ASSERT (!d6 || d == d6);
	    d6 = d;
	  }
	}
      }
      if (d6)
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d6) && (!d5 || d5 == d6));
    }
    else { /* !FTT_CELL_IS_LEAF (n2) */
      FttCellChildren child;
      FttCell * d;

      ftt_cell_children_direction (n2, 1, &child);
      face[len].c1 = child.c[0];
      face[len].c2 = child.c[1];
      face[len++].c3 = child.c[3];
      face[len].c1 = child.c[0];
      face[len].c2 = child.c[3];
      face[len++].c3 = child.c[2];
      n2 = child.c[0];

      d2 = ftt_cell_neighbor_not_periodic (n2, 2);
      if (d2) {
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d2));
	face[len].c1 = n2;
	face[len].c2 = child.c[2];
	face[len++].c3 = d2;
	d6 = ftt_cell_neighbor_not_periodic (d2, 4);
	if (d6 && !FTT_CELL_IS_LEAF (d6)) {
	  FttCellChildren child;

	  ftt_cell_children_direction (d6, 5, &child);
	  d6 = child.c[2];
	}
      }
      d = ftt_cell_neighbor_not_periodic (child.c[2], 2);
      if (d && d != d2) {
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d));
	face[len].c1 = child.c[2];
	face[len].c2 = d;
	face[len++].c3 = d2;
      }

      d7 = ftt_cell_neighbor_not_periodic (n2, 4);
      if (d7) {
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d7));
	face[len].c1 = n2;
	face[len].c2 = d7;
	face[len++].c3 = child.c[1];
	d = ftt_cell_neighbor_not_periodic (d7, 2);
	if (d) {
	  if (!FTT_CELL_IS_LEAF (d)) {
	    FttCellChildren child;

	    ftt_cell_children_direction (d, 3, &child);
	    d = child.c[2];
	  }
	  CELL_FACES_ASSERT (!d6 || d == d6);
	  d6 = d;
	}
      }
      d = ftt_cell_neighbor_not_periodic (child.c[1], 4);
      if (d && d != d7) {
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d));
	face[len].c1 = child.c[1];
	face[len].c2 = d7;
	face[len++].c3 = d;
      }
      if (d6)
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d6) && (!d5 || d5 == d6));
    }

    if (d2 && d6 && d2 != d6) {
      face[len].c1 = d2;
      face[len].c2 = d6;
      face[len++].c3 = n2;
    }
    if (d7 && d6 && d7 != d6) {
      face[len].c1 = d6;
      face[len].c2 = d7;
      face[len++].c3 = n2;
    }
  }

  /* ----------------------------------------------------------------------*/

  if ((n3 = neighbor.c[4])) {
    if (FTT_CELL_IS_LEAF (n3)) {
      if (ftt_cell_level (n3) == level ||
	  FTT_CELL_ID (cell) == 1 ||
	  FTT_CELL_ID (cell) == 0) {
	d4 = ftt_cell_neighbor_not_periodic (n3, 2);
	if (d4) {
	  if (!FTT_CELL_IS_LEAF (d4)) {
	    FttCellChildren child;
	    
	    ftt_cell_children_direction (d4, 3, &child);
	    if (ftt_cell_level (n3) == level) {
	      face[len].c1 = n3;
	      face[len].c2 = child.c[3];
	      face[len++].c3 = child.c[2];
	    }
	    d4 = child.c[(ftt_cell_level (n3) == level ||
			  FTT_CELL_ID (cell) == 1) ? 3 : 2];
	    d9 = ftt_cell_neighbor_not_periodic (d4, 0);
	  }
	  else if (ftt_cell_level (d4) == level ||
		   FTT_CELL_ID (cell) == 1 ||
		   FTT_CELL_ID (cell) == 5) {
	    d9 = ftt_cell_neighbor_not_periodic (d4, 0);
	    if (d9 && !FTT_CELL_IS_LEAF (d9)) {
	      FttCellChildren child;
	    
	      ftt_cell_children_direction (d9, 1, &child);
	      d9 = child.c[(ftt_cell_level (d4) == level ||
			    FTT_CELL_ID (cell) == 1) ? 3 : 1];
	    }
	  }
	  CELL_FACES_ASSERT (!d3 || d4 == d3);
	}
      }
      if (ftt_cell_level (n3) == level ||
	  FTT_CELL_ID (cell) == 1 ||
	  FTT_CELL_ID (cell) == 3) {
	d8 = ftt_cell_neighbor_not_periodic (n3, 0);
	if (d8) {
	  FttCell * d = NULL;

	  if (!FTT_CELL_IS_LEAF (d8)) {
	    FttCellChildren child;
	    
	    ftt_cell_children_direction (d8, 1, &child);
	    if (ftt_cell_level (n3) == level) {
	      face[len].c1 = n3;
	      face[len].c2 = child.c[3];
	      face[len++].c3 = child.c[2];
	    }
	    d8 = child.c[(ftt_cell_level (n3) == level ||
			  FTT_CELL_ID (cell) == 1) ? 2 : 3];
	    d = ftt_cell_neighbor_not_periodic (d8, 2);
	  }
	  else if (ftt_cell_level (d8) == level ||
		   FTT_CELL_ID (cell) == 1 ||
		   FTT_CELL_ID (cell) == 5) {
	    d = ftt_cell_neighbor_not_periodic (d8, 2);
	    if (d && !FTT_CELL_IS_LEAF (d)) {
	      FttCellChildren child;
	      
	      ftt_cell_children_direction (d, 3, &child);
	      d = child.c[(ftt_cell_level (d8) == level ||
			   FTT_CELL_ID (cell) == 1) ? 2 : 0];
	    }
	  }
	  if (d) {
	    CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d));
	    CELL_FACES_ASSERT (!d9 || d == d9);
	    d9 = d;
	  }
	  CELL_FACES_ASSERT (!d7 || d7 == d8);
	}
      }
      if (d9)
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d9) && 
		  (!d5 || d5 == d9) &&
		  (!d6 || d6 == d9));
    }
    else { /* !FTT_CELL_IS_LEAF (n3) */
      FttCellChildren child;
      FttCell * d;

      ftt_cell_children_direction (n3, 5, &child);
      face[len].c1 = child.c[0];
      face[len].c2 = child.c[2];
      face[len++].c3 = child.c[1];
      face[len].c1 = child.c[2];
      face[len].c2 = child.c[3];
      face[len++].c3 = child.c[1];
      n3 = child.c[1];

      d4 = ftt_cell_neighbor_not_periodic (n3, 2);
      if (d4) {
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d4));
	face[len].c1 = n3;
	face[len].c2 = d4;
	face[len++].c3 = child.c[0];
	d9 = ftt_cell_neighbor_not_periodic (d4, 0);
	if (d9 && !FTT_CELL_IS_LEAF (d9)) {
	  FttCellChildren child;

	  ftt_cell_children_direction (d9, 1, &child);
	  d9 = child.c[3];
	}
	CELL_FACES_ASSERT (!d3 || d4 == d3);
      }
      d = ftt_cell_neighbor_not_periodic (child.c[0], 2);
      if (d && d != d4) {
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d));
	face[len].c1 = child.c[0];
	face[len].c2 = d4;
	face[len++].c3 = d;
      }

      d8 = ftt_cell_neighbor_not_periodic (n3, 0);
      if (d8) {
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d8));
	face[len].c1 = n3;
	face[len].c2 = child.c[3];
	face[len++].c3 = d8;
	d = ftt_cell_neighbor_not_periodic (d8, 2);
	if (d) {
	  if (!FTT_CELL_IS_LEAF (d)) {
	    FttCellChildren child;

	    ftt_cell_children_direction (d, 3, &child);
	    d = child.c[2];
	  }
	  CELL_FACES_ASSERT (!d9 || d == d9);
	  d9 = d;
	}
	CELL_FACES_ASSERT (!d7 || d7 == d8);
      }
      d = ftt_cell_neighbor_not_periodic (child.c[3], 0);
      if (d && d != d8) {
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d));
	face[len].c1 = child.c[3];
	face[len].c2 = d;
	face[len++].c3 = d8;
      }
      if (d9)
	CELL_FACES_ASSERT (FTT_CELL_IS_LEAF (d9) && 
			   (!d5 || d5 == d9) &&
			   (!d6 || d6 == d9));
    }

    if (d4 && d9 && d4 != d9) {
      face[len].c1 = d9;
      face[len].c2 = d4;
      face[len++].c3 = n3;
    }
    if (d8 && d9 && d8 != d9) {
      face[len].c1 = d8;
      face[len].c2 = d9;
      face[len++].c3 = n3;
    }
  }

  return len;
}

#endif /* not FTT_2D */

typedef struct {
  Array data;
  Field f;
  guint i, dim;
  GfsVariable ** v;
  gchar * name;
} GfsArray;

static GfsArray * gfs_array_new (GfsVariable ** v,
				 gchar * name,
				 Type t, 
				 guint dim,
				 gboolean field)
{
  GfsArray * a;

  a = g_malloc0 (sizeof (GfsArray));
  if (!(a->data = DXNewArray (t, CATEGORY_REAL, 1, dim))) {
    g_free (a);
    return NULL;
  }
  if (field) {
    if (!(a->f = DXNewField ())) {
      DXDelete ((Object) a->data);
      g_free (a);
      return NULL;
    }
    if (!DXSetStringAttribute ((Object) a->data, "dep", "positions")) {
      DXDelete ((Object) a->data);
      DXDelete ((Object) a->f);
      g_free (a);
      return NULL;
    }
    if (!DXSetComponentValue (a->f, "data", (Object) a->data)) {
      DXDelete ((Object) a->data);
      DXDelete ((Object) a->f);
      g_free (a);
      return NULL;
    }
  }
  a->dim = dim;
  if (v) {
    a->v = g_malloc (dim*sizeof (GfsVariable *));
    memcpy (a->v, v, dim*sizeof (GfsVariable *));
  }
  if (name)
    a->name = g_strdup (name);
  return a;
}

static void gfs_array_add (GfsArray * a, gpointer val)
{
  DXAddArrayData (a->data, a->i++, 1, val);
}

static void gfs_array_destroy (GfsArray * a, gboolean free_seg)
{
  if (free_seg) {
    if (a->f)
      DXDelete ((Object) a->f);
    else
      DXDelete ((Object) a->data);
  }
  if (a->v)
    g_free (a->v);
  if (a->name)
    g_free (a->name);
  g_free (a);
}

static gint cell_index (FttCell * cell,
			GfsArray * pos,
			GPtrArray * data,
			FttVector * lambda)
{
  if (GFS_STATE (cell)->dp <= 0.) {
    gfloat fp[3];
    FttVector p;
    guint i;

    for (i = 0; i < data->len; i++) {
      GfsArray * a = data->pdata[i];
      guint j;
      gfloat val[3];

      if (a->dim == FTT_DIMENSION)
	for (j = 0; j < a->dim; j++)
	  val[j] = GFS_VARIABLE (cell, a->v[j]->i)/((gdouble *) lambda)[j];
      else if (a->dim == 1)
	val[0] = GFS_VARIABLE (cell, a->v[0]->i);
      else
	g_assert_not_reached ();
      gfs_array_add (a, val);
    }
    gfs_cell_cm (cell, &p);
    fp[0] = p.x/lambda->x;
    fp[1] = p.y/lambda->y;
    fp[2] = p.z/lambda->z;
    gfs_array_add (pos, fp);
    GFS_STATE (cell)->dp = pos->i;
  }
  return GFS_STATE (cell)->dp - 1.;
}

#if FTT_2D

static void add_cell (FttCell * cell, gpointer * par)
{
  GfsArray * pos = par[0];
  GfsArray * con = par[1];
  GPtrArray * data = par[2];
  FttVector * lambda = par[3];
  Edge e[12];
  guint i, len;
  gint index[3];

  len = cell_edges (cell, e);
  index[0] = cell_index (cell, pos, data, lambda);
  for (i = 0; i < len; i++) {
    index[1] = cell_index (e[i].c1, pos, data, lambda);
    index[2] = cell_index (e[i].c2, pos, data, lambda);
    gfs_array_add (con, index);
  }
}

#else /* not FTT_2D */

static void add_cell (FttCell * cell, gpointer * par)
{
  GfsArray * pos = par[0];
  GfsArray * con = par[1];
  GPtrArray * data = par[2];
  FttVector * lambda = par[3];
  Face f[24];
  guint i, len;
  gint index[4];

  len = cell_faces (cell, f);
  index[0] = cell_index (cell, pos, data, lambda);
  for (i = 0; i < len; i++) {
    index[1] = cell_index (f[i].c1, pos, data, lambda);
    index[2] = cell_index (f[i].c2, pos, data, lambda);
    index[3] = cell_index (f[i].c3, pos, data, lambda);
    gfs_array_add (con, index);
  }
}

#endif /* not FTT_2D */

static void add_solid_vertex (GtsPoint * p, gpointer * par)
{
  gfloat fp[3];
  GfsArray * pos = par[0];
  FttVector * lambda = par[1];
  
  GTS_OBJECT (p)->reserved = GUINT_TO_POINTER (pos->i);
  fp[0] = p->x/lambda->x;
  fp[1] = p->y/lambda->y;
  fp[2] = p->z/lambda->z;
  gfs_array_add (pos, fp);
}

static void add_solid_face (GtsTriangle * t, GfsArray * con)
{
  gint index[3];
  GtsVertex * v1, * v2, * v3;

  gts_triangle_vertices (t, &v1, &v2, &v3);
  index[0] = GPOINTER_TO_UINT (GTS_OBJECT (v1)->reserved);
  index[1] = GPOINTER_TO_UINT (GTS_OBJECT (v2)->reserved);
  index[2] = GPOINTER_TO_UINT (GTS_OBJECT (v3)->reserved);
  gfs_array_add (con, index);
}

#if FTT_2D
Error m_ImportGfs2D (Object * in, Object * out);

Error m_ImportGfs2D (Object * in, Object * out)
#else /* not FTT_2D */
Error m_ImportGfs3D (Object * in, Object * out);

Error m_ImportGfs3D (Object * in, Object * out)
#endif /* not FTT_2D */
{
  char * filename;
  FILE * fp = NULL;
  GtsFile * f = NULL;
  GfsSimulation * sim = NULL;
  GfsDomain * domain;
  gpointer par[4];
  GfsArray * pos = NULL, * con = NULL;
  GPtrArray * data = NULL;
  Group group = NULL;
  GfsVariable * v;
  guint i;
  gboolean pos_used = FALSE, con_used = FALSE;
  Field solid = NULL;

  /* extract the file name from in[0] */
  if (!in[0]) {
    DXSetError (ERROR_BAD_PARAMETER, "missing filename");
    goto error;
  }
  else if (!DXExtractString (in[0], &filename)) {
    DXSetError (ERROR_BAD_PARAMETER, "filename must be a string");
    goto error;
  }

  /* check to see that the file is accessible, and is a GFS file */
  fp = fopen (filename, "rt");
  if (!fp) {
    DXSetError (ERROR_BAD_PARAMETER,
		"file \"%s\" is not accessible", 
		filename);
     goto error;
  }
  f = gts_file_new (fp);

  gfs_init (NULL, NULL);

  if ((sim = gfs_simulation_read (f)) == NULL) {
    DXSetError (ERROR_BAD_PARAMETER,
		"file \"%s\" is not a valid GFS file\n"
		"%s:%d:%d: %s",
		filename, filename, f->line, f->pos, f->error);
    goto error;
  }
  gts_file_destroy (f);
  f = NULL;
  fclose (fp);
  fp = NULL;

  domain = GFS_DOMAIN (sim);
    
  gfs_domain_match (domain);
  v = domain->variables;
  while (v) {
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, v);
    v = v->next;
  }

  /* positions */
  if (!(pos = gfs_array_new (NULL, NULL, TYPE_FLOAT, FTT_DIMENSION, FALSE)))
    goto error;

  /* connections */
  if (!(con = gfs_array_new (NULL, NULL, TYPE_INT, FTT_DIMENSION + 1, FALSE)))
    goto error;

  /* data */
  data = g_ptr_array_new ();
  v = domain->variables_io;
  while (v) {
    GfsVariable * vv[FTT_DIMENSION], * v1;

    if (!strcmp (v->name, "U") && v->next &&
	!strcmp (v->next->name, "V")
#if (!FTT_2D)
	&& v->next->next && !strcmp (v->next->next->name, "W")
#endif /* not FTT_2D */
	) {
      vv[0] = v; v = v->next;
      vv[1] = v; 
#if (!FTT_2D)
      v = v->next;
      vv[2] = v;
#endif /* not FTT_2D */
      g_ptr_array_add (data, gfs_array_new (vv, "U", TYPE_FLOAT, 
					    FTT_DIMENSION, TRUE));
    }
    else if (v->name[strlen (v->name) - 1] == 'x' &&
	     (vv[0] = v) && (v1 = v->next) &&
	     v1->name[strlen (v1->name) - 1] == 'y' &&
	     (vv[1] = v1) 
#if (!FTT_2D)
	     && (v1 = v1->next) &&
	     v1->name[strlen (v1->name) - 1] == 'z' &&
	     (vv[2] = v1)
#endif /* not FTT_2D */
	     ) {
      gchar * name = g_strndup (v->name, strlen (v->name) - 1);
      g_ptr_array_add (data, gfs_array_new (vv, name, TYPE_FLOAT, 
					    FTT_DIMENSION, TRUE));
      g_free (name);
      v = v1;
    }
    else 
      g_ptr_array_add (data, gfs_array_new (&v, v->name, TYPE_FLOAT, 1, TRUE));
    v = v->next;
  }

  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, gfs_dp);
  par[0] = pos;
  par[1] = con;
  par[2] = data;
  par[3] = &domain->lambda;
  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) add_cell, par);

  if (!(group = DXNewGroup ()))
    goto error;

  if (!DXSetStringAttribute ((Object) con->data, "element type", 
			     FTT_DIMENSION == 2 ? "triangles" : "tetrahedra"))
    goto error;

  for (i = 0; i < data->len; i++) {
    GfsArray * a = data->pdata[i];

    if (!DXSetComponentValue (a->f, "positions", (Object) pos->data))
      goto error;
    pos_used = TRUE;
    if (!DXSetComponentValue (a->f, "connections", (Object) con->data))
      goto error;
    con_used = TRUE;
    if (!DXSetMember (group, a->name, (Object) a->f))
      goto error;
    gfs_array_destroy (a, FALSE);
    data->pdata[i] = NULL;
  }
  g_ptr_array_free (data, TRUE);
  data = NULL;
  gfs_array_destroy (pos, FALSE);
  pos = NULL;
  gfs_array_destroy (con, FALSE);
  con = NULL;

  if (sim->surface) {
    if (!(solid = DXNewField ()))
      goto error;
    if (!(pos = gfs_array_new (NULL, NULL, TYPE_FLOAT, 3, FALSE)))
      goto error;
    pos_used = FALSE;
    if (!(con = gfs_array_new (NULL, NULL, TYPE_INT, 3, FALSE)))
      goto error;
    con_used = FALSE;
    par[0] = pos;
    par[1] = &domain->lambda;
    gts_surface_foreach_vertex (sim->surface, (GtsFunc) add_solid_vertex, par);
    gts_surface_foreach_face (sim->surface, (GtsFunc) add_solid_face, con);
    gts_surface_foreach_vertex (sim->surface, 
				(GtsFunc) gts_object_reset_reserved, NULL);
    if (!DXSetStringAttribute ((Object) con->data, 
			       "element type", "triangles"))
      goto error;
    if (!DXSetComponentValue (solid, "positions", (Object) pos->data))
      goto error;
    gfs_array_destroy (pos, FALSE);
    pos = NULL;
    if (!DXSetComponentValue (solid, "connections", (Object) con->data))
      goto error;
    gfs_array_destroy (con, FALSE);
    con = NULL;
    if (!DXSetMember (group, "solid", (Object) solid))
      goto error;
    solid = NULL;
  }
  gts_object_destroy (GTS_OBJECT (sim));
  sim = NULL;

  if (!DXEndObject ((Object) group))
    goto error;

  out[0] = (Object) group;
  return OK;

 error:
  if (fp) fclose (fp);
  if (f) gts_file_destroy (f);
  if (sim) gts_object_destroy (GTS_OBJECT (sim));
  if (data) {
    guint i;

    for (i = 0; i < data->len; i++) {
      GfsArray * a = data->pdata[i];

      if (a) gfs_array_destroy (a, TRUE);
    }
    g_ptr_array_free (data, TRUE);
  }
  if (pos) gfs_array_destroy (pos, !pos_used);
  if (con) gfs_array_destroy (con, !con_used);
  if (solid) DXDelete ((Object) solid);
  if (group) DXDelete ((Object) group);
  return ERROR;
}


