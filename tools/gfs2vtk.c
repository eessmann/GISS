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
#include "config.h"
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif /* HAVE_UNISTD_H */

#include "init.h"
#include "simulation.h"
#include "graphic.h"
#include "solid.h"

#if (!FTT_2D)
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

static void set_cell_index (FttCell * cell, guint * np)
{
  GFS_STATE (cell)->dp = (*np)++;
}

static void write_pos (FttCell * cell, FILE * fp)
{
  FttVector p;

  gfs_cell_cm (cell, &p);
  fprintf (fp, "%g %g %g\n", p.x, p.y, p.z);
}

#define cell_index(c) ((guint) GFS_STATE (c)->dp)

static void add_cell_tetras (FttCell * cell, GArray * tetra)
{
  Face f[24];
  guint i, len;
  guint index[4];

  len = cell_faces (cell, f);
  index[0] = cell_index (cell);
  for (i = 0; i < len; i++) {
    index[1] = cell_index (f[i].c1);
    index[2] = cell_index (f[i].c2);
    index[3] = cell_index (f[i].c3);
    g_array_append_val (tetra, (index[0]));
  }
}

static void write_vector (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  FILE * fp = data[1];

  fprintf (fp, "%g %g %g\n",
	   GFS_VARIABLE (cell, v->i),
	   GFS_VARIABLE (cell, v->next->i),
	   GFS_VARIABLE (cell, v->next->next->i));
}

static void write_scalar (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  FILE * fp = data[1];

  fprintf (fp, "%g\n", GFS_VARIABLE (cell, v->i));
}

static void gfs_write_vtk_3D (GfsDomain * domain, FILE * fp)
{
  guint np = 0, i;
  GArray * tetra;
  GfsVariable * v;

  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) set_cell_index, &np);
  fprintf (fp, 
	   "# vtk DataFile Version 2.0\n"
	   "Generated by Gerris\n"
           "ASCII\n"
	   "DATASET UNSTRUCTURED_GRID\n"
	   "POINTS %u float\n",
	   np);
  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) write_pos, fp);
  tetra = g_array_new (FALSE, FALSE, 4*sizeof (guint));
  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) add_cell_tetras, tetra);
  fprintf (fp, "CELLS %u %u\n", tetra->len, 5*tetra->len);
  for (i = 0; i < tetra->len; i++) {
    guint * index = &g_array_index (tetra, guint, 4*i);

    fprintf (fp, "4 %u %u %u %u\n", index[0], index[1], index[2], index[3]);
  }
  fprintf (fp, "CELL_TYPES %u\n", tetra->len);
  for (i = 0; i < tetra->len; i++)
    fputs ("10\n", fp);
  g_array_free (tetra, TRUE);
  fprintf (fp, "POINT_DATA %u\n", np);
  v = domain->variables;
  while (v) {
    gpointer data[2];

    if (!strcmp (v->name, "U") && v->next &&
	!strcmp (v->next->name, "V") && v->next->next &&
	!strcmp (v->next->next->name, "W")) {
      fputs ("VECTORS U float\n", fp);
      data[0] = v;
      data[1] = fp;
      gfs_domain_cell_traverse (domain,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) write_vector, data);
      v = v->next->next->next;
    }
    else if (v->name[strlen (v->name) - 1] == 'x' && v->next &&
	     v->next->name[strlen (v->name) - 1] == 'y' && v->next->next &&
	     v->next->next->name[strlen (v->name) - 1] == 'z') {
      gchar * name = g_strndup (v->name, strlen (v->name) - 1);

      fprintf (fp, "VECTORS %s float\n", name);
      g_free (name);
      data[0] = v;
      data[1] = fp;
      gfs_domain_cell_traverse (domain,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) write_vector, data);
      v = v->next->next->next;
    }
    else {
      fprintf (fp, 
	       "SCALARS %s float 1\n"
	       "LOOKUP_TABLE default\n",
	       v->name);
      data[0] = v;
      data[1] = fp;
      gfs_domain_cell_traverse (domain,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) write_scalar, data);
      v = v->next;
    }
  }
}


#endif /* not FTT_2D */

int main (int argc, char * argv[])
{
  int c = 0;
  gboolean verbose = FALSE;
  GtsFile * fp;
  GfsSimulation * simulation;
  GfsDomain * domain;
  GfsVariable * v;
  FILE * field, * bound;
  gchar * fname;

  gfs_init (&argc, &argv);

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"help", no_argument, NULL, 'h'},
      {"verbose", no_argument, NULL, 'v'},
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "hv", 
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "hv"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'v': /* verbose */
      verbose = TRUE;
      break;
    case 'h': /* help */
      fprintf (stderr,
     "Usage: gfs2vtk [OPTION] FILE < GFS_FILE\n"
     "Converts a Gerris simulation file to VTK format. Two files are generated:\n"
     "  FILE_field.vtk: contains field data\n"
     "  FILE_bound.vtk: contains solid boundaries\n"
     "\n"
     "  -v      --verbose     display statistics and other info\n"
     "  -h      --help        display this help and exit\n"
     "\n"
     "Reports bugs to %s\n",
	       FTT_MAINTAINER);
      return 0; /* success */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `gfs2vtk --help' for more information.\n");
      return 1; /* failure */
    }
  }

  if (optind >= argc) {
    fprintf (stderr, 
	     "gfs2vtk: missing FILE\n"
	     "Try `gfs2vtk --help' for more information.\n");
    return 1;
  }
  fname = g_strconcat (argv[optind], "_field.vtk", NULL);
  field = fopen (fname, "wt");
  g_free (fname);
  if (field == NULL) {
    fprintf (stderr, 
	     "gfs2vtk: cannot open file `%s_field.vtk'\n",
	     argv[optind]);
    perror ("");
    return 1;
  }

  fp = gts_file_new (stdin);
  if (!(simulation = gfs_simulation_read (fp))) {
    fprintf (stderr, 
	     "gfs2vtk: file on standard input is not a valid simulation file\n"
	     "<stdin>:%d:%d: %s\n",
	     fp->line, fp->pos, fp->error);
    return 1;
  }

  domain = GFS_DOMAIN (simulation);

  gfs_domain_match (domain);

  v = domain->variables;
  while (v) {
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, v);
    v = v->next;
  }

#if FTT_2D
  gfs_write_vtk (domain, FTT_TRAVERSE_LEAFS, -1, field);
#else /* 3D */
  gfs_write_vtk_3D (domain, field);
#endif /* 3D */

  fclose (field);

  if (simulation->surface) {
    fname = g_strconcat (argv[optind], "_bound.vtk", NULL);
    bound = fopen (fname, "wt");
    g_free (fname);
    if (bound == NULL) {
      fprintf (stderr, 
	       "gfs2vtk: cannot open file `%s_bound.vtk'\n",
	       argv[optind]);
      perror ("");
      return 1;
    }
    gts_surface_write_vtk (simulation->surface, bound);
    fclose (bound);
  }

  gts_object_destroy (GTS_OBJECT (simulation));
  gts_file_destroy (fp);

  return 0;
}
