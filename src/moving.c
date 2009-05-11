/* Gerris - The GNU Flow Solver
 * Copyright (C) 2005-2009 National Institute of Water and Atmospheric Research
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

#include <stdlib.h>
#include <math.h>
#include <gts.h>
#include "moving.h"
#include "simulation.h"
#include "domain.h"
#include "utils.h"
#include "ftt.h"
#include "refine.h"
#include "adaptive.h"
#include "solid.h"
#include "vof.h"
#include "surface.h"
#include "advection.h"
#include "source.h"

typedef struct {
  GfsSimulation * sim;
  GfsSolidMoving * s;
  GfsVariable * old_solid_v, ** sold2, ** v;
} SolidInfo;

static double surface_value (FttCell * cell, GfsVariable * v, FttVector * ca)
{
  gdouble val = 0.;
  if (!v->surface_bc)
    /* default surface BC for velocity is zero */
    return 0.;
  else if (GFS_STATE (cell)->solid) {
    FttVector oldca;
    if (ca) {
      oldca = GFS_STATE (cell)->solid->ca;
      GFS_STATE (cell)->solid->ca = *ca;
    }
    (* GFS_SURFACE_GENERIC_BC_CLASS (GTS_OBJECT (v->surface_bc)->klass)->bc) (cell, v->surface_bc);
    if (ca)
      GFS_STATE (cell)->solid->ca = oldca;
    val = GFS_STATE (cell)->solid->fv;
  }
  else {
    GfsSolidVector solid;
    if (ca)
      solid.ca = *ca;
    else
      ftt_cell_pos (cell, &solid.ca);
    solid.cm = solid.ca;
    GFS_STATE (cell)->solid = &solid;
    (* GFS_SURFACE_GENERIC_BC_CLASS (GTS_OBJECT (v->surface_bc)->klass)->bc) (cell, v->surface_bc);
    GFS_STATE (cell)->solid = NULL;
    val = solid.fv;
  }
  if (!(cell->flags & GFS_FLAG_DIRICHLET))
    g_assert_not_implemented ();
  return val;
}

static void init_new_cell_velocity_from_solid (FttCell * cell, SolidInfo * p)
{
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    GFS_VALUE (cell, p->v[c]) = surface_value (cell, p->v[c], NULL);
}

static void update_neighbors (FttCell * cell)
{
  gint i;
  FttCellNeighbors neighbor;
  g_assert (cell);

  ftt_cell_neighbors (cell, &neighbor);
  for (i = 0; i < FTT_NEIGHBORS; i++)
    if (neighbor.c[i] && neighbor.c[i]->children) {
      ftt_cell_neighbors_not_cached (neighbor.c[i], &(neighbor.c[i]->children->neighbors));}
}

static gboolean refine_maxlevel (FttCell * cell, gint * maxlevel)
{
  return (ftt_cell_level (cell) < *maxlevel);
}

static void moving_cell_coarse_fine (FttCell * cell, GfsVariable * v)
{
  FttCell * parent = ftt_cell_parent (cell);

  GFS_VALUE (cell, v) = GFS_VALUE (parent, v);
  if (!GFS_CELL_IS_BOUNDARY (parent)) {
    FttVector p;
    FttComponent c;

    ftt_cell_relative_pos (cell, &p);    
    for (c = 0; c < FTT_DIMENSION; c++)
      GFS_VALUE (cell, v) += (&p.x)[c]*gfs_center_van_leer_gradient (parent, c, v->i);
  }
}

#define OLD_SOLID(c) (*((GfsSolidVector **) &(GFS_VALUE (c, old_solid_v))))
#define SOLD2(c, d)  (GFS_VALUE (c, sold2[d]))

static void sold2_fine_init (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  guint n;

  ftt_cell_children (parent, &child);
  for (n = 0; n < FTT_CELLS; n++)
    if (child.c[n])
      GFS_VALUE (child.c[n], v) = 1.;
}

static void moving_cell_init (FttCell * cell, SolidInfo * solid_info)
{
  GSList * i;
  gint k;
  GfsDomain * domain = GFS_DOMAIN (solid_info->sim);
  GfsVariable * old_solid_v = GFS_MOVING_SIMULATION (domain)->old_solid;

  gfs_cell_init (cell, domain);

  i = domain->variables;
  while (i) {
    GfsVariable * v = i->data;
    
    if (v->coarse_fine == (GfsVariableFineCoarseFunc) gfs_cell_coarse_fine)
      moving_cell_coarse_fine (cell, v);
    /* these are variables specific to moving solid boudaries */
    else if (v->coarse_fine != sold2_fine_init && v != old_solid_v)
      g_assert_not_implemented ();
    i = i->next;
  }

  g_assert (OLD_SOLID (cell) == NULL);
  OLD_SOLID (cell) = g_malloc0 (sizeof (GfsSolidVector));
  OLD_SOLID (cell)->a = 0.;
  GfsVariable ** sold2 = solid_info->sold2;
  if (sold2)
    for (k = 0; k < FTT_NEIGHBORS; k++)
      SOLD2 (cell, k) = OLD_SOLID (cell)->s[k] = 0.;

  init_new_cell_velocity_from_solid (cell, solid_info);
}

static void moving_cell_fine_init (FttCell * cell, SolidInfo * solid_info)
{
  GfsDomain * domain = GFS_DOMAIN(solid_info->sim);
  GfsVariable * old_solid_v = GFS_MOVING_SIMULATION (domain)->old_solid;
  GfsVariable ** sold2 = solid_info->sold2;
  FttCellChildren child;
  guint n;
 
  gfs_cell_fine_init (cell, domain);

  /* need to update the neighbors of the "undestroyed" parent cell */
  update_neighbors (cell);
 
  ftt_cell_children (cell, &child);
  for (n = 0; n < FTT_CELLS; n++) {
    GfsSolidVector * solid = OLD_SOLID (child.c[n]);
    gint k;
    g_assert (!solid);
    solid = OLD_SOLID (child.c[n]) = g_malloc0 (sizeof (GfsSolidVector));
    solid->a = 0.;
    if (sold2)
      for (k = 0; k < FTT_NEIGHBORS; k++)
	SOLD2 (child.c[n], k) = solid->s[k] = 0.;
  }
}

static void create_new_cells (FttCell * cell, GfsSurface * s, SolidInfo * solid_info)
{
  GfsSolidMoving * solid = solid_info->s;
  gint maxlevel = gfs_function_value (solid->level, cell);

  if (FTT_CELL_IS_DESTROYED (cell) && ftt_cell_level (cell) <= maxlevel) {
    cell->flags &= ~FTT_FLAG_DESTROYED;
    moving_cell_init (cell, solid_info);
    if (ftt_cell_level (cell) < maxlevel)
      ftt_cell_refine (cell,
		       (FttCellRefineFunc) refine_maxlevel, &maxlevel,
		       (FttCellInitFunc) moving_cell_fine_init, solid_info);
  }
  else if (ftt_cell_level (cell) < maxlevel)
    ftt_cell_refine (cell,
		     (FttCellRefineFunc) refine_maxlevel, &maxlevel,
		     (FttCellInitFunc) gfs_cell_fine_init, solid_info->sim);
}

static void refine_cell_corner (FttCell * cell, GfsDomain * domain)
{
  if (ftt_refine_corner (cell))
    ftt_cell_refine_single (cell, (FttCellInitFunc) gfs_cell_fine_init, domain);
}

static void remesh_surface_moving (GfsSimulation * sim, GfsSolidMoving * s)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  SolidInfo solid_info;
  guint depth;
  gint l;

  solid_info.sim = sim;
  solid_info.s = s;
  solid_info.sold2 = GFS_MOVING_SIMULATION (sim)->sold2;
  solid_info.v = gfs_domain_velocity (domain);
  gfs_domain_traverse_cut (domain, GFS_SOLID (s)->s,
			   FTT_POST_ORDER, FTT_TRAVERSE_LEAFS | FTT_TRAVERSE_DESTROYED,
			   (FttCellTraverseCutFunc) create_new_cells, &solid_info);
  /* fixme: can this be moved outside the solid loop in move_solids()? */
  depth = gfs_domain_depth (domain);
  for (l = depth - 2; l >= 0; l--)
    gfs_domain_cell_traverse (domain,
    			      FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, l,
    			      (FttCellTraverseFunc) refine_cell_corner, domain);
}

static void solid_moving_destroy (GtsObject * object)
{
  gts_object_destroy (GTS_OBJECT (GFS_SOLID_MOVING (object)->level));
  (* GTS_OBJECT_CLASS (gfs_solid_moving_class ())->parent_class->destroy) (object);
}

static void solid_moving_read (GtsObject ** o, GtsFile * fp)
{
  GfsSolidMoving * solid = GFS_SOLID_MOVING (*o);
  
  if (GTS_OBJECT_CLASS (gfs_solid_moving_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_solid_moving_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (!GFS_IS_SURFACE (GFS_SOLID (solid)->s) || !GFS_SURFACE (GFS_SOLID (solid)->s)->s) {
    gts_file_error (fp, "moving implicit surfaces are not implemented yet");
    return;
  }

  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return;
  }
  fp->scope_max++;
  gts_file_next_token (fp);
  
  while (fp->type != GTS_ERROR && fp->type != '}') {
    if (fp->type == '\n') {
      gts_file_next_token (fp);
      continue;
    }
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a keyword");
      return;
    }
    if (!strcmp (fp->token->str, "level")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      gfs_function_read (solid->level, gfs_object_simulation (*o), fp);
    }
    else {
      gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
      return;
    }
  }
  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);
}

static void solid_moving_write (GtsObject * object, FILE * fp)
{
  GfsSolidMoving * solid = GFS_SOLID_MOVING (object);
  
  if (GTS_OBJECT_CLASS (gfs_solid_moving_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_solid_moving_class ())->parent_class->write) 
      (object, fp);
  fputs (" { level =", fp);
  gfs_function_write (solid->level, fp);
  fputc ('}', fp);
}

static int cell_is_corner (FttCell * cell)
{
  FttDirection d, d1, d2;
  gdouble  norm;
  FttCellNeighbors neighbors;
  FttVector n1, n2;


  g_assert (cell);

  ftt_cell_neighbors (cell,&neighbors);

  d1 = d2 = -1;

  if (!GFS_IS_MIXED(cell))
    return 0;

  for (d = 0; d < FTT_NEIGHBORS; d ++)
    if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0. && d1 == -1 && d2 == -1)
      d1 = d;
    else if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0 && d2 == -1)
      d2 = d;
    else if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0.)
      g_assert_not_reached ();

  if ( d1 == -1 || d2 == -1) {
    FttVector pos;
    ftt_cell_pos (cell,&pos);
    
    printf("REA: %f, %f \n", pos.x, pos.y);
    printf("d1: %i d2: %i  \n", d1,d2);
    
    g_assert_not_reached ();
  }

 


  gfs_solid_normal (neighbors.c[d1], &n1);
  norm= sqrt(n1.x*n1.x+n1.y*n1.y);
  n1.x /= norm;
  n1.y /= norm;
  gfs_solid_normal (neighbors.c[d2], &n2);
  norm= sqrt(n2.x*n2.x+n2.y*n2.y);
  n2.x /= norm;
  n2.y /= norm;


  if (d1/2 == d2/2)
    return 0;
  else {
    if (neighbors.c[d2])
      if ( neighbors.c[d1])
	if (GFS_IS_MIXED (neighbors.c[d2]) && GFS_IS_MIXED (neighbors.c[d1]))
	  if (fabs(n1.x*n2.x+n1.y*n2.y) < 0.70) {
	    if (GFS_STATE(neighbors.c[d1])->solid->s[d1] > 0 && GFS_STATE(neighbors.c[d1])->solid->s[d1] < 1)
	      return 1;
	    if (GFS_STATE(neighbors.c[d2])->solid->s[d2] > 0 && GFS_STATE(neighbors.c[d2])->solid->s[d2] < 1)
	      return 1;
	  }
    return 0;
  }
}

static int cell_was_corner (FttCell * cell, GfsVariable * old_solid_v, GfsVariable ** sold2)
{
  FttDirection d, d1, d2;
  
  g_assert (cell);

  d1 = d2 = -1;

  if (!OLD_SOLID (cell))
    return 0;

  for (d = 0; d < FTT_NEIGHBORS; d ++)
    if (SOLD2 (cell, d) != 1. && SOLD2 (cell, d) != 0. && d1 == -1 && d2 == -1)
      d1 = d;
    else if (SOLD2 (cell, d) != 1. && SOLD2 (cell, d) != 0 && d2 == -1)
      d2 = d;
    else if (SOLD2 (cell, d) != 1. && SOLD2 (cell, d) != 0.)
      g_assert_not_reached (); 

  if (d1/2 == d2/2)
    return 0;
  else {
    FttCellNeighbors neighbors;
    FttVector n1,n2;
    FttComponent c;
    gdouble norm;

    ftt_cell_neighbors (cell,&neighbors);

    if (neighbors.c[d1])
      for (c = 0; c < FTT_DIMENSION; c++)
	(&n1.x)[c] = (SOLD2 (neighbors.c[d1], 2*c + 1) - SOLD2 (neighbors.c[d1], 2*c));
    
    if (neighbors.c[d2])
      for (c = 0; c < FTT_DIMENSION; c++)
	(&n2.x)[c] = (SOLD2 (neighbors.c[d2], 2*c + 1) - SOLD2 (neighbors.c[d2], 2*c));
  
    norm= sqrt(n1.x*n1.x+n1.y*n1.y);
    n1.x /= norm;
    n1.y /= norm;
    norm= sqrt(n2.x*n2.x+n2.y*n2.y);
    n2.x /= norm;
    n2.y /= norm;

    
    
    if (neighbors.c[d1] && neighbors.c[d2])
      if (fabs(n1.x*n2.x+n1.y*n2.y) < 0.70) {
	if (SOLD2 (neighbors.c[d1], d1) > 0 && SOLD2 (neighbors.c[d1], d1) < 1)
	  return 1.;
	else if (SOLD2 (neighbors.c[d2], d2) > 0 && SOLD2 (neighbors.c[d2], d2) < 1)
	  return 1;
      }
    return 0;
  }
}

static double new_fluid_old_solid (FttCell * cell, FttDirection d1, 
				   GfsVariable * old_solid,
				   GfsVariable ** sold2) 
{  
  FttDirection d;
  FttCellNeighbors neighbors;
  double s1, s2;

  g_assert(cell);

  s1 = 1.-SOLD2 (cell, d1);
  ftt_cell_neighbors (cell,&neighbors);
  for  (d = 0; d < 2*FTT_DIMENSION;d++)
    if (d != 2*(d1/2) && d != 2*(d1/2)+1)
      if (neighbors.c[d])
	if(GFS_IS_MIXED(neighbors.c[d]))
	  if (!cell_is_corner (neighbors.c[d]) && 
	      !cell_was_corner (neighbors.c[d], old_solid, sold2)) {
	    if (GFS_STATE(neighbors.c[d])->solid->s[d1] != 1.) {
	      if (SOLD2 (neighbors.c[d], d1) == 0.)
		{
		  s2 = GFS_STATE(neighbors.c[d])->solid->s[d1];
		  return s2/(s1+s2);
		}		  
	    }
	  } 
  return -1.;
}

static double new_solid_old_fluid (FttCell * cell, FttDirection d1, 
				   GfsVariable * old_solid,
				   GfsVariable ** sold2) 
{  
  FttDirection d;
  FttCellNeighbors neighbors;
  double s1, s2;

  g_assert(cell);

  s1 = 1.-GFS_STATE (cell)->solid->s[d1];
		    
  ftt_cell_neighbors (cell,&neighbors);
  for  (d = 0; d < 2*FTT_DIMENSION;d++)
    if (d != 2*(d1/2) && d != 2*(d1/2)+1)
      if (neighbors.c[d])
	if (!cell_is_corner(neighbors.c[d]) && 
	    !cell_was_corner(neighbors.c[d], old_solid, sold2))
	  if (GFS_STATE(neighbors.c[d])->solid)	 
	    if (GFS_STATE(neighbors.c[d])->solid->s[d1] == 0. && SOLD2 (neighbors.c[d], d1) != 1.) {
	      
	      s2 = SOLD2 (neighbors.c[d], d1);
	      return s1/(s1+s2);
	    }
  return -1.;
}

static double new_solid_old_solid (FttCell * cell, FttDirection d1,
				   GfsVariable * old_solid,
				   GfsVariable ** sold2)
{
  FttDirection d;
  FttCellNeighbors neighbors;
  double s1, s2;

  g_assert(cell);

  s1 = GFS_STATE (cell)->solid->s[d1];
  ftt_cell_neighbors (cell,&neighbors);
  for  (d = 0; d < 2*FTT_DIMENSION;d++)
    if (d != 2*(d1/2) && d != 2*(d1/2)+1)
      if (neighbors.c[d] &&
	  !cell_is_corner(neighbors.c[d]) && 
	  !cell_was_corner(neighbors.c[d], old_solid, sold2)) {
	if ((GFS_IS_MIXED(neighbors.c[d]) && GFS_STATE(neighbors.c[d])->solid->s[d1] == 1.) || !GFS_IS_MIXED(neighbors.c[d])) {
	  if (SOLD2 (neighbors.c[d], d1) != 1.){
	    s2 = 1.-SOLD2 (neighbors.c[d], d1);
	    return s1/(s1+s2);
	  }
	}
	else if ((GFS_STATE(cell)->solid->s[d1] == 0. && GFS_IS_MIXED(neighbors.c[d])) ) {
	  s1 = SOLD2 (cell, d1);
	  s2 = 1.-GFS_STATE(neighbors.c[d])->solid->s[d1];
	  return s2/(s1+s2);
	}
      }
  return -1.;
}

static void second_order_face_fractions (FttCell * cell, GfsMovingSimulation * sim)
{
#ifndef FTT_2D /* 3D */
  g_assert_not_implemented ();
#endif

  GfsVariable * old_solid_v = sim->old_solid;
  GfsVariable ** sold2 = sim->sold2;
  gdouble dt1, dt2, dto1, dto2, s1, s2;
  gint d1, d2, d, do1, do2;
  FttCellNeighbors neighbors;

  dt1 = dt2 = dto1 = dto2 = -2;
  d1 = d2 = do1 = do2 = -1;
  s1 = s2 = -1;

  g_assert(cell);
      
  ftt_cell_neighbors (cell,&neighbors);

  if (!OLD_SOLID (cell) && !GFS_IS_MIXED(cell))
    return;

  if (!OLD_SOLID (cell)) {
    FttDirection c;
    OLD_SOLID (cell) = g_malloc0 (sizeof (GfsSolidVector));

    OLD_SOLID (cell)->a = 1.;
    for (c = 0; c < FTT_NEIGHBORS; c++)
      OLD_SOLID (cell)->s[c] = 1.;
  }

  /* Find directions of intersection */
  if (GFS_IS_MIXED(cell))
    for (d = 0; d < FTT_NEIGHBORS; d ++) {
      if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0. && d1 == -1 && d2 == -1)
	d1 = d;
      else if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0 && d2 == -1)
	d2 = d;
      else if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0.)
	g_assert_not_reached (); 
    }
  
  for (d = 0; d < FTT_NEIGHBORS; d ++) {
    if (SOLD2 (cell, d) != 1. && SOLD2 (cell, d) != 0. && do1 == -1 && do2 == -1)
      do1 = d;
    else if (SOLD2 (cell, d) != 1. && SOLD2 (cell, d) != 0 && do2 == -1)
      do2 = d;
    else if (SOLD2 (cell, d) != 1. && SOLD2 (cell, d) != 0.)
      g_assert_not_reached (); 
  }

  /* Treats easy cases */
  if (d1 != -1 && d1 == do1)
    OLD_SOLID (cell)->s[d1] = SOLD2 (cell, d1);
  if (d2 != -1 && d2 == do2)
    OLD_SOLID (cell)->s[d2] = SOLD2 (cell, d2);

  if (d1 == do1 && d2 == do2)
    return;
    
  /* Finds timescale for d1/do1 */
  if (d1 != -1) {
    if (SOLD2 (cell, d1) == 1.) {
      dt1 = new_solid_old_fluid (cell, d1, old_solid_v, sold2);
      if (dt1 == -1)
	if (neighbors.c[d1]){
	  FttDirection dop = ftt_opposite_direction[d1];
	  dt1 = new_solid_old_fluid (neighbors.c[d1], dop, old_solid_v, sold2);
	}
    }   
    else if (SOLD2 (cell, d1) == 0.){
      dt1 = new_solid_old_solid (cell, d1, old_solid_v, sold2);
    }  
  }
  
  if (do1 != -1 && do1 != d1 && do1 != d2) {
    if (GFS_IS_MIXED(cell) && GFS_STATE(cell)->solid->s[do1] == 0.)
      dto1 = new_solid_old_solid (cell, do1, old_solid_v, sold2);
    else
      dto1 = new_fluid_old_solid (cell, do1, old_solid_v, sold2);
  }

  /* Finds timescale for d2/do2 */

  if (d2 != -1) {
    if (SOLD2 (cell, d2) == 1.) {
      dt2 = new_solid_old_fluid (cell, d2, old_solid_v, sold2);
      if (dt2 == -1 && neighbors.c[d2]) {
	FttDirection dop = ftt_opposite_direction[d2];
	dt2 = new_solid_old_fluid (neighbors.c[d2], dop, old_solid_v, sold2);
      }
    }
    else if (SOLD2 (cell, d2) == 0.)
      dt2 = new_solid_old_solid (cell, d2, old_solid_v, sold2);
  }
 
  if (do2 != -1 && do2 != d1 && do2 != d2) {
    if (GFS_IS_MIXED(cell) && GFS_STATE(cell)->solid->s[do2] == 0.)
      dto2 = new_solid_old_solid (cell, do2, old_solid_v, sold2);
    else
      dto2 = new_fluid_old_solid (cell, do2, old_solid_v, sold2);
  }
 
  /* Uses time-scale from other faces if one is missing */
  if (dt1 == -1) {
    if (dto1 != -2)
      dt1 = dto1;
    else if (dt2 != -2)
      dt1 = dt2;
    else if (dto2 != -2)
      dt1 = dto2;
  }

  if (dt2 == -1) {
    if (dt1 != -2)
      dt2 = dt1;
    else if (dto2 != -2)
      dt2 = dto2;
    else if (dto1 != -2)
      dt2 = dto1;
  }

  if (dto1 == -1) {
    if (dt1 != -2)
      dto1 = dt1;
    else if (dt2 != -2)
      dto1 = dt2;
    else if (dto2 != -2)
      dto1 = dto2;
  }

  if (dto2 == -1) {
    if (dt1 != -2)
      dto2 = dt1;
    else if (dt2 != -2)
      dto2 = dt2;
    else if (dto1 != -2)
      dto2 = dto1;
  }

  /* Treats cell is corner */
  if (dt1 != -2 && dt2 != -2) {
    if (dt1 != dt2 && d1/2 != d2/2) {
      if (cell_is_corner (cell)) {
	if (dt1 < dt2)
	  dt2 = dt1;
	else
	  dt1 = dt2;
      }}}

  /* Treats cell was corner */
  if (dto1 != -2 && dto2 != -2 && 
      dto1 != dto2 && do1/2 != do2/2 &&
      cell_was_corner (cell, old_solid_v, sold2)) {
    if (dto1 < dto2)
      dto2 = dto1;
    else
      dto1 = dto2;
  }
  
  /* Compute the t^n+1/2 contribution of the face */
  if (do1 > -1)
    if (do1 != d1 && do1 != d2) {
      OLD_SOLID (cell)->s[do1]=SOLD2 (cell, do1)*(1-dto1)+dto1;
      if (neighbors.c[do1])
	if (!OLD_SOLID (neighbors.c[do1]) || !GFS_IS_MIXED(neighbors.c[do1])) {
	  if (!OLD_SOLID (neighbors.c[do1])) {
	    FttDirection c;
	    OLD_SOLID (neighbors.c[do1]) = g_malloc0 (sizeof (GfsSolidVector));
	    OLD_SOLID (neighbors.c[do1])->a = 1.;
	    for (c = 0; c < FTT_NEIGHBORS; c++)
	      OLD_SOLID (neighbors.c[do1])->s[c] = 1.;
	  }	  
	  OLD_SOLID (neighbors.c[do1])->s[ftt_opposite_direction[do1]] = SOLD2 (cell, do1)*(1-dto1)+dto1;
	}
    }

  if (do2 > -1)
    if (do2 != d1 && do2 != d2) {
      OLD_SOLID (cell)->s[do2]=SOLD2 (cell, do2)*(1-dto2)+dto2;
      if (neighbors.c[do2])
	if (!OLD_SOLID (neighbors.c[do2]) || !GFS_IS_MIXED(neighbors.c[do2])) {
	  if (!OLD_SOLID (neighbors.c[do2])) {
	    FttDirection c;
	    OLD_SOLID (neighbors.c[do2]) = g_malloc0 (sizeof (GfsSolidVector));
	    OLD_SOLID (neighbors.c[do2])->a = 1.;
	    for (c = 0; c < FTT_NEIGHBORS; c++)
	      OLD_SOLID (neighbors.c[do2])->s[c] = 1.;
	  }	  
	  OLD_SOLID (neighbors.c[do2])->s[ftt_opposite_direction[do2]] = SOLD2 (cell, do2)*(1-dto2)+dto2;
	}
    }


  if (d1 > -1) {
    if (SOLD2 (cell, d1) == 0.)
      OLD_SOLID (cell)->s[d1] = GFS_STATE(cell)->solid->s[d1]*(dt1-1.); 
    else if (SOLD2 (cell, d1) == 1.)
      OLD_SOLID (cell)->s[d1] = (dt1-1.)*GFS_STATE(cell)->solid->s[d1]+2.-dt1;
  }

  if (d2 > -1) {
    if (SOLD2 (cell, d2) == 0.)
      OLD_SOLID (cell)->s[d2] = GFS_STATE(cell)->solid->s[d2]*(dt2-1.); 
    else if (SOLD2 (cell, d2) == 1.)
      OLD_SOLID (cell)->s[d2] = (dt2-1.)*GFS_STATE(cell)->solid->s[d2]+2.-dt2;
  }

  if (d1/2 == d2/2 && do1 == -1 && do2 == -1)  /* third face has to be treated for the timescale determined on the other faces */  
    for (d = 0; d < FTT_NEIGHBORS; d ++)
      if (d/2 != d1/2 && SOLD2 (cell, d) == 0.)
	OLD_SOLID (cell)->s[d] = -1.+dt1+dt2;
    

  if (do1/2 == do2/2 && d1 == -1 && d2 == -1)
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (d/2 != do1/2 && SOLD2 (cell, d) == 0.)
	OLD_SOLID (cell)->s[d] = -1.+dto1+dto2;
}

static void set_old_solid (FttCell * cell, GfsVariable * old_solid_v)
{
  g_free (OLD_SOLID (cell));
  OLD_SOLID (cell) = GFS_STATE (cell)->solid;
  GFS_STATE (cell)->solid = NULL;
  cell->flags &= ~GFS_FLAG_PERMANENT;
}

static void set_sold2 (FttCell * cell, GfsMovingSimulation * sim)
{
  GfsVariable * old_solid_v = sim->old_solid;
  GfsVariable ** sold2 = sim->sold2;
  FttDirection d;

  if (OLD_SOLID (cell))
    for (d = 0; d < FTT_NEIGHBORS; d++)
      SOLD2 (cell, d) = OLD_SOLID (cell)->s[d];
  else
    for (d = 0; d < FTT_NEIGHBORS; d++)
      SOLD2 (cell, d) = 1.;
}

static void check_face (FttCellFace * f, guint * nf)
{
  GfsSolidVector * s = GFS_STATE (f->cell)->solid;

  if (s && !f->neighbor && s->s[f->d] > 0. && s->s[f->d] < 1.)
    (*nf)++;
}

static void check_solid_fractions (GfsBox * box, guint * nf)
{
  FttDirection d;

  gfs_cell_check_solid_fractions (box->root);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    ftt_face_traverse_boundary (box->root, d, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttFaceTraverseFunc) check_face, nf);
}

static void is_diffusion (GfsSource * s, gboolean * diffusion)
{
  *diffusion = (GFS_IS_SOURCE_DIFFUSION (s) != NULL);
}

static void set_permanent (FttCell * cell)
{
  cell->flags |= GFS_FLAG_PERMANENT;
}

static void redistribute_old_face_in_merged (FttCell * cell, 
					     FttCell * merged, FttDirection d, 
					     GfsVariable * old_solid_v)
{  
  gint i;
  gdouble sc, sm;
  
  g_assert (cell != NULL);
  g_assert (merged != NULL);

  sc = ftt_cell_volume(cell);
  sm = ftt_cell_volume(merged);
  
  if (sc != sm)
    printf("Face redistribution not implemented yet for adaptive grid \n");
  g_assert (sc == sm);
  
  for (i = 0; i < FTT_DIMENSION;i++)
    if (i != d/2) {
      FttCellNeighbors neighbors;
      FttVector pos;
      
      ftt_cell_pos(cell,&pos);
      ftt_cell_neighbors (merged,&neighbors);

      GfsSolidVector * old_solid_merged = OLD_SOLID (merged);
      if (!old_solid_merged) {
	FttDirection c;
	OLD_SOLID (merged) = old_solid_merged = g_malloc0 (sizeof (GfsSolidVector));
	old_solid_merged->a = 1.;
	for (c = 0; c < FTT_NEIGHBORS; c++)
	  old_solid_merged->s[c] = 1.;
      }
          
      old_solid_merged->s[2*i] += OLD_SOLID (cell)->s[2*i];

      if (neighbors.c[2*i]) {
	GfsSolidVector * old_solid = OLD_SOLID (neighbors.c[2*i]);
	if (!old_solid) {
	  FttDirection c;
	  OLD_SOLID (neighbors.c[2*i]) = old_solid = g_malloc0 (sizeof (GfsSolidVector));
	  old_solid->a = 1.;
	  for (c = 0; c < FTT_NEIGHBORS; c++)
	    old_solid->s[c] = 1.;
	}
	old_solid->s[ftt_opposite_direction[2*i]] += OLD_SOLID (cell)->s[2*i];
      }
      
      old_solid_merged->s[2*i+1] += OLD_SOLID (cell)->s[2*i+1];
      
      if (neighbors.c[2*i+1]) {
	GfsSolidVector * old_solid = OLD_SOLID (neighbors.c[2*i+1]);
	if (!old_solid) {
	  FttDirection c;
	  OLD_SOLID (neighbors.c[2*i+1]) = old_solid = g_malloc0 (sizeof (GfsSolidVector));
	  old_solid->a = 1.;
	  for (c = 0; c < FTT_NEIGHBORS; c++)
	    old_solid->s[c] = 1.;
	}
	old_solid->s[ftt_opposite_direction[2*i+1]] += OLD_SOLID (cell)->s[2*i+1];	
      }
    }
}

static void redistribute_old_face (FttCell * cell, FttCell * merged, GfsVariable * old_solid) 
{
  FttCellNeighbors neighbors;
  FttDirection d;

  ftt_cell_neighbors (cell,&neighbors);
  for (d = 0; d< FTT_NEIGHBORS; d++)
    if (neighbors.c[d])
      redistribute_old_face_in_merged (cell, neighbors.c[d], d, old_solid);
}

typedef struct {
  GfsDomain * domain;
  GfsVariable * status;
  GfsVariable ** v;
} ReInitParams;

static void redistribute_destroyed_cells_content (FttCell * cell, ReInitParams * p)
{
  if (GFS_VALUE (cell, p->status) != 1.)
    return;

  GfsDomain * domain = p->domain;
  GfsVariable * old_solid_v = GFS_MOVING_SIMULATION (domain)->old_solid;
  GSList * i;
  FttCell * merged, * next;
  gdouble s1, s2;
  gint c;

  if (!OLD_SOLID (cell) || !(merged = OLD_SOLID (cell)->merged))
    return;
  while (OLD_SOLID (merged) && (next = OLD_SOLID (merged)->merged))
    merged = next;

  s1 = ftt_cell_volume (cell);
  s2 = ftt_cell_volume (merged);

  /* redistribution of the velocity */
  for (c = 0; c < FTT_DIMENSION; c++) {
    gdouble a = OLD_SOLID (merged) ? OLD_SOLID (merged)->a : 1.;
    GfsVariable * var = p->v[c];
    GFS_VALUE (merged, var) = (s1*OLD_SOLID (cell)->a*GFS_VALUE (cell, var) +
			       s2*a*GFS_VALUE (merged, var))
      /(s1*OLD_SOLID (cell)->a + s2*a);
  }
  
  /* redistribution of tracers */
  i = domain->variables;
  while (i) {
    if (GFS_IS_VARIABLE_TRACER (i->data)) {
      gdouble a = OLD_SOLID (merged) ? OLD_SOLID (merged)->a : 1.;
      GfsVariableTracer * t = GFS_VARIABLE_TRACER(i->data);
      GfsVariable * var = t->advection.v;
      GFS_VALUE (merged, var) = (s1*OLD_SOLID (cell)->a*GFS_VALUE (cell, var) +
				 s2*a*GFS_VALUE (merged, var))
	/(s1*OLD_SOLID (cell)->a + s2*a);
    }
    i = i->next;
  }
    
  if (!OLD_SOLID (merged)) {
    OLD_SOLID (merged) = g_malloc0 (sizeof (GfsSolidVector));
    OLD_SOLID (merged)->a = 1.; 
  }
  OLD_SOLID (merged)->a += s1/s2*OLD_SOLID (cell)->a;
  redistribute_old_face (cell, merged, GFS_MOVING_SIMULATION (domain)->old_solid);
}

/**
 * gfs_domain_reinit_solid_fractions:
 * @domain: a #GfsDomain.
 * @i: a list of #GfsSolids.
 *
 * Reinitializes the solid fractions of all the cells of @domain.
 *
 * If @destroy_solid is set to %TRUE, the cells entirely contained in
 * the solid are destroyed using @cleanup as cleanup function.  
 *
 * Destroy the cells that are not fluid cells anymore when the solid 
 * have moved.
 *
 * The fluid fractions of the destroyed is redistributed.
 *
 * Returns: the number of thin cells.
 */
static guint domain_reinit_solid_fractions (GfsSimulation * sim,
					    GSList * i)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsVariable * status;

  g_return_val_if_fail (sim != NULL, 0);

  status = gfs_temporary_variable (domain);
  guint thin = gfs_init_solid_fractions_leaves (domain, i, status);

  if (sim->time.t != 0.) {
    ReInitParams rp;
    rp.domain = domain;
    rp.status = status;
    rp.v = gfs_domain_velocity (domain);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) redistribute_destroyed_cells_content, &rp);
  }

  gfs_init_solid_fractions_from_children (domain, TRUE,
					  (FttCellCleanupFunc) gfs_cell_cleanup, domain, 
					  status);
  gts_object_destroy (GTS_OBJECT (status));
  return thin;
}

/**
 * reinit_solid_fractions:
 * @sim: a #GfsSimulation.
 *
 * Calls the domain_reinit_solid_fractions(). Matches the
 * boundaries by calling gfs_domain_match().
 */
static void reinit_solid_fractions (GfsSimulation * sim)
{
  guint nf = 0;
  GfsDomain * domain = GFS_DOMAIN (sim);;
  GSList * solids = gfs_simulation_get_solids (sim);
  if (solids) {
    gfs_domain_timer_start (domain, "solid_fractions");
    sim->thin = domain_reinit_solid_fractions (sim, solids);
    g_slist_free (solids);
    gfs_domain_match (domain);
    gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			       (FttCellTraverseFunc) set_permanent, NULL);
    gfs_domain_timer_stop (domain, "solid_fractions");
  }
  gts_container_foreach (GTS_CONTAINER (sim), (GtsFunc) check_solid_fractions, &nf);
  if (nf > 0) {
    GSList * i = domain->variables;
    gboolean diffusion = FALSE;
    
    while (i && !diffusion) {
      GfsVariable * v = i->data;

      if (v->sources)
	gts_container_foreach (v->sources, (GtsFunc) is_diffusion, &diffusion);
      i = i->next;
    }
    if (diffusion)
      g_warning ("the solid surface cuts %d boundary cells,\n"
		 "this may cause errors for diffusion terms\n", nf);
  }
}

/* see gfs_advection_update() for a description of what this function does */
static void moving_advection_update (GSList * merged, const GfsAdvectionParams * par)
{
  GfsVariable * old_solid_v = GFS_MOVING_SIMULATION (par->v->domain)->old_solid;

  if (merged->next == NULL) { /* cell is not merged */
    FttCell * cell = merged->data;
    gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
    gdouble olda = OLD_SOLID (cell) ? OLD_SOLID (cell)->a : 1.;

    if (GFS_IS_MIXED (cell))
      g_assert (!gfs_cell_is_small (cell));

    GFS_VALUE (cell, par->v) = (olda*GFS_VALUE (cell, par->v) + GFS_VALUE (cell, par->fv))/a;
  }
  else if (1 /* par->average */) {
    /* average value */
    GSList * i = merged;
    gdouble w = 0., total_vol = 0.;

    while (i) {
      FttCell * cell = i->data;
      gdouble vol = ftt_cell_volume (cell);
      gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
      gdouble olda = OLD_SOLID (cell) ? OLD_SOLID (cell)->a : 1.;
      
      total_vol += vol*a;
      w += vol*(olda*GFS_VALUE (cell, par->v) + GFS_VALUE (cell, par->fv));
      i = i->next;
    }
    w /= total_vol;

    i = merged;
    while (i) {
      FttCell * cell = i->data;
      GFS_VALUE (cell, par->v) = w;
      i = i->next;
    }
  }
  else {
    GSList * i = merged;
    gdouble w = 0., total_vol = 0.;

    while (i) {
      FttCell * cell = i->data;
      gdouble vol = ftt_cell_volume (cell);
      gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
      gdouble olda = OLD_SOLID (cell) ? OLD_SOLID (cell)->a : 1.;

      total_vol += vol*a;
      if (a < GFS_SMALL) {
	GFS_VALUE (cell, par->v) = olda*GFS_VALUE (cell, par->v)/a + 
	  GFS_VALUE (cell, par->fv)/GFS_SMALL;
	w += vol*GFS_VALUE (cell, par->fv)*(1. - a/GFS_SMALL);   
      }
      else
	GFS_VALUE (cell, par->v) = (olda*GFS_VALUE (cell, par->v) + GFS_VALUE (cell, par->fv))/a;
      i = i->next;
    }
    w /= total_vol;
    
    i = merged;
    while (i) {
      FttCell * cell = i->data;
      /* fixme: small cells should be excluded here?? 
	 (with corresponding modification in total_vol) */
      GFS_VALUE (cell, par->v) += w;
      i = i->next;
    }
  }
}

static double face_fraction_half (const FttCellFace * face, const GfsAdvectionParams * par)
{
  GfsVariable * old_solid_v = GFS_MOVING_SIMULATION (par->v->domain)->old_solid;
  if (face->cell && OLD_SOLID (face->cell))
    return OLD_SOLID (face->cell)->s[face->d];
  return 1.;
}

/* see gfs_face_advection_flux() for the initial implementation with static boundaries */
static void moving_face_advection_flux (const FttCellFace * face,
					const GfsAdvectionParams * par)
{
  gdouble flux;
  
  /* fixme: what's up with face mapping? */
  flux = face_fraction_half (face, par)*GFS_FACE_NORMAL_VELOCITY (face)*par->dt*
    gfs_face_upwinded_value (face, GFS_FACE_UPWINDING, NULL)/ftt_cell_size (face->cell);
  if (!FTT_FACE_DIRECT (face))
    flux = - flux;
  GFS_VARIABLE (face->cell, par->fv->i) -= flux;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_VARIABLE (face->neighbor, par->fv->i) += flux;
    break;
  case FTT_FINE_COARSE:
    GFS_VARIABLE (face->neighbor, par->fv->i) += flux/FTT_CELLS;
    break;
  default:
    g_assert_not_reached ();
  }
}

/* see gfs_face_velocity_advection_flux() for the initial implementation with static boundaries */
static void moving_face_velocity_advection_flux (const FttCellFace * face,
						 const GfsAdvectionParams * par)
{
  gdouble flux;
  FttComponent c = par->v->component;

  g_return_if_fail (c >= 0 && c < FTT_DIMENSION);

  /* fixme: what's up with face mapping? */
  flux = face_fraction_half (face, par)*GFS_FACE_NORMAL_VELOCITY (face)*
    par->dt/ftt_cell_size (face->cell);
#if 0
  if (c == face->d/2) /* normal component */
    flux *= GFS_FACE_NORMAL_VELOCITY (face);
  else /* tangential component */
#else
    flux *= gfs_face_upwinded_value (face, par->upwinding, par->u)
      /* pressure correction */
      - gfs_face_interpolated_value (face, par->g[c]->i)*par->dt/2.;
#endif
  if (!FTT_FACE_DIRECT (face))
    flux = - flux;
  GFS_VARIABLE (face->cell, par->fv->i) -= flux;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_VARIABLE (face->neighbor, par->fv->i) += flux;
    break;
  case FTT_FINE_COARSE:
    GFS_VARIABLE (face->neighbor, par->fv->i) += flux/FTT_CELLS;
    break;
  default:
    g_assert_not_reached ();
  }
}

static void moving_init (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN(sim);
  GSList * i = domain->variables;

  if (sim->advection_params.moving_order == 2)
    sim->advection_params.flux = moving_face_velocity_advection_flux;
  else
    sim->advection_params.flux = gfs_face_velocity_advection_flux;
  sim->advection_params.update = (GfsMergedTraverseFunc) moving_advection_update;

  while (i) {
    if (GFS_IS_VARIABLE_TRACER_VOF (i->data))
      g_assert_not_implemented ();
    else if (GFS_IS_VARIABLE_TRACER (i->data)) {
      GfsAdvectionParams * par = &GFS_VARIABLE_TRACER (i->data)->advection;
      if (sim->advection_params.moving_order == 2)
	par->flux = moving_face_advection_flux;
      else
	par->flux = gfs_face_advection_flux;
      par->update = sim->advection_params.update;
      par->moving_order = sim->advection_params.moving_order;
    }
    i = i->next;
  }
}

static gboolean solid_moving_event (GfsEvent * event, GfsSimulation * sim)
{
  return (GFS_SOLID_MOVING (event)->active = 
	  (* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_solid_moving_class ())->parent_class)->event) 
	  (event, sim));
}

static void solid_moving_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = solid_moving_destroy;
  GTS_OBJECT_CLASS (klass)->read = solid_moving_read;
  GTS_OBJECT_CLASS (klass)->write = solid_moving_write;
  klass->event = solid_moving_event;
}

static void solid_moving_init (GfsSolidMoving * solid)
{
  gfs_event_set (GFS_EVENT (solid),
		 0., G_MAXDOUBLE/2., -1.,
		 0, G_MAXINT/2, 1);
  solid->level = gfs_function_new (gfs_function_class (), 0.);
}

GfsEventClass * gfs_solid_moving_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo solid_moving_info = {
      "GfsSolidMoving",
      sizeof (GfsSolidMoving),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) solid_moving_class_init,
      (GtsObjectInitFunc) solid_moving_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_solid_class ()), &solid_moving_info);
  }

  return klass;
}

/* GfsMovingRun: Object */

#define MOVING_CFL 0.45

static void set_dtmax (FttCell * cell, SolidInfo * p)
{
  gdouble size = ftt_cell_size (cell);
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++) {
    gdouble v = fabs (surface_value (cell, p->v[c], NULL));
    if (v != 0.) {
      gdouble dt = size*MOVING_CFL/v;
      if (dt < p->sim->time.dtmax)
	p->sim->time.dtmax = dt;
    }
  }
}

static void moving_simulation_set_timestep (GfsSimulation * sim)
{
  gdouble dtmax = sim->time.dtmax;
  SolidInfo p;
  p.sim = sim;
  p.v = gfs_domain_velocity (GFS_DOMAIN (sim));
  gfs_domain_traverse_mixed (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseFunc) set_dtmax, &p);
  gfs_simulation_set_timestep (sim);
  sim->time.dtmax = dtmax;
}

static void swap_fractions (FttCell * cell, GfsVariable * old_solid_v) {
  FttDirection c;

  g_assert (cell);

  if (FTT_CELL_IS_LEAF(cell)) {
    if (OLD_SOLID (cell)) {
      GfsSolidVector * solid_old = OLD_SOLID (cell);
      
      if (GFS_STATE (cell)->solid) {
	GfsSolidVector * solid = GFS_STATE (cell)->solid;
	
	for (c = 0; c < 2*FTT_DIMENSION; c++) 
	  if (solid->s[c] == 0.)
	    solid_old->s[c] = 0;
	  else
	    solid_old->s[c] = (solid_old->s[c]+solid->s[c])/2. ;	
      }
      else
	for (c = 0; c < 2*FTT_DIMENSION; c++)
	  solid_old->s[c] = (solid_old->s[c]+1.)/2. ;
    }
    else if (GFS_STATE (cell)->solid) {
      GfsSolidVector * solid = GFS_STATE (cell)->solid;
      GfsSolidVector * solid_old = OLD_SOLID (cell) = g_malloc0 (sizeof (GfsSolidVector));
      OLD_SOLID (cell)->a= 1.;

      for (c = 0; c < 2*FTT_DIMENSION; c++) 
	solid_old->s[c] = 1.;	

      for (c = 0; c < 2*FTT_DIMENSION; c++) 
	if (solid->s[c] == 0.)
	  solid_old->s[c] = 0;
	else
	  solid_old->s[c] = (solid_old->s[c]+solid->s[c])/2. ;
    }
  }
  
  if (OLD_SOLID (cell)) {
    if (GFS_STATE(cell)->solid) {
      GfsSolidVector * tmp = OLD_SOLID (cell);
      OLD_SOLID (cell) = GFS_STATE(cell)->solid;
      GFS_STATE(cell)->solid = tmp;
      tmp = NULL;
    }
    else {
      GFS_STATE(cell)->solid = OLD_SOLID (cell);
      OLD_SOLID (cell) = NULL;
    }
  }
  else if (GFS_STATE(cell)->solid) {
    OLD_SOLID (cell) = GFS_STATE(cell)->solid;
    GFS_STATE(cell)->solid = NULL;
  }
}

static void old_solid_fractions_from_children (FttCell * cell)
{
  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren child;
    guint i;
    
    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i])
	old_solid_fractions_from_children (child.c[i]);
    
      gfs_cell_init_solid_fractions_from_children (cell);
  }
}

static void foreach_box (GfsBox * box, gpointer data)
{
  old_solid_fractions_from_children (box->root);
}

static void swap_face_fractions (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) swap_fractions, 
			    GFS_MOVING_SIMULATION (sim)->old_solid);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) foreach_box, NULL);
}

static void swap_fractions_back (FttCell * cell, GfsVariable * old_solid_v) 
{
  if (OLD_SOLID (cell))
    if (GFS_STATE(cell)->solid) {
      GfsSolidVector * tmp = OLD_SOLID (cell);
      OLD_SOLID (cell) = GFS_STATE(cell)->solid;
      GFS_STATE(cell)->solid = tmp;
      tmp = NULL;
    }
    else {
      GFS_STATE(cell)->solid = OLD_SOLID (cell);
      OLD_SOLID (cell) = NULL;
    }
  else
    if (GFS_STATE(cell)->solid) {
      OLD_SOLID (cell) = GFS_STATE(cell)->solid;
      GFS_STATE(cell)->solid = NULL;
    }
}

static void swap_face_fractions_back (GfsSimulation * sim) 
{
  gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) swap_fractions_back,
			    GFS_MOVING_SIMULATION (sim)->old_solid);
}

static void move_vertex (GtsPoint * p, SolidInfo * par)
{ 
  FttVector pos = *((FttVector *) &p->x);
  FttCell * cell = gfs_domain_locate (GFS_DOMAIN (par->sim), pos, -2);
  if (cell) {
    gdouble dt = par->sim->advection_params.dt;
    FttComponent c;
    for (c = 0; c < FTT_DIMENSION; c++)
      (&p->x)[c] += surface_value (cell, par->v[c], &pos)*dt;
  }
  else
    g_warning ("point %g,%g not in domain", pos.x, pos.y);
}

static void solid_move_remesh (GfsSolidMoving * solid, GfsSimulation * sim)
{
  GfsSurface * surface = GFS_SURFACE (GFS_SOLID (solid)->s);
  if (surface->s) {
    SolidInfo p;
    p.sim = sim;
    p.s = solid;
    p.v = gfs_domain_velocity (GFS_DOMAIN (sim));
    gts_surface_foreach_vertex (surface->s, (GtsFunc) move_vertex, &p);
  }
  else
    /* implicit surface */
    g_assert_not_implemented ();
  remesh_surface_moving (sim, solid);
}

static void move_solids (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsVariable * old_solid = GFS_MOVING_SIMULATION (sim)->old_solid;
  GfsVariable * sold2[FTT_NEIGHBORS];

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) set_old_solid, old_solid);

  if (sim->advection_params.moving_order == 2) {
    FttDirection d;
    for (d = 0; d < FTT_NEIGHBORS; d++) {
      sold2[d] = gfs_domain_add_variable (domain, NULL, NULL);
      sold2[d]->coarse_fine = sold2_fine_init;
    }
    GFS_MOVING_SIMULATION (sim)->sold2 = sold2;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) set_sold2, sim);
  }

  GSList * solids = gfs_simulation_get_solids (sim), * s = solids;
  while (s) {
    if (GFS_IS_SOLID_MOVING (s->data) && GFS_SOLID_MOVING (s->data)->active)
      solid_move_remesh (s->data, sim);
    s = s->next;
  }
  g_slist_free (solids);
  reinit_solid_fractions (sim);
  gfs_set_merged (domain);

  if (sim->advection_params.moving_order == 2) {
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) second_order_face_fractions, sim);
    FttDirection d;
    for (d = 0; d < FTT_NEIGHBORS; d++)
      gts_object_destroy (GTS_OBJECT (sold2[d]));    
    GFS_MOVING_SIMULATION (sim)->sold2 = NULL;
  }
}

typedef struct {
  GfsDomain * domain;
  gdouble dt;
  FttComponent c;
  GfsVariable * div;
  GfsVariable * v;
} DivergenceData;

static void moving_divergence_approx (FttCell * cell, DivergenceData * p)
{
  GFS_VALUE (cell, p->div) += 
    GFS_STATE (cell)->solid->fv*(GFS_STATE (cell)->solid->s[2*p->c + 1] -
				 GFS_STATE (cell)->solid->s[2*p->c])*ftt_cell_size (cell);
}

static void moving_approximate_projection (GfsDomain * domain,
					   GfsMultilevelParams * par,
					   GfsAdvectionParams * apar,
					   GfsVariable * p,
					   GfsFunction * alpha,
					   GfsVariable * res,
					   GfsVariable ** g)
{
  DivergenceData q;
  GfsVariable ** v = gfs_domain_velocity (domain);
  
  q.div = gfs_temporary_variable (domain);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, q.div);
  for (q.c = 0; q.c < FTT_DIMENSION; q.c++) {
    gfs_domain_surface_bc (domain, v[q.c]);
    gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			       (FttCellTraverseFunc) moving_divergence_approx, &q);
  }
  gfs_approximate_projection (domain, par, apar, p, alpha, q.div, res, g);
  gts_object_destroy (GTS_OBJECT (q.div));
}

static void moving_divergence_mac (FttCell * cell, DivergenceData * p)
{
  GfsVariable * old_solid_v = GFS_MOVING_SIMULATION (p->domain)->old_solid;
  gdouble size = ftt_cell_size (cell);
  gdouble a = GFS_STATE (cell)->solid ? GFS_STATE (cell)->solid->a : 1.;
  gdouble olda = OLD_SOLID (cell) ? OLD_SOLID (cell)->a : 1.;
  
  GFS_VALUE (cell, p->div) = (olda - a)*size*size/p->dt;
}

static void moving_mac_projection (GfsSimulation * sim,
				   GfsMultilevelParams * par,
				   GfsAdvectionParams * apar,
				   GfsVariable * p,
				   GfsFunction * alpha,
				   GfsVariable ** g)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  DivergenceData q;

  if (apar->moving_order == 2) {
    q.dt = apar->dt;
    swap_face_fractions (sim);
  }
  else /* first order */
    q.dt = - apar->dt;

  q.div = gfs_temporary_variable (domain);
  q.domain = domain;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) moving_divergence_mac, &q);
  gfs_mac_projection (domain, par, apar, p, alpha, q.div, g);
  gts_object_destroy (GTS_OBJECT (q.div));

  if (apar->moving_order == 2)
    swap_face_fractions_back (sim);
}

static void moving_simulation_run (GfsSimulation * sim)
{
  GfsVariable * p, * pmac, * res = NULL, * g[FTT_DIMENSION], * gmac[FTT_DIMENSION];
  GfsVariable ** gc = sim->advection_params.gc ? g : NULL;
  GfsDomain * domain;
  GSList * i;

  domain = GFS_DOMAIN (sim);

  p = gfs_variable_from_name (domain->variables, "P");
  g_assert (p);
  pmac = gfs_variable_from_name (domain->variables, "Pmac");
  g_assert (pmac);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++) {
    gmac[c] = gfs_temporary_variable (domain);
    gfs_variable_set_vector (gmac[c], c);
    if (sim->advection_params.gc) {
      g[c] = gfs_temporary_variable (domain);
      gfs_variable_set_vector (g[c], c);
    }
    else
      g[c] = gmac[c];
  }

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  i = domain->variables;
  while (i) {
    if (GFS_IS_VARIABLE_RESIDUAL (i->data))
      res = i->data;
    i = i->next;
  }

  moving_init (sim);

  moving_simulation_set_timestep (sim);
  if (sim->time.i == 0)
    moving_approximate_projection (domain,
				   &sim->approx_projection_params,
				   &sim->advection_params,
				   p, sim->physical_params.alpha, res, g);
  else if (sim->advection_params.gc)
    gfs_update_gradients (domain, p, sim->physical_params.alpha, g);

  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    move_solids (sim);
   
    gfs_predicted_face_velocities (domain, FTT_DIMENSION, &sim->advection_params);
    
    gfs_variables_swap (p, pmac);

    moving_mac_projection (sim,
			   &sim->projection_params, 
			   &sim->advection_params,
			   p, sim->physical_params.alpha, gmac);

    gfs_variables_swap (p, pmac);
    
    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_half_do, sim);
   
    gfs_centered_velocity_advection_diffusion (domain,
					       FTT_DIMENSION,
					       &sim->advection_params,
					       gmac,
					       sim->time.i > 0 || !gc ? gc : gmac,
					       sim->physical_params.alpha);
        
    gfs_advance_tracers (domain, sim->advection_params.dt);

    if (gc) {
      gfs_source_coriolis_implicit (domain, sim->advection_params.dt);
      gfs_correct_centered_velocities (domain, FTT_DIMENSION, sim->time.i > 0 ? gc : gmac, 
				       -sim->advection_params.dt);
    }
    else if (gfs_has_source_coriolis (domain)) {
      gfs_correct_centered_velocities (domain, FTT_DIMENSION, gmac, sim->advection_params.dt);
      gfs_source_coriolis_implicit (domain, sim->advection_params.dt);
      gfs_correct_centered_velocities (domain, FTT_DIMENSION, gmac, -sim->advection_params.dt);
    }

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gfs_simulation_adapt (sim);

    moving_approximate_projection (domain,
				   &sim->approx_projection_params, 
				   &sim->advection_params, p, sim->physical_params.alpha, res, g);

    sim->time.t = sim->tnext;
    sim->time.i++;

    moving_simulation_set_timestep (sim);

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);

  for (c = 0; c < FTT_DIMENSION; c++) {
    gts_object_destroy (GTS_OBJECT (gmac[c]));
    if (sim->advection_params.gc)
      gts_object_destroy (GTS_OBJECT (g[c]));
  }
}

static void moving_simulation_class_init (GfsSimulationClass * klass)
{
  klass->run = moving_simulation_run;
}

static void old_solid_cleanup (FttCell * cell, GfsVariable * old_solid_v)
{
  g_free (OLD_SOLID (cell));
  OLD_SOLID (cell) = NULL;
}

static void none (void) {}

static void moving_simulation_init (GfsDomain * domain)
{
  gfs_domain_add_variable (domain, "div", "Divergence")->centered = TRUE;

  /* old_solid will hold a pointer to a GfsSolidVector */
  GfsVariable * old_solid = gfs_domain_add_variable (domain, NULL, NULL); 
  GFS_MOVING_SIMULATION (domain)->old_solid = old_solid;
  /* pointers need to be "interpolated" correctly (i.e. not at all) */
  old_solid->coarse_fine = (GfsVariableFineCoarseFunc) none;
  old_solid->fine_coarse = (GfsVariableFineCoarseFunc) none;
  /* the memory needs to be freed when the cell is cleaned up */
  old_solid->cleanup = (FttCellCleanupFunc) old_solid_cleanup;
  /* switch off boundary conditions */
  GfsBc * bc = gfs_bc_new (gfs_bc_class (), old_solid, FALSE);
  bc->bc = bc->homogeneous_bc = bc->face_bc = (FttFaceTraverseFunc) none;
  gfs_variable_set_default_bc (old_solid, bc);
}

GfsSimulationClass * gfs_moving_simulation_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_moving_simulation_info = {
      "GfsMovingSimulation",
      sizeof (GfsMovingSimulation),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) moving_simulation_class_init,
      (GtsObjectInitFunc) moving_simulation_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()), 
				  &gfs_moving_simulation_info);
  }

  return klass;
}
