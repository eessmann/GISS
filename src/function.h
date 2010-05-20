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

#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#define NODATA GFS_NODATA

static double Dirichlet = 1.;
static double Neumann = 0.;
static GfsSimulation * _sim = NULL;
static FttCell * _cell = NULL;

static double dd (const gchar * name, FttComponent c) {
  GfsVariable * v = gfs_variable_from_name (GFS_DOMAIN (_sim)->variables, name);
  if (v == NULL)
    return 0.;
  g_return_val_if_fail (_cell != NULL, 0.);
  return gfs_dimensional_value (v, gfs_center_gradient (_cell, c, v->i)/
				(_sim->physical_params.L*ftt_cell_size (_cell)));
}

static double dx (const gchar * name) { return dd (name, FTT_X); }
static double dy (const gchar * name) { return dd (name, FTT_Y); }
#if !FTT_2D
static double dz (const gchar * name) { return dd (name, FTT_Z); }
#endif /* 3D */

static double area (const gchar * name)
{
  GfsVariable * v = gfs_variable_from_name (GFS_DOMAIN (_sim)->variables, name);
  if (v == NULL || !GFS_IS_VARIABLE_TRACER_VOF (v))
    return 0.;
  g_return_val_if_fail (_cell != NULL, 0.);
  GfsVariableTracerVOF * t = GFS_VARIABLE_TRACER_VOF (v);
  FttVector m, p;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&m.x)[c] = GFS_VALUE (_cell, t->m[c]);
  return gfs_plane_area_center (&m, GFS_VALUE (_cell, t->alpha), &p)/
    (_sim->physical_params.L*ftt_cell_size (_cell));
}

static double correctness (const gchar * name)
{
  GfsVariable * v = gfs_variable_from_name (GFS_DOMAIN (_sim)->variables, name);
  if (v == NULL || !GFS_IS_VARIABLE_TRACER_VOF (v))
    return 0.;
  g_return_val_if_fail (_cell != NULL, 0.);
  return gfs_vof_correctness (_cell, GFS_VARIABLE_TRACER_VOF (v));
}

static double distance (double xo, double yo, double zo)
{
  /* fixme: this doesn't take mapping into account properly */
  GtsPoint o;
  o.x = xo; o.y = yo; o.z = zo;
  gfs_simulation_map (_sim, (FttVector *) &o.x);
  GtsBBox bb;
  ftt_cell_bbox (_cell, &bb);
  gdouble min, max;
  gts_bbox_point_distance2 (&bb, &o, &min, &max);
  return sqrt (min)*_sim->physical_params.L;
}

#endif /* __FUNCTION_H__ */
