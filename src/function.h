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

static double Dirichlet = 1.;
static double Neumann = 0.;
static GfsSimulation * _sim = NULL;
static FttCell * _cell = NULL;

static double dd (const gchar * name, FttComponent c) {
  GfsVariable * v = gfs_variable_from_name (GFS_DOMAIN (_sim)->variables, name);
  if (v == NULL)
    return 0.;
  g_return_val_if_fail (_cell != NULL, 0.);
  return gfs_center_gradient (_cell, c, v->i)/ftt_cell_size (_cell);
}

static double dx (const gchar * name) { return dd (name, FTT_X); }
static double dy (const gchar * name) { return dd (name, FTT_Y); }
#if !FTT_2D
static double dz (const gchar * name) { return dd (name, FTT_Z); }
#endif /* 3D */

#endif /* __FUNCTION_H__ */
