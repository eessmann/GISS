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

#ifndef __SOLID_H__
#define __SOLID_H__

#include <gts.h>

#include "fluid.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void         gfs_cell_fluid                              (FttCell * cell);
void         gfs_cell_init_solid_fractions        (FttCell * root, 
						   GtsSurface * s,
						   GNode * stree,
						   gboolean is_open,
						   gboolean destroy_solid,
						   FttCellCleanupFunc cleanup,
						   gpointer data);
void         gfs_cell_init_solid_fractions_from_children (FttCell * cell);
gboolean     gfs_cell_check_solid_fractions              (FttCell * root,
							  GtsSurface * solid,
							  GNode * stree,
							  gboolean is_open);
gboolean     gfs_refine_mixed                       (const FttCell * cell);
void         gfs_cell_init_fraction                 (FttCell * root, 
						     GtsSurface * s,
						     GNode * stree,
						     gboolean is_open,
						     GfsVariable * c);
void         gfs_cell_cm                            (const FttCell * cell, 
						     FttVector * cm);
void         gfs_face_ca                            (const FttCellFace * face, 
						     FttVector * ca);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __SOLID_H__ */
