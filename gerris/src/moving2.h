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

/*
 * Code only used by second-order moving solid boundaries
 */

#ifndef GISS_MOVING2_H
#define GISS_MOVING2_H

#include "moving.h"

#define SOLD2(c, d)  (GFS_VALUE (c, sold2[d]))

void sold2_fine_init(FttCell *parent, GfsVariable *v);

int cell_is_corner(FttCell *cell);

int cell_was_corner(FttCell *cell, GfsVariable *old_solid_v, GfsVariable **sold2);

double new_fluid_old_solid(FttCell *cell, FttDirection d1,
                                  GfsVariable *old_solid,
                                  GfsVariable **sold2);

double new_solid_old_fluid(FttCell *cell, FttDirection d1,
                                  GfsVariable *old_solid,
                                  GfsVariable **sold2);

double new_solid_old_solid(FttCell *cell, FttDirection d1,
                                  GfsVariable *old_solid,
                                  GfsVariable **sold2);

void second_order_face_fractions(FttCell *cell, GfsSimulationMoving *sim);

void set_sold2(FttCell *cell, GfsSimulationMoving *sim);

void redistribute_old_face_in_merged(FttCell *cell,
                                            FttCell *merged, FttDirection d,
                                            GfsVariable *old_solid_v);

void redistribute_old_face(FttCell *cell, FttCell *merged, GfsVariable *old_solid);

double face_fraction_half(const FttCellFace *face, const GfsAdvectionParams *par);

void moving_face_advection_flux(const FttCellFace *face,
                                       const GfsAdvectionParams *par);

void moving_face_velocity_advection_flux(const FttCellFace *face,
                                                const GfsAdvectionParams *par);

void swap_fractions(FttCell *cell, GfsVariable *old_solid_v);

void old_solid_fractions_from_children(FttCell *cell);

void foreach_box(GfsBox *box, gpointer data);

void swap_face_fractions(GfsSimulation *sim);

void swap_fractions_back(FttCell *cell, GfsVariable *old_solid_v);

void swap_face_fractions_back(GfsSimulation *sim);

void moving_divergence_distribution_second_order(GSList *merged, DivergenceData *p);



#endif //GISS_MOVING2_H
