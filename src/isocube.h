/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2004 National Institute of Water and Atmospheric
 * Research
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

/* isocube adapted from GTS (see gts/src/iso.c and gts/doc/isocube.fig) */
static guint edge1[12][2] = {
  {0, 4}, {1, 5}, {3, 7}, {2, 6},
  {0, 2}, {1, 3}, {5, 7}, {4, 6},
  {0, 1}, {4, 5}, {6, 7}, {2, 3}
};
static FttVector vertex[8] = {
  {0.,0.,0.},{0.,0.,1.},{0.,1.,0.},{0.,1.,1.},
  {1.,0.,0.},{1.,0.,1.},{1.,1.,0.},{1.,1.,1.}
};
static guint face[6][4][2] = {
  {{7,0},{10,0},{6,1},{9,1}}, /* right */
  {{4,0},{11,0},{5,1},{8,1}}, /* left */
  {{3,0},{10,0},{2,1},{11,1}},/* top */
  {{0,0},{9,0},{1,1},{8,1}},  /* bottom */
  {{1,0},{6,0},{2,1},{5,1}},  /* front */
  {{0,0},{7,0},{3,1},{4,1}}   /* back */
};
static guint face_v[6][4] = {
  {4,6,7,5},/* right */
  {0,2,3,1},/* left */
  {2,6,7,3},/* top */
  {0,4,5,1},/* bottom */
  {1,5,7,3},/* front */
  {0,4,6,2} /* back */
};
