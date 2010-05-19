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

#ifndef __INIT_H__
#define __INIT_H__

#include <gts.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

GtsObjectClass ** gfs_classes             (void);
void              gfs_init                (int * argc, 
					   char *** argv);
void gfs_catch_floating_point_exceptions   (void);
int  gfs_restore_floating_point_exceptions (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __INIT_H__ */
