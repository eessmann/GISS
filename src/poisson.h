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

#ifndef __POISSON_H__
#define __POISSON_H__

#include <glib.h>
#include <gts.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "domain.h"

void                  gfs_relax                      (GfsDomain * domain,
						      guint d,
						      gint max_depth,
						      GfsVariable * u,
						      GfsVariable * rhs,
						      GfsVariable * dia);
void                  gfs_residual                   (GfsDomain * domain,
						      guint d,
						      FttTraverseFlags flags,
						      gint max_depth,
						      GfsVariable * u,
						      GfsVariable * rhs,
						      GfsVariable * dia,
						      GfsVariable * res);
void                  gfs_poisson_coefficients       (GfsDomain * domain,
						      GfsVariable * c,
						      gdouble rho);
void                  gfs_poisson_cycle              (GfsDomain * domain,
						      guint d,
						      guint levelmin,
						      guint depth,
						      guint nrelax,
						      GfsVariable * u,
						      GfsVariable * rhs,
						      GfsVariable * dia,
						      GfsVariable * res);

void                  gfs_diffusion_coefficients     (GfsDomain * domain,
						      GfsSourceDiffusion * d,
						      gdouble dt,
						      GfsVariable * dia);
void                  gfs_viscosity_coefficients     (GfsDomain * domain,
						      GfsSourceDiffusion * d,
						      gdouble dt,
						      GfsVariable * c,
						      gdouble rho,
						      GfsVariable * dia);
void                  gfs_diffusion_rhs              (GfsDomain * domain,
						      GfsVariable * v,
						      GfsVariable * rhs,
						      GfsVariable * dia);
void                  gfs_diffusion_residual         (GfsDomain * domain,
						      GfsVariable * u,
						      GfsVariable * rhs,
						      GfsVariable * dia,
						      GfsVariable * res);
void                  gfs_diffusion_cycle            (GfsDomain * domain,
						      guint levelmin,
						      guint depth,
						      guint nrelax,
						      GfsVariable * u,
						      GfsVariable * rhs,
						      GfsVariable * dia,
						      GfsVariable * res);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __POISSON_H__ */
