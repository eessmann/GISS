/* Gerris - The GNU Flow Solver			(-*-C-*-)
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
 * Author of the object: Gaurav Tomar <gaurav.tomar03@gmail.com>
 */
#include "particulates.h"

/* DropletToParticle */
typedef struct _GfsDropletToParticle GfsDropletToParticle;
struct _GfsDropletToParticle{
  /*< private >*/
  GfsParticleList parent;
  GfsVariable * v;
  
  /*< public >*/
  GfsFunction * fc;
  GfsVariable * c;
  gint min;
  gdouble resetwith;
  gdouble density;
};

#define DROPLET_TO_PARTICLE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsDropletToParticle,\
					         gfs_droplet_to_particle_class ())

#define IS_DROPLET_TO_PARTICLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_droplet_to_particle_class ()))

GfsEventClass * gfs_droplet_to_particle_class  (void);
