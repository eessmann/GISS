/* Gerris - The GNU Flow Solver			(-*-C-*-)
 * Copyright (C) 2010 National Institute of Water and Atmospheric Research
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

#include "simulation.h"
#include "source.h"
#include "adaptive.h"
#include "output.h"
#include "solid.h"

#include <stdlib.h>

/* GfsSourceVolume: Header */

typedef struct _GfsSourceVolume         GfsSourceVolume;

struct _GfsSourceVolume {
  /*< private >*/
  GfsSourceGeneric parent;

  /*< public >*/
  GfsFunction * intensity;
};

#define GFS_SOURCE_VOLUME(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceVolume,\
					         gfs_source_volume_class ())
#define GFS_IS_SOURCE_VOLUME(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_volume_class ()))

GfsSourceGenericClass * gfs_source_volume_class  (void);

/* GfsSourcevolume: Object */

static void gfs_source_volume_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_SOURCE_VOLUME (o)->intensity));

  (* GTS_OBJECT_CLASS (gfs_source_volume_class ())->parent_class->destroy) (o) ;
}

typedef struct {                                                              
  GfsVariable * v, * div;                                                      
} SourceVolumePar;   

static void source_divergence_mac (FttCell * cell, SourceVolumePar * p)
{
  gdouble sum = 0.;                                                            
  GSList * i = GTS_SLIST_CONTAINER (p->v->sources)->items;                    
                                                                              
  while (i) {                                                                 
    GfsSourceGeneric * s = i->data;                                           
                                                                              
    if (s->centered_value)                                                    
      sum += (* s->centered_value) (s, cell, p->v);                           
    i = i->next;                                                              
  }                                                                           

  /* div*h^2 (units source m^3/s/vol) */
  GFS_VALUE (cell, p->div) -= sum*ftt_cell_volume (cell);
}


static void divergence_source (GfsDomain * domain, GfsAdvectionParams * apar, GfsVariable * div) 
{
  SourceVolumePar p;
  p.div = div;
  p.v = gfs_variable_from_name (domain->variables, "P");
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                            (FttCellTraverseFunc) source_divergence_mac, &p); 
}

static gdouble source_volume_value (GfsSourceGeneric * s,
				    FttCell * cell,
				    GfsVariable * v)
{
  return gfs_function_value (GFS_SOURCE_VOLUME (s)->intensity, cell);
}

static void gfs_source_volume_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_volume_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  
  GfsSourceVolume * sv = GFS_SOURCE_VOLUME (*o);
  GfsSimulation * sim = gfs_object_simulation (sv);
  
  gfs_function_read (sv->intensity, sim, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsSourceGeneric * s = GFS_SOURCE_GENERIC (*o);
  s->mac_value = s->centered_value = source_volume_value;
  s->face_value = NULL;
  sim->divergence_hook = divergence_source;

  GfsVariable * v = gfs_variable_from_name (GFS_DOMAIN (sim)->variables, "P");
  g_assert (v);
  if (v->sources == NULL)
    v->sources = gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
  gts_container_add (v->sources, GTS_CONTAINEE (s));  
}

static void gfs_source_volume_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_volume_class ())->parent_class->write) (o, fp); 
  gfs_function_write (GFS_SOURCE_VOLUME (o)->intensity, fp); 
}

static void gfs_source_volume_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_source_volume_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_source_volume_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_source_volume_destroy;
}

static void gfs_source_volume_init ( GfsSourceVolume *s )
{
  s->intensity = gfs_function_new (gfs_function_class (), 0.);
}

GfsSourceGenericClass * gfs_source_volume_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_volume_info = {
      "GfsSourceVolume",
      sizeof (GfsSourceVolume),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_volume_class_init,
      (GtsObjectInitFunc) gfs_source_volume_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()),
				  &gfs_source_volume_info);
  }

  return klass;
}

/* Initialize module */

const gchar gfs_module_name[] = "sourcevolume";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
  gfs_source_volume_class ();
  return NULL; 
}
