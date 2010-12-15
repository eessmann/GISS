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

/* GfsSourceMass: Header */

typedef struct _GfsSourceMass         GfsSourceMass;

struct _GfsSourceMass {
  /*< private >*/
  GfsSourceVelocity parent;

  /*< public >*/
  GfsFunction * intensity;
};

typedef struct {                                                              
  GfsVariable * v, * div;                                                      
} SourceMassPar;   

#define GFS_SOURCE_MASS(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceMass,\
					         gfs_source_mass_class ())
#define GFS_IS_SOURCE_MASS(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_mass_class ()))

GfsSourceGenericClass * gfs_source_mass_class  (void);

/* GfsSourceMass: Object */

static void gfs_source_mass_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_SOURCE_MASS (o)->intensity));

  (* GTS_OBJECT_CLASS (gfs_source_mass_class ())->parent_class->destroy) (o) ;
}

static void source_divergence_mac (FttCell * cell, SourceMassPar * p)
{
  gdouble sum = 0.;                                                            
  GSList * i = GTS_SLIST_CONTAINER (p->v->sources)->items;                    
                                                                              
  while (i) {                                                                 
    GfsSourceGeneric * s = i->data;                                           
                                                                              
    if (s->centered_value)                                                    
      sum += (* s->centered_value) (s, cell, p->v);                           
    i = i->next;                                                              
  }                                                                           

  GFS_VALUE (cell, p->div) = sum*ftt_cell_volume (cell);
}


static void divergence_source (GfsDomain * domain, GfsAdvectionParams * apar, GfsVariable * div) 
{
  SourceMassPar p;
  p.div = div;
  p.v = gfs_variable_from_name (domain->variables, "P");
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                            (FttCellTraverseFunc) source_divergence_mac, &p); 
}

static gdouble source_mass_value (GfsSourceGeneric * s,
                                  FttCell * cell,
                                  GfsVariable * v)
{
  return gfs_function_value (GFS_SOURCE_MASS (s)->intensity, cell);
}

static void gfs_source_mass_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_mass_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  
  GfsSourceMass * s = GFS_SOURCE_MASS (*o);
  GfsSimulation * sim = gfs_object_simulation (s);
  
  sim->divergence_hook = divergence_source ;

  gfs_function_read (s->intensity, gfs_object_simulation (*o), fp);
  if (fp->type != GTS_ERROR) {
    GfsSourceGeneric * s = GFS_SOURCE_GENERIC (*o);
    s->mac_value = s->centered_value = source_mass_value;
    s->face_value = NULL;
  }

  GfsVariable * v = gfs_variable_from_name (GFS_DOMAIN(sim)->variables, "P");
  if (v->sources == NULL)
    v->sources = gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
  gts_container_add (v->sources, GTS_CONTAINEE (s));
  
}

static void gfs_source_mass_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_source_mass_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_source_mass_class ())->parent_class->write) 
      (o, fp); 
  GfsSourceMass * s = GFS_SOURCE_MASS (o);
  gfs_function_write (s->intensity, fp); 
}

static void gfs_source_mass_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_source_mass_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_source_mass_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_source_mass_destroy;
}

static void gfs_source_mass_init ( GfsSourceMass *s )
{
  s->intensity = gfs_function_new (gfs_function_class (), 0.);  
 // gfs_function_set_units (s->intensity, 0.); //check units (s^(-1)) Is it really required? Is it not divergence automatically rescaled?
}

GfsSourceGenericClass * gfs_source_mass_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_mass_info = {
      "GfsSourceMass",
      sizeof (GfsSourceMass),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_mass_class_init,
      (GtsObjectInitFunc) gfs_source_mass_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()),
				  &gfs_source_mass_info);
  }

  return klass;
}


/* GfsPointSourceMass: Header */

typedef struct _GfsPointSourceMass         GfsPointSourceMass;

struct _GfsPointSourceMass {
  /*< private >*/
  GfsSourceMass parent;

  /*< public >*/
  FttVector pos;
};

#define GFS_POINT_SOURCE_MASS(obj)            GTS_OBJECT_CAST (obj,\
					         GfsPointSourceMass,\
					         gfs_point_source_mass_class ())
#define GFS_IS_POINT_SOURCE_MASS(obj)         (gts_object_is_from_class (obj,\
						 gfs_point_source_mass_class ()))

GfsSourceGenericClass * gfs_point_source_mass_class  (void);

/* GfsPointSourceMass: Object */

static void point_divergence_source (GfsDomain * domain, GfsAdvectionParams * apar, GfsVariable * div) 
{
  GfsVariable * v = gfs_variable_from_name (domain->variables, "P");

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                           (FttCellTraverseFunc) gfs_cell_reset, div); //is it required?

  GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;                    
                                                                              
  while (i) {                                                                 
    GfsSourceGeneric * s = i->data; 
    GfsSourceMass * sm = GFS_SOURCE_MASS(s);                                          
    GfsPointSourceMass * ps = GFS_POINT_SOURCE_MASS(s);                                          
    FttCell * cell = gfs_domain_locate (domain, ps->pos, -1, NULL);
    GFS_VALUE (cell, div)  += gfs_function_value (sm->intensity, cell) ;
                                                                              
    i = i->next;                                                              
  }                                                                           

}

static void gfs_point_source_mass_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_point_source_mass_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  
  GfsPointSourceMass * s = GFS_POINT_SOURCE_MASS (*o);
  GfsSimulation * sim = gfs_object_simulation (s);

  sim->divergence_hook = point_divergence_source ;

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
      gts_file_error (fp, "expecting a number (p.x)");      
      return;
  }
  s->pos.x = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
      gts_file_error (fp, "expecting a number (p.y)");
      return;
  }
  s->pos.y = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
      gts_file_error (fp, "expecting a number (p.z)");
      return;
  }
  s->pos.z = atof (fp->token->str);
  gts_file_next_token (fp);

}

static void gfs_point_source_mass_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_point_source_mass_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_point_source_mass_class ())->parent_class->write) 
      (o, fp); 
  GfsPointSourceMass * s = GFS_POINT_SOURCE_MASS (o);
  fprintf (fp, " %g %g %g", s->pos.x, s->pos.y, s->pos.z);
}

static void gfs_point_source_mass_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_point_source_mass_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_point_source_mass_write;
}

GfsSourceGenericClass * gfs_point_source_mass_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_point_source_mass_info = {
      "GfsPointSourceMass",
      sizeof (GfsPointSourceMass),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_point_source_mass_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_mass_class ()),
				  &gfs_point_source_mass_info);
  }

  return klass;
}
 
/* Initialize module */

const gchar gfs_module_name[] = "sourcemass";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
  gfs_source_mass_class ();
  gfs_point_source_mass_class ();
  return NULL; 
}
