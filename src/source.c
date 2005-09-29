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

#include <stdlib.h>
#include <math.h>
#include "source.h"
#include "simulation.h"
#include "solid.h"

/**
 * gfs_variable_mac_source:
 * @v: a #GfsVariable.
 * @cell: a #FttCell.
 *
 * Returns: the sum of all the sources for variable @v in @cell.
 */
gdouble gfs_variable_mac_source (GfsVariable * v, FttCell * cell)
{
  gdouble sum;
  GSList * i;

  g_return_val_if_fail (v != NULL, 0.);
  g_return_val_if_fail (cell != NULL, 0.);

  if (v->sources == NULL)
    return 0.;

  sum = 0.;
  i = GTS_SLIST_CONTAINER (v->sources)->items;
  while (i) {
    GtsObject * o = i->data;

    if (GFS_SOURCE_GENERIC_CLASS (o->klass)->mac_value)
      sum += (* GFS_SOURCE_GENERIC_CLASS (o->klass)->mac_value) 
	(GFS_SOURCE_GENERIC (o), cell, v);
    i = i->next;
  }
  return sum;
}

static void add_sources (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsVariable * sv = data[1];
  gdouble * dt = data[2];
  GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
  gdouble sum = 0;
  
  i = GTS_SLIST_CONTAINER (v->sources)->items;
  while (i) {
    GtsObject * o = i->data;

    if (GFS_SOURCE_GENERIC_CLASS (o->klass)->centered_value)
      sum += (* GFS_SOURCE_GENERIC_CLASS (o->klass)->centered_value)
	(GFS_SOURCE_GENERIC (o), cell, v);
    i = i->next;
  }
  GFS_VARIABLE (cell, sv->i) += (*dt)*sum;
}

/**
 * gfs_domain_variable_centered_sources:
 * @domain: a #GfsDomain.
 * @v: a #GfsVariable.
 * @sv: a #GfsVariable.
 * @dt: the timestep.
 *
 * Adds the source terms for @v to @sv.
 */
void gfs_domain_variable_centered_sources (GfsDomain * domain, 
					   GfsVariable * v,
					   GfsVariable * sv,
					   gdouble dt)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);
  g_return_if_fail (sv != NULL);

  if (v->sources) {
    gpointer data[3];
    
    data[0] = v;
    data[1] = sv;
    data[2] = &dt;
    gfs_domain_cell_traverse (domain, 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) add_sources, data);    
  }
}

/* GfsSourceGeneric: Object */

static void source_generic_init (GfsSourceGeneric * s)
{
  GFS_EVENT (s)->istep = 1;
}

GfsSourceGenericClass * gfs_source_generic_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_generic_info = {
      "GfsSourceGeneric",
      sizeof (GfsSourceGeneric),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) source_generic_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_source_generic_info);
  }

  return klass;
}

/* GfsSourceScalar: Object */

static void source_scalar_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_source_scalar_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_source_scalar_class ())->parent_class->write) 
      (o, fp);

  g_assert (GFS_SOURCE_SCALAR (o)->v);
  fprintf (fp, " %s", GFS_SOURCE_SCALAR (o)->v->name);
}

static void source_scalar_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceScalar * source;
  GfsDomain * domain;

  if (GTS_OBJECT_CLASS (gfs_source_scalar_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_scalar_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  source = GFS_SOURCE_SCALAR (*o);
  domain =  GFS_DOMAIN (gfs_object_simulation (source));
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsVariable)");
    return;
  }
  source->v = gfs_variable_from_name (domain->variables, 
				      fp->token->str);
  if (source->v == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  if (source->v->sources == NULL)
    source->v->sources = 
      gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
  gts_container_add (source->v->sources, GTS_CONTAINEE (source));
  
  gts_file_next_token (fp);
}

static void source_scalar_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read =  source_scalar_read;
  GTS_OBJECT_CLASS (klass)->write = source_scalar_write;
}

GfsSourceGenericClass * gfs_source_scalar_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_scalar_info = {
      "GfsSourceScalar",
      sizeof (GfsSourceScalar),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_scalar_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()),
				  &gfs_source_scalar_info);
  }

  return klass;
}

/* GfsSourceVelocity: Object */

static void source_velocity_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceVelocity * source;
  GfsDomain * domain;
  FttComponent c;

  if (GTS_OBJECT_CLASS (gfs_source_velocity_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_velocity_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  source = GFS_SOURCE_VELOCITY (*o);
  domain =  GFS_DOMAIN (gfs_object_simulation (source));
  if (!(source->v = gfs_domain_velocity (domain))) {
    gts_file_error (fp, "cannot find velocity components");
    return;
  }
  for (c = 0; c < FTT_DIMENSION; c++) {
    if (source->v[c]->sources == NULL)
      source->v[c]->sources = 
	gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
    gts_container_add (source->v[c]->sources, GTS_CONTAINEE (source));
  }
}

static void source_velocity_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read =  source_velocity_read;
}

GfsSourceGenericClass * gfs_source_velocity_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_info = {
      "GfsSourceGeneric",
      sizeof (GfsSourceVelocity),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_velocity_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()),
				  &source_info);
  }

  return klass;
}

/* GfsSource: Object */

static void source_destroy (GtsObject * o)
{
  if (GFS_SOURCE (o)->intensity)
    gts_object_destroy (GTS_OBJECT (GFS_SOURCE (o)->intensity));

  (* GTS_OBJECT_CLASS (gfs_source_class ())->parent_class->destroy) (o);
}

static void source_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_source_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GFS_SOURCE (*o)->intensity = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_read (GFS_SOURCE (*o)->intensity, gfs_object_simulation (*o), fp);
}

static void source_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_source_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_source_class ())->parent_class->write) 
      (o, fp);
  gfs_function_write (GFS_SOURCE (o)->intensity, fp);
}

static gdouble source_value (GfsSourceGeneric * s, 
			     FttCell * cell, 
			     GfsVariable * v)
{
  return gfs_function_value (GFS_SOURCE (s)->intensity, cell);
}

static void source_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = source_destroy;
  GTS_OBJECT_CLASS (klass)->read = source_read;
  GTS_OBJECT_CLASS (klass)->write = source_write;

  klass->mac_value = klass->centered_value = source_value;
}

GfsSourceGenericClass * gfs_source_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_info = {
      "GfsSource",
      sizeof (GfsSource),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_scalar_class ()),
				  &source_info);
  }

  return klass;
}

/* GfsSourceControl: Object */

static gboolean source_control_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* gfs_event_class ()->event) (event, sim)) {
    GfsSourceControl * s = GFS_SOURCE_CONTROL (event);
    GtsRange r = gfs_domain_stats_variable (GFS_DOMAIN (sim), GFS_SOURCE_SCALAR (event)->v,
					    FTT_TRAVERSE_LEAFS, -1);
    s->s = (gfs_function_value (GFS_SOURCE (s)->intensity, NULL) - r.mean)/
      sim->advection_params.dt;
    return TRUE;
  }
  return FALSE;
}

static gdouble source_control_value (GfsSourceGeneric * s, 
				     FttCell * cell, 
				     GfsVariable * v)
{
  return GFS_SOURCE_CONTROL (s)->s;
}

static void source_control_class_init (GfsSourceGenericClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = source_control_event;
  klass->mac_value = klass->centered_value = source_control_value;
}

GfsSourceGenericClass * gfs_source_control_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_control_info = {
      "GfsSourceControl",
      sizeof (GfsSourceControl),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_control_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_class ()),
				  &source_control_info);
  }

  return klass;
}

/* GfsDiffusion: Object */

static void diffusion_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_DIFFUSION (o)->val));

  (* GTS_OBJECT_CLASS (gfs_diffusion_class ())->parent_class->destroy) (o);
}

static void diffusion_read (GtsObject ** o, GtsFile * fp)
{
  gfs_function_read (GFS_DIFFUSION (*o)->val, gfs_object_simulation (*o), fp);
}

static void diffusion_write (GtsObject * o, FILE * fp)
{
  gfs_function_write (GFS_DIFFUSION (o)->val, fp);
}

static gdouble diffusion_face (GfsDiffusion * d, FttCellFace * f)
{
  return gfs_function_face_value (d->val, f);
}

static gdouble diffusion_cell (GfsDiffusion * d, FttCell * cell)
{
  return gfs_function_value (d->val, cell);
}

static void diffusion_class_init (GfsDiffusionClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = diffusion_destroy;
  GTS_OBJECT_CLASS (klass)->read = diffusion_read;
  GTS_OBJECT_CLASS (klass)->write = diffusion_write;
  GFS_EVENT_CLASS (klass)->event = NULL;
  klass->face = diffusion_face;
  klass->cell = diffusion_cell;
}

static void diffusion_init (GfsDiffusion * d)
{
  d->val = gfs_function_new (gfs_function_class (), 0.);
}

GfsDiffusionClass * gfs_diffusion_class (void)
{
  static GfsDiffusionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo diffusion_info = {
      "GfsDiffusion",
      sizeof (GfsDiffusion),
      sizeof (GfsDiffusionClass),
      (GtsObjectClassInitFunc) diffusion_class_init,
      (GtsObjectInitFunc) diffusion_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &diffusion_info);
  }

  return klass;
}

gdouble gfs_diffusion_face (GfsDiffusion * d, FttCellFace * f)
{
  return (* GFS_DIFFUSION_CLASS (GTS_OBJECT (d)->klass)->face) (d, f);
}

gdouble gfs_diffusion_cell (GfsDiffusion * d, FttCell * cell)
{
  return (* GFS_DIFFUSION_CLASS (GTS_OBJECT (d)->klass)->cell) (d, cell);
}

/* GfsDiffusionMulti: Object */

static void diffusion_multi_destroy (GtsObject * object)
{
  g_slist_foreach (GFS_DIFFUSION_MULTI (object)->d, (GFunc) gts_object_destroy, NULL);
  g_slist_free (GFS_DIFFUSION_MULTI (object)->d);

  (* GTS_OBJECT_CLASS (gfs_diffusion_multi_class ())->parent_class->destroy) 
    (object);
}

static void diffusion_multi_read (GtsObject ** o, GtsFile * fp)
{
  GfsDiffusionMulti * m = GFS_DIFFUSION_MULTI (*o);
  gboolean constant = TRUE;
  
  while (fp->type != GTS_ERROR && fp->type != '\n') {
    GtsObject * d = gts_object_new (GTS_OBJECT_CLASS (gfs_diffusion_class ()));

    (* d->klass->read) (&d, fp);
    if (m->d || d->klass != GTS_OBJECT_CLASS (gfs_diffusion_class ()))
      constant = FALSE;
    m->d = g_slist_prepend (m->d, d);
  }
  m->d = g_slist_reverse (m->d);
  if (fp->type != GTS_ERROR && m->d->next)
    m->c = gfs_variable_from_name (GFS_DOMAIN (gfs_object_simulation (m))->variables, "C");
  if (!constant) {
    m->mu = gfs_domain_add_variable (GFS_DOMAIN (gfs_object_simulation (m)), 
				     "_gfs_diffusion_multi");
    g_assert (m->mu);
  }
}

static void diffusion_multi_write (GtsObject * o, FILE * fp)
{
  GSList * i = GFS_DIFFUSION_MULTI (o)->d;

  while (i) {
    (* GTS_OBJECT (i->data)->klass->write) (i->data, fp);
    i = i->next;
  }
}

static gdouble diffusion_multi_face (GfsDiffusion * d, FttCellFace * f)
{
  gdouble mu1 = gfs_diffusion_face (GFS_DIFFUSION_MULTI (d)->d->data, f);

  if (!GFS_DIFFUSION_MULTI (d)->d->next)
    return mu1;
  else {
    /* fixme: c should be evaluated at t or t + dt/2 */
    gdouble c = gfs_face_interpolated_value (f, GFS_DIFFUSION_MULTI (d)->c->i);
    gdouble mu2 = gfs_diffusion_face (GFS_DIFFUSION_MULTI (d)->d->next->data, f);

    return mu1 + c*(mu2 - mu1);
  }
}

static gdouble diffusion_multi_cell (GfsDiffusion * d, FttCell * cell)
{
  gdouble mu1 = gfs_diffusion_cell (GFS_DIFFUSION_MULTI (d)->d->data, cell);

  if (!GFS_DIFFUSION_MULTI (d)->d->next)
    return mu1;
  else {
    /* fixme: c should be evaluated at t or t + dt/2 */
    gdouble c = GFS_VARIABLE (cell, GFS_DIFFUSION_MULTI (d)->c->i);
    gdouble mu2 = gfs_diffusion_cell (GFS_DIFFUSION_MULTI (d)->d->next->data, cell);
    
    return mu1 + c*(mu2 - mu1);
  }
}

static void set_mu (FttCell * cell, GfsDiffusion * d)
{
  GFS_VARIABLE (cell, GFS_DIFFUSION_MULTI (d)->mu->i) = diffusion_multi_cell (d, cell);
}

static gboolean diffusion_multi_event (GfsEvent * event, GfsSimulation * sim)
{
  GfsDiffusionMulti * m = GFS_DIFFUSION_MULTI (event);
  GSList * i = m->d;

  while (i) {
    if (GFS_EVENT_CLASS (GTS_OBJECT (i->data)->klass)->event)
      (* GFS_EVENT_CLASS (GTS_OBJECT (i->data)->klass)->event) (event, sim);
    i = i->next;
  }

  if (m->mu) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) set_mu, m);
    gfs_domain_bc (GFS_DOMAIN (sim), FTT_TRAVERSE_LEAFS, -1, m->mu);
  }
  
  return TRUE;
}

static void diffusion_multi_class_init (GfsDiffusionClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = diffusion_multi_read;
  GTS_OBJECT_CLASS (klass)->write = diffusion_multi_write;
  GTS_OBJECT_CLASS (klass)->destroy = diffusion_multi_destroy;
  GFS_EVENT_CLASS (klass)->event = diffusion_multi_event;
  klass->face = diffusion_multi_face;
  klass->cell = diffusion_multi_cell;
}

GfsDiffusionClass * gfs_diffusion_multi_class (void)
{
  static GfsDiffusionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo diffusion_multi_info = {
      "GfsDiffusionMulti",
      sizeof (GfsDiffusionMulti),
      sizeof (GfsDiffusionClass),
      (GtsObjectClassInitFunc) diffusion_multi_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_diffusion_class ()),
				  &diffusion_multi_info);
  }

  return klass;
}

/* GfsSourceDiffusion: Object */

static void source_diffusion_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_SOURCE_DIFFUSION (o)->D));

  (* GTS_OBJECT_CLASS (gfs_source_diffusion_class ())->parent_class->destroy) (o);
}

static GfsSourceDiffusion * previous_diffusion_source (GfsVariable * v,
						       GfsSourceDiffusion * d)
{
  GSList * i;

  i = GTS_SLIST_CONTAINER (v->sources)->items;
  while (i) {
    if (i->data != d && GFS_IS_SOURCE_DIFFUSION (i->data))
      return i->data;
    i = i->next;
  }
  return NULL;
}

static void source_diffusion_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceDiffusion * d;

  if (GTS_OBJECT_CLASS (gfs_source_diffusion_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_diffusion_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  d = GFS_SOURCE_DIFFUSION (*o);
  if (previous_diffusion_source (GFS_SOURCE_SCALAR (d)->v, d)) {
    gts_file_error (fp, "only one diffusion source can be specified");
    return;
  }

  gfs_object_simulation_set (d->D, gfs_object_simulation (d));
  (* GTS_OBJECT (d->D)->klass->read) ((GtsObject **) &d->D, fp);
}

static void source_diffusion_write (GtsObject * o, FILE * fp)
{
  GfsSourceDiffusion * d = GFS_SOURCE_DIFFUSION (o);

  (* GTS_OBJECT_CLASS (gfs_source_diffusion_class ())->parent_class->write) (o, fp);
  (* GTS_OBJECT (d->D)->klass->write) (GTS_OBJECT (d->D), fp);
}

static gboolean source_diffusion_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* gfs_event_class ()->event) (event, sim)) {
    GfsSourceDiffusion * d = GFS_SOURCE_DIFFUSION (event);

    if ((* GFS_EVENT_CLASS (GTS_OBJECT (d->D)->klass)->event))
      (* GFS_EVENT_CLASS (GTS_OBJECT (d->D)->klass)->event) (GFS_EVENT (d->D), sim);
    return TRUE;
  }
  return FALSE;
}

static gdouble source_diffusion_value (GfsSourceGeneric * s, 
				       FttCell * cell,
				       GfsVariable * v)
{
  FttCellFace f;
  FttCellNeighbors n;
  GfsGradient g = { 0., 0. };
  FttComponent c;
  gdouble v0, h;

  if (GFS_IS_MIXED (cell)) /* this improves results for channel test */
    return 0.;

  c = v->component;

  v0 = GFS_VARIABLE (cell, v->i);
  f.cell = cell;
  ftt_cell_neighbors (cell, &n);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    gdouble D;

    f.neighbor = n.c[f.d];
    D = gfs_source_diffusion_face (GFS_SOURCE_DIFFUSION (s), &f);
    if (f.neighbor) {
      GfsGradient e;

      gfs_face_gradient (&f, &e, v->i, -1);
      g.a += D*e.a;
      g.b += D*e.b;
    }
    else if (f.d/2 == c) {
      g.a += D;
      g.b -= D*v0;
    }
  }
  h = ftt_cell_size (cell);
  return (g.b - g.a*v0)/(h*h);
}

static void source_diffusion_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = source_diffusion_destroy;
  GTS_OBJECT_CLASS (klass)->read = source_diffusion_read;
  GTS_OBJECT_CLASS (klass)->write = source_diffusion_write;

  GFS_EVENT_CLASS (klass)->event = source_diffusion_event;

  klass->mac_value = source_diffusion_value;
}

static void source_diffusion_init (GfsSourceDiffusion * d)
{
  d->D = GFS_DIFFUSION (gts_object_new (GTS_OBJECT_CLASS (gfs_diffusion_class ())));
}

GfsSourceGenericClass * gfs_source_diffusion_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_diffusion_info = {
      "GfsSourceDiffusion",
      sizeof (GfsSourceDiffusion),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_diffusion_class_init,
      (GtsObjectInitFunc) source_diffusion_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_scalar_class ()),
				  &source_diffusion_info);
  }

  return klass;
}

gdouble gfs_source_diffusion_face (GfsSourceDiffusion * d, FttCellFace * f)
{
  g_return_val_if_fail (d != NULL, 0.);
  g_return_val_if_fail (f != NULL, 0.);

  return gfs_diffusion_face (d->D, f);
}

gdouble gfs_source_diffusion_cell (GfsSourceDiffusion * d, FttCell * cell)
{
  g_return_val_if_fail (d != NULL, 0.);
  g_return_val_if_fail (cell != NULL, 0.);
  g_return_val_if_fail (GFS_IS_MIXED (cell), 0.);

  return gfs_diffusion_cell (d->D, cell);
}

/* GfsSourceDiffusionExplicit: Object */

static void explicit_diffusion (FttCell * cell, GfsSourceGeneric * s)
{
  GFS_VARIABLE (cell, GFS_SOURCE_DIFFUSION_EXPLICIT (s)->s->i) = 
    source_diffusion_value (s, cell, s->v);
}

static gboolean gfs_source_diffusion_explicit_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_source_diffusion_explicit_class ())->parent_class)->event) (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) explicit_diffusion, event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_source_diffusion_explicit_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_source_diffusion_explicit_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_diffusion_explicit_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GFS_SOURCE_DIFFUSION_EXPLICIT (*o)->s = 
    gfs_temporary_variable (GFS_DOMAIN (gfs_object_simulation (*o)));
}

static void gfs_source_diffusion_explicit_destroy (GtsObject * o)
{
  if (GFS_SOURCE_DIFFUSION_EXPLICIT (o)->s)
    gts_object_destroy (GTS_OBJECT (GFS_SOURCE_DIFFUSION_EXPLICIT (o)->s));

  (* GTS_OBJECT_CLASS (gfs_source_diffusion_explicit_class ())->parent_class->destroy) (o);
}

static gdouble source_diffusion_explicit_value (GfsSourceGeneric * s, 
					     FttCell * cell,
					     GfsVariable * v)
{
  return GFS_VARIABLE (cell, GFS_SOURCE_DIFFUSION_EXPLICIT (s)->s->i);
}

static void gfs_source_diffusion_explicit_class_init (GfsSourceGenericClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_source_diffusion_explicit_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_source_diffusion_explicit_read;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_source_diffusion_explicit_destroy;
  klass->mac_value = klass->centered_value = source_diffusion_explicit_value;
}

GfsSourceGenericClass * gfs_source_diffusion_explicit_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_diffusion_explicit_info = {
      "GfsSourceDiffusionExplicit",
      sizeof (GfsSourceDiffusionExplicit),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_diffusion_explicit_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_diffusion_class ()),
				  &gfs_source_diffusion_explicit_info);
  }

  return klass;
}

/* GfsSourceViscosity: Object */

static void source_viscosity_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_SOURCE_VISCOSITY (o)->D));

  (* GTS_OBJECT_CLASS (gfs_source_viscosity_class ())->parent_class->destroy) (o);
}

static void source_viscosity_read (GtsObject ** o, GtsFile * fp)
{
  FttComponent c;
  GfsSourceViscosity * s;

  (* GTS_OBJECT_CLASS (gfs_source_viscosity_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  for (c = 0; c < FTT_DIMENSION; c++) {
    GfsVariable * v = GFS_SOURCE_VELOCITY (*o)->v[c];

    if (v->sources) {
      GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
      
      while (i) {
	if (i->data != *o && GFS_IS_SOURCE_DIFFUSION (i->data)) {
	  gts_file_error (fp, "variable '%s' cannot have multiple diffusion source terms", v->name);
	  return;
	}
	i = i->next;
      }
    }
  }

  s = GFS_SOURCE_VISCOSITY (*o);
  gfs_object_simulation_set (s->D, gfs_object_simulation (s));
  (* GTS_OBJECT (s->D)->klass->read) ((GtsObject **) &s->D, fp);
}

static void source_viscosity_write (GtsObject * o, FILE * fp)
{
  GfsSourceViscosity * s = GFS_SOURCE_VISCOSITY (o);

  (* GTS_OBJECT_CLASS (gfs_source_viscosity_class ())->parent_class->write) (o, fp);
  (* GTS_OBJECT (s->D)->klass->write) (GTS_OBJECT (s->D), fp);
}

static gdouble source_viscosity_non_diffusion_value (GfsSourceGeneric * s,
						     FttCell * cell,
						     GfsVariable * v)
{
  GfsVariable * mu = GFS_DIFFUSION_MULTI (GFS_SOURCE_DIFFUSION (s)->D)->mu;

  if (mu == NULL)
    return 0.;
  else {
    GfsVariable ** u = GFS_SOURCE_VELOCITY (s)->v;
    FttComponent c = v->component, j;
    gdouble rho = 1.
      /* fixme: + GFS_STATE (cell)->c*(gfs_object_simulation (s)->physical_params.rho - 1.)*/;
    gdouble h = ftt_cell_size (cell);
    gdouble a = 0.;

    for (j = 0; j < FTT_DIMENSION; j++)
      a += (gfs_center_gradient (cell, c, u[j]->i)*
	    gfs_center_gradient (cell, j, mu->i));
    return a/(rho*h*h);
  }
}

static gdouble source_viscosity_value (GfsSourceGeneric * s,
				       FttCell * cell,
				       GfsVariable * v)
{
  gdouble rho = 1.
    /* fixme: + GFS_STATE (cell)->c*(gfs_object_simulation (s)->physical_params.rho - 1.)*/;

  return (source_diffusion_value (s, cell, v)/rho +
	  source_viscosity_non_diffusion_value (s, cell, v));
}

static gboolean source_viscosity_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* gfs_event_class ()->event) (event, sim)) {
    GfsSourceViscosity * s = GFS_SOURCE_VISCOSITY (event);

    if ((* GFS_EVENT_CLASS (GTS_OBJECT (s->D)->klass)->event))
      (* GFS_EVENT_CLASS (GTS_OBJECT (s->D)->klass)->event) (GFS_EVENT (s->D), sim);
    return TRUE;
  }
  return FALSE;
}

static void source_viscosity_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = source_viscosity_destroy;
  GTS_OBJECT_CLASS (klass)->read = source_viscosity_read;
  GTS_OBJECT_CLASS (klass)->write = source_viscosity_write;
  
  GFS_EVENT_CLASS (klass)->event = source_viscosity_event;

  klass->mac_value = source_viscosity_value;
  klass->centered_value = source_viscosity_non_diffusion_value;
}

static void source_viscosity_init (GfsSourceViscosity * s)
{
  s->D = GFS_DIFFUSION (gts_object_new (GTS_OBJECT_CLASS (gfs_diffusion_class ())));
}

GfsSourceGenericClass * gfs_source_viscosity_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_viscosity_info = {
      "GfsSourceViscosity",
      sizeof (GfsSourceViscosity),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_viscosity_class_init,
      (GtsObjectInitFunc) source_viscosity_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
				  &source_viscosity_info);
  }

  return klass;
}

/* GfsSourceCoriolis: Object */

static void source_coriolis_destroy (GtsObject * o)
{
  FttComponent c;

  if (GFS_SOURCE_CORIOLIS (o)->omegaz)
    gts_object_destroy (GTS_OBJECT (GFS_SOURCE_CORIOLIS (o)->omegaz));
  if (GFS_SOURCE_CORIOLIS (o)->drag)
    gts_object_destroy (GTS_OBJECT (GFS_SOURCE_CORIOLIS (o)->drag));

  for (c = 0; c <  2; c++)
    if (GFS_SOURCE_CORIOLIS (o)->u[c])
      gts_object_destroy (GTS_OBJECT (GFS_SOURCE_CORIOLIS (o)->u[c]));

  (* GTS_OBJECT_CLASS (gfs_source_class ())->parent_class->destroy) (o);
}

static void gfs_source_coriolis_read (GtsObject ** o, GtsFile * fp)
{
  FttComponent c;
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  (* GTS_OBJECT_CLASS (gfs_source_coriolis_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  for (c = 0; c < FTT_DIMENSION; c++) {
    GfsVariable * v = GFS_SOURCE_VELOCITY (*o)->v[c];

    if (v->sources) {
      GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
      
      while (i) {
	if (i->data != *o && GFS_IS_SOURCE_CORIOLIS (i->data)) {
	  gts_file_error (fp, "variable '%s' cannot have multiple Coriolis source terms", v->name);
	  return;
	}
	i = i->next;
      }
    }
  }

  GFS_SOURCE_CORIOLIS (*o)->omegaz = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_read (GFS_SOURCE_CORIOLIS (*o)->omegaz, gfs_object_simulation (*o), fp);

  if (fp->type != '\n') {
    GFS_SOURCE_CORIOLIS (*o)->drag = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (GFS_SOURCE_CORIOLIS (*o)->drag, gfs_object_simulation (*o), fp);
  }

#if (!FTT_2D)
  gts_container_remove (GFS_SOURCE_VELOCITY (*o)->v[FTT_Z]->sources, GTS_CONTAINEE (*o));
#endif /* 3D */ 
 
  for (c = 0; c <  2; c++)
    GFS_SOURCE_CORIOLIS (*o)->u[c] = gfs_temporary_variable (domain);
}

static void gfs_source_coriolis_write (GtsObject * o, FILE * fp)
{
  GfsSourceCoriolis * s = GFS_SOURCE_CORIOLIS (o);

  (* GTS_OBJECT_CLASS (gfs_source_coriolis_class ())->parent_class->write) (o, fp);
  gfs_function_write (s->omegaz, fp);
  if (s->drag)
    gfs_function_write (s->drag, fp);
}

static gdouble gfs_source_coriolis_mac_value (GfsSourceGeneric * s,
					      FttCell * cell,
					      GfsVariable * v)
{
  GfsSourceVelocity * sv = GFS_SOURCE_VELOCITY (s);
  gdouble f;

  f = gfs_function_value (GFS_SOURCE_CORIOLIS (s)->omegaz, cell);
  switch (v->component) {
  case FTT_X: return   f*GFS_VARIABLE (cell, sv->v[1]->i);
  case FTT_Y: return - f*GFS_VARIABLE (cell, sv->v[0]->i);
  default: g_assert_not_reached ();
  }
  return 0.;
}

static void save_coriolis (FttCell * cell, GfsSourceCoriolis * s)
{
  GfsSourceVelocity * sv = GFS_SOURCE_VELOCITY (s);
  FttComponent c;
  gdouble f;

  f = gfs_function_value (s->omegaz, cell)/2.;
  for (c = 0; c < 2; c++)
    GFS_VARIABLE (cell, s->u[c]->i) = c == FTT_X ?
      f*GFS_VARIABLE (cell, sv->v[1]->i) :
      - f*GFS_VARIABLE (cell, sv->v[0]->i);
}

static gboolean gfs_source_coriolis_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class)->event) (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) save_coriolis, event);
    return TRUE;
  }
  return FALSE;
}

static gdouble gfs_source_coriolis_centered_value (GfsSourceGeneric * s,
						   FttCell * cell,
						   GfsVariable * v)
{
  return GFS_VARIABLE (cell, GFS_SOURCE_CORIOLIS (s)->u[v->component]->i);
}

static void gfs_source_coriolis_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = source_coriolis_destroy;
  GTS_OBJECT_CLASS (klass)->read = gfs_source_coriolis_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_source_coriolis_write;

  GFS_EVENT_CLASS (klass)->event = gfs_source_coriolis_event;
  klass->mac_value = gfs_source_coriolis_mac_value;
  klass->centered_value = gfs_source_coriolis_centered_value;
}

GfsSourceGenericClass * gfs_source_coriolis_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_coriolis_info = {
      "GfsSourceCoriolis",
      sizeof (GfsSourceCoriolis),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_coriolis_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
				  &gfs_source_coriolis_info);
  }

  return klass;
}

/**
 * gfs_has_source_coriolis:
 * @domain: a #GfsDomain.
 *
 * Returns: the #GfsSourceCoriolis associated with @domain or %NULL.
 */
GfsSourceCoriolis * gfs_has_source_coriolis (GfsDomain * domain)
{
  GfsVariable * v;

  g_return_val_if_fail (domain != NULL, NULL);

  v = gfs_variable_from_name (domain->variables, "U");
  g_return_val_if_fail (v != NULL, NULL);

  if (v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;

    while (i) {
      if (GFS_IS_SOURCE_CORIOLIS (i->data))
	return i->data;
      i = i->next;
    }
  }
  return NULL;
}

static void implicit_coriolis (FttCell * cell, GfsSourceCoriolis * s)
{
  GfsSourceVelocity * sv = GFS_SOURCE_VELOCITY (s);
  gdouble c, u, v;
  GfsSimulation * sim = gfs_object_simulation (s);

  c = sim->advection_params.dt*gfs_function_value (s->omegaz, cell)/2.;
  u = GFS_VARIABLE (cell, sv->v[0]->i);
  v = GFS_VARIABLE (cell, sv->v[1]->i);
  if (s->drag) {
    gdouble e = sim->advection_params.dt*gfs_function_value (s->drag, cell)/2.;
    GFS_VARIABLE (cell, sv->v[0]->i) = (u + c*v/(1. + e))/((1. + e) + c*c/(1. + e));
    GFS_VARIABLE (cell, sv->v[1]->i) = (v - c*u/(1. + e))/((1. + e) + c*c/(1. + e));
  }
  else {
    GFS_VARIABLE (cell, sv->v[0]->i) = (u + c*v)/(1. + c*c);
    GFS_VARIABLE (cell, sv->v[1]->i) = (v - c*u)/(1. + c*c);
  }
}

/**
 * gfs_source_coriolis_implicit:
 * @domain: a #GfsDomain.
 * @dt: the timestep.
 *
 * Applies the implicit part of the Coriolis source term of @domain.
 */
void gfs_source_coriolis_implicit (GfsDomain * domain,
				   gdouble dt)
{
  GfsSourceCoriolis * s;

  g_return_if_fail (domain != NULL);

  if ((s = gfs_has_source_coriolis (domain)))
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) implicit_coriolis, s);
}
