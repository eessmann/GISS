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

static void source_generic_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_source_generic_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_source_generic_class ())->parent_class->write) 
      (o, fp);

  g_assert (GFS_SOURCE_GENERIC (o)->v);
  fprintf (fp, " %s", GFS_SOURCE_GENERIC (o)->v->name);
}

static void source_generic_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceGeneric * source;
  GfsDomain * domain;

  if (GTS_OBJECT_CLASS (gfs_source_generic_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_generic_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  source = GFS_SOURCE_GENERIC (*o);
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

static void source_generic_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read =  source_generic_read;
  GTS_OBJECT_CLASS (klass)->write = source_generic_write;
}

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
      (GtsObjectClassInitFunc) source_generic_class_init,
      (GtsObjectInitFunc) source_generic_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_source_generic_info);
  }

  return klass;
}

/* GfsSourceVector: Object */

static void source_vector_write (GtsObject * o, FILE * fp)
{
  FttComponent c;

  if (GTS_OBJECT_CLASS (gfs_source_vector_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_source_vector_class ())->parent_class->write) 
      (o, fp);

  for (c = 1; c < FTT_DIMENSION; c++) {
    g_assert (GFS_SOURCE_VECTOR (o)->v[c]);
    fprintf (fp, " %s", GFS_SOURCE_VECTOR (o)->v[c]->name);
  }
}

static void source_vector_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceVector * source;
  GfsDomain * domain;
  FttComponent c;

  if (GTS_OBJECT_CLASS (gfs_source_vector_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_vector_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  source = GFS_SOURCE_VECTOR (*o);
  domain =  GFS_DOMAIN (gfs_object_simulation (source));
  for (c = 1; c < FTT_DIMENSION; c++) {
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a string (GfsVariable)");
      return;
    }
    source->v[c] = gfs_variable_from_name (domain->variables, 
					   fp->token->str);
    if (source->v[c] == NULL) {
      gts_file_error (fp, "unknown variable `%s'", fp->token->str);
      return;
    }
    if (source->v[c]->sources == NULL)
      source->v[c]->sources = 
	gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
    gts_container_add (source->v[c]->sources, GTS_CONTAINEE (source));
  
    gts_file_next_token (fp);
  }
}

static void source_vector_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read =  source_vector_read;
  GTS_OBJECT_CLASS (klass)->write = source_vector_write;
}

GfsSourceGenericClass * gfs_source_vector_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_info = {
      "GfsSourceGeneric",
      sizeof (GfsSourceVector),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_vector_class_init,
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
  gfs_object_simulation (GFS_SOURCE (*o)->intensity) = gfs_object_simulation (*o);
  gfs_function_read (GFS_SOURCE (*o)->intensity, fp);
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
  FttVector p;

  if (v->centered)
    ftt_cell_pos (cell, &p);
  else
    gfs_cell_cm (cell, &p);
  return gfs_function_value (GFS_SOURCE (s)->intensity, cell, &p,
			     gfs_object_simulation (s)->time.t);
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
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()),
				  &source_info);
  }

  return klass;
}

/* GfsDiffusion: Object */

static void diffusion_read (GtsObject ** o, GtsFile * fp)
{
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (D)");
    return;
  }
  GFS_DIFFUSION (*o)->val = atof (fp->token->str);
  gts_file_next_token (fp);
}

static void diffusion_write (GtsObject * o, FILE * fp)
{
  fprintf (fp, " %g", GFS_DIFFUSION (o)->val);
}

static gdouble diffusion_face (GfsDiffusion * d, FttCellFace * f)
{
  return d->val;
}

static gdouble diffusion_cell (GfsDiffusion * d, FttCell * cell)
{
  return d->val;
}

static void diffusion_class_init (GfsDiffusionClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = diffusion_read;
  GTS_OBJECT_CLASS (klass)->write = diffusion_write;
  GFS_EVENT_CLASS (klass)->event = NULL;
  klass->face = diffusion_face;
  klass->cell = diffusion_cell;
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
      (GtsObjectInitFunc) NULL,
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
  if (previous_diffusion_source (GFS_SOURCE_GENERIC (d)->v, d)) {
    gts_file_error (fp, "only one diffusion source can be specified");
    return;
  }

  gfs_object_simulation (d->D) = gfs_object_simulation (d);
  (* GTS_OBJECT (d->D)->klass->read) ((GtsObject **) &d->D, fp);
}

static void source_diffusion_write (GtsObject * o, FILE * fp)
{
  GfsSourceDiffusion * d = GFS_SOURCE_DIFFUSION (o);

  if (GTS_OBJECT_CLASS (gfs_source_diffusion_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_source_diffusion_class ())->parent_class->write) 
      (o, fp);
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

  c = GFS_VELOCITY_COMPONENT (v->i);

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
  d->D = GFS_DIFFUSION (gts_object_new (GTS_OBJECT_CLASS (gfs_diffusion_multi_class ())));
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
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()),
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
    gfs_domain_add_variable (GFS_DOMAIN (gfs_object_simulation (*o)), NULL);
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

static void source_viscosity_read (GtsObject ** o, GtsFile * fp)
{
  GfsVariable * v;
  FttComponent c;

  if (GTS_OBJECT_CLASS (gfs_source_viscosity_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_viscosity_class ())->parent_class->read)
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  
  v = GFS_SOURCE_GENERIC (*o)->v->next;
  for (c = 1; c < FTT_DIMENSION; c++, v = v->next)
    if (!v) {
      gts_file_error (fp, "not enough velocity components");
      return;
    }
    else {
      if (v->sources == NULL)
	v->sources = 
	  gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
      else if (previous_diffusion_source (v, NULL)) {
	gts_file_error (fp, "only one diffusion source can be specified for a given variable");
	return;
      }
      gts_container_add (v->sources, GTS_CONTAINEE (*o));
    }
}

static gdouble source_viscosity_non_diffusion_value (GfsSourceGeneric * s,
						     FttCell * cell,
						     GfsVariable * v)
{
  GfsVariable * mu = GFS_DIFFUSION_MULTI (GFS_SOURCE_DIFFUSION (s)->D)->mu;

  if (mu == NULL)
    return 0.;
  else {
    FttComponent c = GFS_VELOCITY_COMPONENT (v->i), i;
    gdouble rho = 1.
      /* fixme: + GFS_STATE (cell)->c*(gfs_object_simulation (s)->physical_params.rho - 1.)*/;
    gdouble h = ftt_cell_size (cell);
    gdouble a = 0.;

    for (i = 0; i < FTT_DIMENSION; i++)
      a += (gfs_center_gradient (cell, c, GFS_VELOCITY_INDEX (i))*
	    gfs_center_gradient (cell, i, mu->i));
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

static void source_viscosity_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = source_viscosity_read;
  klass->mac_value = source_viscosity_value;
  klass->centered_value = source_viscosity_non_diffusion_value;
}

GfsSourceGenericClass * gfs_source_viscosity_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_viscosity_info = {
      "GfsSourceViscosity",
      sizeof (GfsSourceDiffusion),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_viscosity_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_diffusion_class ()),
				  &source_viscosity_info);
  }

  return klass;
}

/* GfsSourceCoriolis: Object */

static void source_coriolis_destroy (GtsObject * o)
{
  if (GFS_SOURCE_CORIOLIS (o)->omegaz)
    gts_object_destroy (GTS_OBJECT (GFS_SOURCE_CORIOLIS (o)->omegaz));

  (* GTS_OBJECT_CLASS (gfs_source_class ())->parent_class->destroy) (o);
}

static void gfs_source_coriolis_read (GtsObject ** o, GtsFile * fp)
{
  FttComponent c;
  GfsVariable * v;
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  if (GTS_OBJECT_CLASS (gfs_source_coriolis_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_coriolis_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GFS_SOURCE_CORIOLIS (*o)->omegaz = gfs_function_new (gfs_function_class (), 0.);
  gfs_object_simulation (GFS_SOURCE_CORIOLIS (*o)->omegaz) = gfs_object_simulation (*o);
  gfs_function_read (GFS_SOURCE_CORIOLIS (*o)->omegaz, fp);

  v = GFS_SOURCE_GENERIC (*o)->v->next;
  for (c = 1; c < 2; c++, v = v->next) {
    if (!v) {
      gts_file_error (fp, "not enough velocity components");
      return;
    }
    else {
      if (v->sources == NULL)
	v->sources = gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
      gts_container_add (v->sources, GTS_CONTAINEE (*o));
    }
  }
  for (c = 0; c <  2; c++) {
    GFS_SOURCE_CORIOLIS (*o)->u[c] = gfs_domain_add_variable (domain, NULL);
    g_assert (GFS_SOURCE_CORIOLIS (*o)->u[c]);
  }
}

static void gfs_source_coriolis_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_source_coriolis_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_source_coriolis_class ())->parent_class->write) (o, fp);
  gfs_function_write (GFS_SOURCE_CORIOLIS (o)->omegaz, fp);
}

static gdouble gfs_source_coriolis_mac_value (GfsSourceGeneric * s,
					      FttCell * cell,
					      GfsVariable * v)
{
  FttVector p;
  gdouble f;

  gfs_cell_cm (cell, &p);
  f = gfs_function_value (GFS_SOURCE_CORIOLIS (s)->omegaz, NULL, &p, 
			  gfs_object_simulation (s)->time.t);
  switch (GFS_VELOCITY_COMPONENT (v->i)) {
  case FTT_X: return   f*GFS_STATE (cell)->v;
  case FTT_Y: return - f*GFS_STATE (cell)->u;
  default: g_assert_not_reached ();
  }
  return 0.;
}

static void save_coriolis (FttCell * cell, GfsSourceCoriolis * s)
{
  FttComponent c;
  FttVector p;
  gdouble f;

  gfs_cell_cm (cell, &p);
  f = gfs_function_value (s->omegaz, NULL, &p, gfs_object_simulation (s)->time.t)/2.;
  for (c = 0; c < 2; c++)
    GFS_VARIABLE (cell, s->u[c]->i) = c == FTT_X ? f*GFS_STATE (cell)->v : -f*GFS_STATE (cell)->u;
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
  FttComponent c = GFS_VELOCITY_COMPONENT (v->i);
  GfsSourceCoriolis * b = GFS_SOURCE_CORIOLIS (s);

  return GFS_VARIABLE (cell, b->u[c]->i);
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
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()),
				  &gfs_source_coriolis_info);
  }

  return klass;
}

static void implicit_coriolis (FttCell * cell, GfsSourceCoriolis * s)
{
  FttVector p;
  gdouble c, u, v;
  GfsSimulation * sim = gfs_object_simulation (s);

  gfs_cell_cm (cell, &p);
  c = sim->advection_params.dt*gfs_function_value (s->omegaz, NULL, &p, sim->time.t)/2.;
  u = GFS_STATE (cell)->u;
  v = GFS_STATE (cell)->v;
  GFS_STATE (cell)->u = (u + c*v)/(1. + c*c);
  GFS_STATE (cell)->v = (v - c*u)/(1. + c*c);
}

/**
 * gfs_source_coriolis_implicit:
 * @sim: a #GfsSimulation.
 * @apar: the advection parameters.
 * @p: the pressure at time n.
 *
 * Applies the implicit part of the Coriolis source term of @sim.
 *
 * Returns: %TRUE if the implicit part was applied, %FALSE otherwise.
 */
gboolean gfs_source_coriolis_implicit (GfsSimulation * sim,
				       GfsAdvectionParams * apar,
				       GfsVariable * p)
{
  GfsVariable * v;

  g_return_val_if_fail (sim != NULL, FALSE);
  g_return_val_if_fail (p != NULL, FALSE);

  v = gfs_variable_from_name (GFS_DOMAIN (sim)->variables, "U");
  g_assert (v);
  if (v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
    guint implicit_sources = 0;
    GfsSourceCoriolis * s = NULL;

    while (i) {
      if (GFS_IS_SOURCE_CORIOLIS (i->data)) {
	if (s != NULL) {
	  g_warning ("Multiple Coriolis source terms are not consistent");
	  return FALSE;
	}
	else if (implicit_sources > 0) {
	  g_warning ("Multiple implicit source terms are not consistent");
	  return FALSE;
	}
	else
	  s = i->data;
      }
      else if (!GFS_SOURCE_GENERIC_CLASS (GTS_OBJECT (i->data)->klass)->centered_value) {
	implicit_sources++;
	if (s != NULL || implicit_sources > 1) {
	  g_warning ("Multiple implicit source terms are not consistent");
	  return FALSE;
	}
      }
      i = i->next;
    }

    if (s != NULL) {
      gfs_poisson_coefficients (GFS_DOMAIN (sim), apar->c, apar->rho);
      gfs_correct_normal_velocities (GFS_DOMAIN (sim), 2, p, apar->dt);
      gfs_correct_centered_velocities (GFS_DOMAIN (sim), 2, apar->dt);
      gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) implicit_coriolis, s);
      return TRUE;
    }
  }
  return FALSE;
}
