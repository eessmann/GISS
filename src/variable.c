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
#include "variable.h"
#include "vof.h"

/* GfsVariable: Object */

static void variable_init_domain (GfsVariable * v, GfsDomain * domain)
{
  v->i = gfs_domain_alloc (domain);
  v->domain = domain;
  GTS_OBJECT (v)->reserved = domain;
}

static void gfs_variable_read (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain;
  GfsVariable * v, * old;

  if (GTS_OBJECT_CLASS (gfs_variable_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_variable_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (name)");
    return;
  }
  domain = (*o)->reserved;
  if (gfs_derived_variable_from_name (domain->derived_variables, fp->token->str)) {
    gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
    return;
  }
  v = GFS_VARIABLE1 (*o);
  v->name = g_strdup (fp->token->str);
  gts_file_next_token (fp);

  if ((old = gfs_variable_from_name (domain->variables, v->name))) {
    GSList * i;
    if ((i = g_slist_find (domain->variables_io, old)))
      i->data = v;
    domain->variables = g_slist_remove (domain->variables, old);
    v->units = old->units;
    gts_object_destroy (GTS_OBJECT (old));
  }
  variable_init_domain (v, domain);
  domain->variables = g_slist_append (domain->variables, v);
}

static void gfs_variable_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_variable_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_variable_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s", GFS_VARIABLE1 (o)->name);
}

static void gfs_variable_destroy (GtsObject * object)
{
  GfsVariable * v = GFS_VARIABLE1 (object);

  g_free (v->name);
  g_free (v->description);
  if (v->sources)
    gts_object_destroy (GTS_OBJECT (v->sources));
  if (v->surface_bc)
    gts_object_destroy (GTS_OBJECT (v->surface_bc));
  if (v->default_bc)
    gts_object_destroy (GTS_OBJECT (v->default_bc));
  if (v->domain) {
    gfs_domain_free (v->domain, v->i);
    v->domain->variables = g_slist_remove (v->domain->variables, v);
  }

  (* GTS_OBJECT_CLASS (gfs_variable_class ())->parent_class->destroy) (object);
}

static void gfs_variable_class_init (GfsVariableClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_variable_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_variable_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_variable_destroy;
}

static void gfs_variable_init (GfsVariable * v)
{
  GFS_EVENT (v)->istep = 1;
  v->centered = FALSE;
  v->component = FTT_DIMENSION;
  v->fine_coarse = (GfsVariableFineCoarseFunc) gfs_get_from_below_intensive;
  v->coarse_fine = (GfsVariableFineCoarseFunc) gfs_cell_coarse_fine;
}

GfsVariableClass * gfs_variable_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_info = {
      "GfsVariable",
      sizeof (GfsVariable),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) gfs_variable_class_init,
      (GtsObjectInitFunc) gfs_variable_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()), &gfs_variable_info);
  }

  return klass;
}

/**
 * gfs_variable_new:
 * @klass: a #GfsVariableClass.
 * @domain: a #GfsDomain.
 * @name: the name of the variable or %NULL.
 * @description: the variable description or %NULL.
 *
 * Returns: a newly allocated #GfsVariable or %NULL if a variable
 * named @name already exists in @domain.
 */
GfsVariable * gfs_variable_new (GfsVariableClass * klass,
				GfsDomain * domain,
				const gchar * name,
				const gchar * description)
{
  GfsVariable * v;

  g_return_val_if_fail (klass != NULL, NULL);
  g_return_val_if_fail (domain != NULL, NULL);

  if (name &&
      (gfs_variable_from_name (domain->variables, name) ||
       gfs_derived_variable_from_name (domain->derived_variables, name)))
    return NULL;

  v = GFS_VARIABLE1 (gts_object_new (GTS_OBJECT_CLASS (klass)));
  if (name)
    v->name = g_strdup (name);
  if (description)
    v->description = g_strdup (description);
  variable_init_domain (v, domain);

  return v;
}

/**
 * gfs_variable_from_name:
 * @i: the list of available #GfsVariable.
 * @name: the name of the variable to find.
 *
 * Returns: the #GfsVariable @name or %NULL if this variable name does
 * not exist.  
 */
GfsVariable * gfs_variable_from_name (GSList * i,
				      const gchar * name)
{
  g_return_val_if_fail (name != NULL, NULL);

  while (i && (!GFS_VARIABLE1 (i->data)->name || strcmp (name, GFS_VARIABLE1 (i->data)->name)))
    i = i->next;
  return i ? GFS_VARIABLE1 (i->data) : NULL;
}

/**
 * gfs_variables_from_list:
 * @i: the list of available #GfsVariable.
 * @list: a malloc'ed string containing comma separated variable names.
 * @error: where to return the variable name in case of error.
 *
 * Returns: a list of variables or %NULL in case of error, in which
 * case *@error points to the name of the unknown variable.  
 */
GSList * gfs_variables_from_list (GSList * i,
				  gchar * list,
				  gchar ** error)
{
  gchar * s;
  GSList * var = NULL;

  g_return_val_if_fail (i != NULL, NULL);
  g_return_val_if_fail (error != NULL, NULL);

  s = strtok (list, ",");
  while (s) {
    GfsVariable * v = gfs_variable_from_name (i, s);

    if (v == NULL) {
      *error = s;
      g_slist_free (var);
      return NULL;
    }
    var = g_slist_append (var, v);
    s = strtok (NULL, ",");
  }
  return var;
}

/**
 * gfs_variables_swap:
 * @v1: a #GfsVariable.
 * @v2: a #GfsVariable.
 *
 * Swaps the values of @v1 and @v2, belonging to the same #GfsDomain.
 */
void gfs_variables_swap (GfsVariable * v1, GfsVariable * v2)
{
  guint i;

  g_return_if_fail (v1 != NULL);
  g_return_if_fail (v2 != NULL);
  g_return_if_fail (v1->domain == v2->domain);

  i = v1->i; v1->i = v2->i; v2->i = i;
}

/**
 * gfs_variable_set_vector:
 * @v: the components of the vector.
 * @n: the vector dimension.
 *
 * Sets @v[0],...,@v[n-1] as components of a vector quantity.
 */
void gfs_variable_set_vector (GfsVariable ** v, guint n)
{
  guint i, j;

  g_return_if_fail (v != NULL);
  g_return_if_fail (n <= FTT_DIMENSION);

  for (i = 0; i < n; i++) {
    g_return_if_fail (v[i] != NULL);
    v[i]->component = i;
    for (j = 0; j < n; j++)
      v[i]->vector[j] = v[j];
  }    
}

/* GfsVariableBoolean: object */

static void boolean_fine_coarse (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  gint i;

  ftt_cell_children (parent, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i] && GFS_VALUE (child.c[i], v) < 0.) {
      GFS_VALUE (parent, v) = -1.;
      return;
    }

  gfs_get_from_below_intensive (parent, v);
}

static void boolean_coarse_fine (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  gdouble val = GFS_VALUE (parent, v);
  gint i;

  ftt_cell_children (parent, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i])
      GFS_VALUE (child.c[i], v) = val;
}

static void variable_boolean_init (GfsVariable * v)
{
  v->fine_coarse = boolean_fine_coarse;
  v->coarse_fine = boolean_coarse_fine;
  v->description = g_strdup ("Boolean");
}

GfsVariableClass * gfs_variable_boolean_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_boolean_info = {
      "GfsVariableBoolean",
      sizeof (GfsVariable),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) variable_boolean_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_boolean_info);
  }

  return klass;
}

/* GfsVariableTracer: object */

static void variable_tracer_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_tracer_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type == '{')
    gfs_advection_params_read (&GFS_VARIABLE_TRACER (*o)->advection, fp);
  if (fp->type != GTS_ERROR && fp->type == '{')
    g_warning ("%d:%d: specifying diffusion parameters is not done here anymore!",
	       fp->line, fp->pos);
}

static void variable_tracer_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_tracer_class ())->parent_class->write) (o, fp);

  fputc (' ', fp);
  gfs_advection_params_write (&GFS_VARIABLE_TRACER (o)->advection, fp);
}

static void variable_tracer_class_init (GtsObjectClass * klass)
{
  klass->read = variable_tracer_read;
  klass->write = variable_tracer_write;
}

static void variable_tracer_init (GfsVariableTracer * v)
{
  gfs_advection_params_init (&v->advection);
  v->advection.gradient = gfs_center_van_leer_gradient;
  v->advection.flux = gfs_face_advection_flux;
  v->advection.v = GFS_VARIABLE1 (v);
  v->advection.fv = NULL;
  GFS_VARIABLE1 (v)->description = g_strdup ("Tracer");
}

GfsVariableClass * gfs_variable_tracer_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_tracer_info = {
      "GfsVariableTracer",
      sizeof (GfsVariableTracer),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_tracer_class_init,
      (GtsObjectInitFunc) variable_tracer_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_tracer_info);
  }

  return klass;
}

/* GfsVariableResidual: Object */

static void scale_residual (FttCell * cell, GfsVariable * res)
{
  gdouble size = ftt_cell_size (cell);
  gdouble dt = GFS_SIMULATION (res->domain)->advection_params.dt;
  GFS_VARIABLE (cell, res->i) *= dt*dt/(size*size);
}

static gboolean variable_residual_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_variable_class ())->event) (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) scale_residual, event);
    return TRUE;
  }
  return FALSE;
}

static void variable_residual_class_init (GfsEventClass * klass)
{
  klass->event = variable_residual_event;
}

static void variable_residual_init (GfsVariable * v)
{
  v->description = g_strdup ("Residual of the Poisson equation");
}

GfsVariableClass * gfs_variable_residual_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_residual_info = {
      "GfsVariableResidual",
      sizeof (GfsVariable),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_residual_class_init,
      (GtsObjectInitFunc) variable_residual_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_residual_info);
  }

  return klass;
}

/* GfsVariableFiltered: object */

static void variable_filtered_read (GtsObject ** o, GtsFile * fp)
{
  GfsVariableFiltered * v = GFS_VARIABLE_FILTERED (*o);
  GfsDomain * domain;

  (* GTS_OBJECT_CLASS (gfs_variable_filtered_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(v->v = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting a number (niter)");
    return;
  }
  v->niter = atoi (fp->token->str);
  if (v->niter == 0) {
    gts_file_error (fp, "niter must be strictly positive");
    return;
  }
  gts_file_next_token (fp);  

  if (GFS_VARIABLE1 (v)->description)
    g_free (GFS_VARIABLE1 (v)->description);
  GFS_VARIABLE1 (v)->description = g_strjoin (" ", "Variable", v->v->name, "filtered", NULL);

  GFS_VARIABLE1 (v)->units = v->v->units;
}

static void variable_filtered_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_filtered_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s %d", GFS_VARIABLE_FILTERED (o)->v->name, GFS_VARIABLE_FILTERED (o)->niter);
}

static void variable_filtered_event_half (GfsEvent * event, GfsSimulation * sim)
{
  guint n = GFS_VARIABLE_FILTERED (event)->niter;
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsVariable * v = GFS_VARIABLE1 (event);

  gfs_domain_filter (domain, GFS_VARIABLE_FILTERED (event)->v, v);
  while (--n)
    gfs_domain_filter (domain, v, NULL);
  gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) v->fine_coarse, v);
  gfs_domain_bc (domain, FTT_TRAVERSE_NON_LEAFS, -1, v);
}

static gboolean variable_filtered_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_filtered_class ())->parent_class)->event)
      (event, sim)) {
    variable_filtered_event_half (event, sim);
    return TRUE;
  }
  return FALSE;
}

static void variable_filtered_class_init (GtsObjectClass * klass)
{
  klass->read = variable_filtered_read;
  klass->write = variable_filtered_write;
  GFS_EVENT_CLASS (klass)->event = variable_filtered_event;
  GFS_EVENT_CLASS (klass)->event_half = variable_filtered_event_half;
}

static void variable_filtered_init (GfsEvent * v)
{
  /* the variable/event may need to be initialised at the start */
  v->start = -1;
}

GfsVariableClass * gfs_variable_filtered_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_filtered_info = {
      "GfsVariableFiltered",
      sizeof (GfsVariableFiltered),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_filtered_class_init,
      (GtsObjectInitFunc) variable_filtered_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_filtered_info);
  }

  return klass;
}

/* GfsVariableDiagonal: object */

static void unity (FttCell * cell, GfsVariable * v)
{
  GFS_VARIABLE (cell, v->i) = 1.;
}

static void variable_diagonal (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsVariable * tmp = data[1];
  GfsGradient g;
  FttCellNeighbors neighbor;
  FttCellFace f;
  GfsGradient ng;

  GFS_VALUE (cell, tmp) = G_MAXDOUBLE;
  g.a = g.b = 0.;
  f.cell = cell;
  ftt_cell_neighbors (cell, &neighbor);
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = neighbor.c[f.d];
    if (f.neighbor) {
      gfs_face_weighted_gradient (&f, &ng, tmp->i, -1);
      g.a += ng.a;
      g.b += ng.b;
    }
  }
  if (g.a > 0.)
    GFS_VALUE (cell, v) = g.b/g.a;
  else
    GFS_VALUE (cell, v) = G_MAXDOUBLE;
  GFS_VALUE (cell, tmp) = 1.;
}

static gboolean variable_diagonal_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_diagonal_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsVariable * tmp = gfs_temporary_variable (domain);
    gpointer data[2];
    data[0] = event;
    data[1] = tmp;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) unity, tmp);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, tmp);
    gfs_poisson_coefficients (domain, sim->physical_params.alpha);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) variable_diagonal, data);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, GFS_VARIABLE1 (event));
    gts_object_destroy (GTS_OBJECT (tmp));
    return TRUE;
  }
  return FALSE;
}

static void variable_diagonal_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = variable_diagonal_event;
}

GfsVariableClass * gfs_variable_diagonal_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_diagonal_info = {
      "GfsVariableDiagonal",
      sizeof (GfsVariable),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_diagonal_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_diagonal_info);
  }

  return klass;
}

/* GfsVariableFunction: Object */

static void variable_function_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_VARIABLE_FUNCTION (o)->f));

  (* GTS_OBJECT_CLASS (gfs_variable_function_class ())->parent_class->destroy) (o);
}

static void variable_function_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_function_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsVariableFunction * v = GFS_VARIABLE_FUNCTION (*o);
  gfs_function_read (v->f, GFS_DOMAIN (gfs_object_simulation (*o)), fp);
  if (fp->type == GTS_ERROR)
    return;
  gfs_function_set_units (v->f, GFS_VARIABLE1 (*o)->units);
}

static void variable_function_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_function_class ())->parent_class->write) (o, fp);

  gfs_function_write (GFS_VARIABLE_FUNCTION (o)->f, fp);
}

static void variable_function_class_init (GtsObjectClass * klass)
{
  klass->destroy = variable_function_destroy;
  klass->read = variable_function_read;
  klass->write = variable_function_write;
}

static void variable_function_coarse_fine (FttCell * parent, GfsVariable * v)
{
  GfsFunction * f = GFS_VARIABLE_FUNCTION (v)->f;
  FttCellChildren child;
  guint n;

  ftt_cell_children (parent, &child);
  for (n = 0; n < FTT_CELLS; n++)
    if (child.c[n])
      GFS_VALUE (child.c[n], v) = gfs_function_value (f, child.c[n]);
}

static void variable_function_init (GfsVariableFunction * v)
{
  GFS_VARIABLE1 (v)->coarse_fine = variable_function_coarse_fine;
  v->f = gfs_function_new (gfs_function_class (), 0.);
}

GfsVariableClass * gfs_variable_function_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_function_info = {
      "GfsVariableFunction",
      sizeof (GfsVariableFunction),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_function_class_init,
      (GtsObjectInitFunc) variable_function_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_function_info);
  }

  return klass;
}

#if FTT_2D

/* GfsVariableStreamFunction: Object */

static gdouble face_metric (FttCell * cell, FttDirection d, GfsDomain * domain)
{
  if (domain->face_metric) {
    FttCellFace f;
    f.cell = cell;
    f.d = d;
    return (* domain->face_metric) (domain, &f);
  }
  else
    return 1.;
}

static void init_mac_from_stream_function (FttCell * cell,
					   gdouble psi0, gdouble psi1, gdouble psi2, gdouble psi3,
					   gdouble h,
					   GfsDomain * domain,
					   GfsVariable ** u)
{
  GFS_STATE (cell)->f[0].un = (psi2 - psi1)/(h*face_metric (cell, 0, domain));
  GFS_STATE (cell)->f[1].un = (psi3 - psi0)/(h*face_metric (cell, 1, domain));
  GFS_STATE (cell)->f[2].un = (psi3 - psi2)/(h*face_metric (cell, 2, domain));
  GFS_STATE (cell)->f[3].un = (psi0 - psi1)/(h*face_metric (cell, 3, domain));

  GFS_VALUE (cell, u[0]) = (GFS_STATE (cell)->f[0].un + GFS_STATE (cell)->f[1].un)/2.;
  GFS_VALUE (cell, u[1]) = (GFS_STATE (cell)->f[2].un + GFS_STATE (cell)->f[3].un)/2.;
}

static void variable_stream_function_coarse_fine (FttCell * parent, GfsVariable * v)
{
  if (GFS_CELL_IS_BOUNDARY (parent))
    return;

  GfsFunction * f = GFS_VARIABLE_FUNCTION (v)->f;
  FttCellChildren child;
  ftt_cell_children (parent, &child);
  FttVector o, p;
  ftt_cell_pos (parent, &o);
  gdouble h = ftt_cell_size (parent)/2.;
  p.x = o.x - h; p.y = o.y - h;
  gdouble psi0 = gfs_function_spatial_value (f, &p);
  p.x = o.x + h; p.y = o.y - h;
  gdouble psi1 = gfs_function_spatial_value (f, &p);
  p.x = o.x + h; p.y = o.y + h;
  gdouble psi2 = gfs_function_spatial_value (f, &p);
  p.x = o.x - h; p.y = o.y + h;
  gdouble psi3 = gfs_function_spatial_value (f, &p);
  p.x = o.x; p.y = o.y - h;
  gdouble psi4 = gfs_function_spatial_value (f, &p);
  p.x = o.x + h; p.y = o.y;
  gdouble psi5 = gfs_function_spatial_value (f, &p);
  p.x = o.x; p.y = o.y + h;
  gdouble psi6 = gfs_function_spatial_value (f, &p);
  p.x = o.x - h; p.y = o.y;
  gdouble psi7 = gfs_function_spatial_value (f, &p);
  gdouble psi8 = gfs_function_spatial_value (f, &o);
  GfsVariable ** u = gfs_domain_velocity (v->domain);
  init_mac_from_stream_function (child.c[0], psi7, psi8, psi6, psi3, h, v->domain, u);
  init_mac_from_stream_function (child.c[1], psi8, psi5, psi2, psi6, h, v->domain, u);
  init_mac_from_stream_function (child.c[2], psi0, psi4, psi8, psi7, h, v->domain, u);
  init_mac_from_stream_function (child.c[3], psi4, psi1, psi5, psi8, h, v->domain, u);
  guint n;
  for (n = 0; n < FTT_CELLS; n++) {
    ftt_cell_pos (child.c[n], &o);
    GFS_VALUE (child.c[n], v) = gfs_function_spatial_value (f, &o);
  }
}

static void variable_stream_function_fine_coarse (FttCell * cell, GfsVariable * v)
{
  FttCellChildren child;
  ftt_cell_children (cell, &child);
  double s = 2.*face_metric (cell, 0, v->domain);
  GFS_STATE (cell)->f[0].un = 
    (face_metric (child.c[1], 0, v->domain)*GFS_STATE (child.c[1])->f[0].un +
     face_metric (child.c[3], 0, v->domain)*GFS_STATE (child.c[3])->f[0].un)/s;
  s = 2.*face_metric (cell, 1, v->domain);
  GFS_STATE (cell)->f[1].un = 
    (face_metric (child.c[0], 1, v->domain)*GFS_STATE (child.c[0])->f[1].un +
     face_metric (child.c[2], 1, v->domain)*GFS_STATE (child.c[2])->f[1].un)/s;
  s = 2.*face_metric (cell, 2, v->domain);
  GFS_STATE (cell)->f[2].un = 
    (face_metric (child.c[0], 2, v->domain)*GFS_STATE (child.c[0])->f[2].un +
     face_metric (child.c[1], 2, v->domain)*GFS_STATE (child.c[1])->f[2].un)/s;
  s = 2.*face_metric (cell, 3, v->domain);
  GFS_STATE (cell)->f[3].un = 
    (face_metric (child.c[3], 3, v->domain)*GFS_STATE (child.c[3])->f[3].un +
     face_metric (child.c[2], 3, v->domain)*GFS_STATE (child.c[2])->f[3].un)/s;
  GFS_VALUE (cell, v) = (GFS_VALUE (child.c[0], v) + GFS_VALUE (child.c[1], v) + 
			 GFS_VALUE (child.c[2], v) + GFS_VALUE (child.c[3], v))/4.;
}

static void init_streamfunction (FttCell * cell, GfsVariable * v)
{
  GfsFunction * f = GFS_VARIABLE_FUNCTION (v)->f;
  FttVector o, p;
  ftt_cell_pos (cell, &o);
  gdouble h = ftt_cell_size (cell)/2.;
  p.x = o.x - h; p.y = o.y - h;
  gdouble psi0 = gfs_function_spatial_value (f, &p);
  p.x = o.x + h; p.y = o.y - h;
  gdouble psi1 = gfs_function_spatial_value (f, &p);
  p.x = o.x + h; p.y = o.y + h;
  gdouble psi2 = gfs_function_spatial_value (f, &p);
  p.x = o.x - h; p.y = o.y + h;
  gdouble psi3 = gfs_function_spatial_value (f, &p);
  init_mac_from_stream_function (cell, psi0, psi1, psi2, psi3, 2.*h, 
				 v->domain, gfs_domain_velocity (v->domain));
}

static gboolean variable_stream_function_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_variable_function_class ())->event) (event, sim)) {
    gfs_domain_traverse_leaves (GFS_DOMAIN (sim), (FttCellTraverseFunc) init_streamfunction, event);
    return TRUE;
  }
  return FALSE;
}

static void variable_stream_function_class_init (GfsEventClass * klass)
{
  klass->event = variable_stream_function_event;
}

static void variable_stream_function_init (GfsVariable * v)
{
  v->units = 2.;
  v->coarse_fine = variable_stream_function_coarse_fine;
  v->fine_coarse = variable_stream_function_fine_coarse;
  gts_object_destroy (GTS_OBJECT (GFS_VARIABLE_FUNCTION (v)->f));
  GFS_VARIABLE_FUNCTION (v)->f = gfs_function_new (gfs_function_spatial_class (), 0.);
  GFS_EVENT (v)->istep = G_MAXINT/2;
}

GfsVariableClass * gfs_variable_stream_function_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_stream_function_info = {
      "GfsVariableStreamFunction",
      sizeof (GfsVariableFunction),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_stream_function_class_init,
      (GtsObjectInitFunc) variable_stream_function_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_function_class ()), 
				  &gfs_variable_stream_function_info);
  }

  return klass;
}

#endif /* FTT_2D */

/* GfsDerivedVariable: object */

static void gfs_derived_variable_destroy (GtsObject * object)
{
  g_free (GFS_DERIVED_VARIABLE (object)->name);
  g_free (GFS_DERIVED_VARIABLE (object)->description);

  (* GTS_OBJECT_CLASS (gfs_derived_variable_class ())->parent_class->destroy) (object);
}

static void gfs_derived_variable_class_init (GtsObjectClass * klass)
{
  klass->destroy = gfs_derived_variable_destroy;
}

GtsObjectClass * gfs_derived_variable_class (void)
{
  static GtsObjectClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_derived_variable_info = {
      "GfsDerivedVariable",
      sizeof (GfsDerivedVariable),
      sizeof (GtsObjectClass),
      (GtsObjectClassInitFunc) gfs_derived_variable_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()), 
				  &gfs_derived_variable_info);
  }

  return klass;
}

/**
 * gfs_derived_variable_from_name:
 * @i: a list of #GfsDerivedVariable.
 * @name: a name.
 *
 * Returns: the #GfsDerivedVariable @name of @list or %NULL.
 */
GfsDerivedVariable * gfs_derived_variable_from_name (GSList * i, const gchar * name)
{
  g_return_val_if_fail (name != NULL, NULL);

  while (i) {
    GfsDerivedVariable * v = i->data;
    if (!strcmp (v->name, name))
      return v;
    i = i->next;
  }
  return NULL;
}
