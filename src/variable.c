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

#include "variable.h"

/* GfsVariable: Object */

static void gfs_variable_read (GtsObject ** o, GtsFile * fp)
{
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsVariableClass)");
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (name)");
    return;
  }
  GFS_VARIABLE1 (*o)->name = g_strdup (fp->token->str);
  gts_file_next_token (fp);
}

static void gfs_variable_write (GtsObject * o, FILE * fp)
{
  fprintf (fp, "%s %s", o->klass->info.name, GFS_VARIABLE1 (o)->name);
}

static void gfs_variable_destroy (GtsObject * object)
{
  GfsVariable * v = GFS_VARIABLE1 (object);

  g_free (v->name);
  if (v->sources)
    gts_object_destroy (GTS_OBJECT (v->sources));
  if (v->surface_bc)
    gts_object_destroy (GTS_OBJECT (v->surface_bc));

  (* GTS_OBJECT_CLASS (gfs_variable_class ())->parent_class->destroy) 
    (object);
}

static void gfs_variable_clone (GtsObject * clone, GtsObject * object)
{
  GfsVariable * c = GFS_VARIABLE1 (clone);
  GfsVariable * v = GFS_VARIABLE1 (object);

  (* GTS_OBJECT_CLASS (gfs_variable_class ())->parent_class->clone) (clone, object);
  if (v->name)
    c->name = g_strdup (v->name);
  c->sources = NULL;
  c->surface_bc = NULL;
}

static void gfs_variable_class_init (GfsVariableClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_variable_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_variable_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_variable_destroy;
  GTS_OBJECT_CLASS (klass)->clone = gfs_variable_clone;
}

static void gfs_variable_init (GfsVariable * v)
{
  v->fine_coarse = (GfsVariableFineCoarseFunc) gfs_get_from_below_intensive;
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
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()), &gfs_variable_info);
  }

  return klass;
}

/**
 * gfs_variable_new:
 * @klass: a #GfsVariableClass.
 * @parent: the parent or %NULL.
 * @name: the name of the variable.
 * @centered: is the variable cell-centered?
 * @i: the variable index.
 *
 * Returns: a newly allocated #GfsVariable,
 */
GfsVariable * gfs_variable_new (GfsVariableClass * klass,
				GtsObject * parent,
				const gchar * name,
				gboolean centered,
				guint i)
{
  GfsVariable * v;

  v = GFS_VARIABLE1 (gts_object_new (GTS_OBJECT_CLASS (klass)));
  if (name)
    v->name = g_strdup (name);
  v->i = i;
  v->centered = centered;
  v->p = parent;
  v->permanent = v;

  return v;
}

/**
 * gfs_variable_list_copy:
 * @v: a #GfsVariable.
 * @parent: the parent of the new list or %NULL.
 *
 * Returns: a new variable list copy of @v.
 */
GfsVariable * gfs_variable_list_copy (GfsVariable * v,
				      GtsObject * parent)
{
  GfsVariable * start = NULL, * prev = NULL;

  while (v) {
    GfsVariable * n = GFS_VARIABLE1 (gts_object_clone (GTS_OBJECT (v)));

    n->p = parent;
    if (prev == NULL)
      start = n;
    else
      prev->next = n;
    prev = n;
    v = v->next;
  }
  return start;
}

/**
 * gfs_variable_list_destroy:
 * @v: a #GfsVariable.
 *
 * Free all the memory allocated for the list starting at @v.
 */
void gfs_variable_list_destroy (GfsVariable * v)
{
  while (v) {
    GfsVariable * next = v->next;

    gts_object_destroy (GTS_OBJECT (v));
    v = next;
  }
}

/**
 * gfs_variable_from_name:
 * @variables: the list of available #GfsVariable.
 * @name: the name of the variable to find.
 *
 * Returns: the #GfsVariable @name or %NULL if this variable name does
 * not exist.  
 */
GfsVariable * gfs_variable_from_name (GfsVariable * variables,
				      const gchar * name)
{
  g_return_val_if_fail (name != NULL, NULL);

  while (variables && (!variables->name || strcmp (name, variables->name)))
    variables = variables->next;
  return variables;
}

/**
 * gfs_variables_from_list:
 * @variables: the list of available #GfsVariable.
 * @list: a malloc'ed string containing comma separated variable names.
 * @error: where to return the variable name in case of error.
 *
 * Returns: a list of variables or %NULL in case of error, in which
 * case *@error points to the name of the unknown variable.  
 */
GfsVariable * gfs_variables_from_list (GfsVariable * variables,
				       gchar * list,
				       gchar ** error)
{
  gchar * s;
  GfsVariable * var = NULL, * prev = NULL;

  g_return_val_if_fail (list != NULL, NULL);
  g_return_val_if_fail (error != NULL, NULL);

  s = strtok (list, ",");
  while (s) {
    GfsVariable * v = gfs_variable_from_name (variables, s), * n;

    if (v == NULL) {
      *error = s;
      gfs_variable_list_destroy (var);
      return NULL;
    }
    n = gfs_variable_new (gfs_variable_class (), v->p, v->name, FALSE, v->i);
    if (prev)
      prev->next = n;
    else
      var = n;
    prev = n;
    s = strtok (NULL, ",");
  }
  return var;
}

/* GfsVariableTracer: object */

static void variable_tracer_init (GfsVariableTracer * v)
{
  gfs_advection_params_init (&v->advection);
  v->advection.gradient = gfs_center_van_leer_gradient;
  v->advection.flux = gfs_face_advection_flux;
  v->advection.v = GFS_VARIABLE1 (v);
  v->advection.fv = gfs_res;

  gfs_multilevel_params_init (&v->diffusion);
  v->diffusion.tolerance = 1e-6;
}

GfsVariableClass * gfs_variable_tracer_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_tracer_info = {
      "GfsVariableTracer",
      sizeof (GfsVariableTracer),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) variable_tracer_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_tracer_info);
  }

  return klass;
}
