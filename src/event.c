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
#include <sys/wait.h>
#include <unistd.h>
#include "event.h"
#include "solid.h"

static gboolean gfs_event_event (GfsEvent * event, GfsSimulation * sim)
{
  if (event->t >= event->end ||
      event->i >= event->iend ||
      sim->time.t > event->end || 
      sim->time.i > event->iend) {
    gts_object_destroy (GTS_OBJECT (event));
    return (event->realised = FALSE);
  }
  if (event->end_event) {
    if (event->n == 0 &&
	(sim->time.t >= sim->time.end ||
	 sim->time.i >= sim->time.iend)) {
      event->n = 1;
      return (event->realised = TRUE);
    }
    else
      return (event->realised = FALSE);
  }
  if (sim->time.t >= event->t) {
    if (event->istep < G_MAXINT) {
      if (event->n == 0) {
	event->i = sim->time.i + event->istep;
	event->n++;
	return (event->realised = TRUE);
      }
    }
    else {
      event->n++;
      event->t = event->start + event->n*event->step;
      return (event->realised = TRUE);
    }
  }
  if (sim->time.i >= event->i) {
    if (event->step < G_MAXDOUBLE) {
      if (event->n == 0) {
	event->start = sim->time.t;
	event->t = event->start + event->step;
	event->n = 1;
	return (event->realised = TRUE);
      }
    }
    else {
      event->n++;
      event->i += event->istep;
      return (event->realised = TRUE);
    }
  }
  return (event->realised = FALSE);
}

static void gfs_event_write (GtsObject * object, FILE * fp)
{
  GfsEvent * event = GFS_EVENT (object);

  fprintf (fp, "%s { ", object->klass->info.name);
  if (event->end_event)
    fputs ("start = end ", fp);
  else {
    if (event->start > 0. && event->start < G_MAXDOUBLE/2.)
      fprintf (fp, "start = %g ", event->start);
    if (event->step < G_MAXDOUBLE)
      fprintf (fp, "step = %g ", event->step);
    if (event->end < G_MAXDOUBLE)
      fprintf (fp, "end = %g ", event->end);
    if (event->istart > 0 && event->istart < G_MAXINT/2)
      fprintf (fp, "istart = %u ", event->istart);
    if (event->istep < G_MAXINT)
      fprintf (fp, "istep = %u ", event->istep);
    if (event->iend < G_MAXINT)
      fprintf (fp, "iend = %u ", event->iend);
  }
  fputc ('}', fp);
}

static void gfs_event_init (GfsEvent * object)
{
  object->t      = 0.;
  object->start  = 0.;
  object->end    = G_MAXDOUBLE;
  object->step   = G_MAXDOUBLE;

  object->i      = 0;
  object->istart = 0;
  object->iend   = G_MAXINT;
  object->istep  = G_MAXINT;

  object->n         = 0;
  object->end_event = FALSE;
}

static void gfs_event_read (GtsObject ** o, GtsFile * fp)
{
  GfsEvent * event = GFS_EVENT (*o);
  GtsObjectClass * klass;
  gboolean class_changed = FALSE;
  GtsFileVariable var[] = {
    {GTS_STRING, "start",  TRUE},
    {GTS_DOUBLE, "end",    TRUE},
    {GTS_DOUBLE, "step",   TRUE},
    {GTS_UINT,   "istart", TRUE},
    {GTS_UINT,   "iend",   TRUE},
    {GTS_UINT,   "istep",  TRUE},
    {GTS_NONE}
  };
  gchar * start = NULL;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsEventClass)");
    return;
  }
  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_event_class ())) {
    gts_file_error (fp, "`%s' is not a GfsEvent", fp->token->str);
    return;
  }
  if (klass != (*o)->klass) {
    *o = gts_object_new (klass);
    gts_object_destroy (GTS_OBJECT (event));
    event = GFS_EVENT (*o);
    class_changed = TRUE;
  }
  gts_file_next_token (fp);

  var[0].data = &start;
  var[1].data = &event->end;
  var[2].data = &event->step;

  var[3].data = &event->istart;
  var[4].data = &event->iend;
  var[5].data = &event->istep;
 
  gts_file_assign_variables (fp, var);

  if (fp->type == GTS_ERROR)
    return;

  if (start) {
    if (!strcmp (start, "end")) {
      event->end_event = TRUE;
      if (var[1].set)
	gts_file_variable_error (fp, var, "end", 
				 "end cannot be set for an `end' event");
      else if (var[2].set)
	gts_file_variable_error (fp, var, "step", 
				 "step cannot be set for an `end' event");
      else if (var[3].set)
	gts_file_variable_error (fp, var, "istart", 
				 "istart cannot be set for an `end' event");
      else if (var[4].set)
	gts_file_variable_error (fp, var, "iend", 
				 "iend cannot be set for an `end' event");
      else if (var[5].set)
	gts_file_variable_error (fp, var, "istep", 
				 "istep cannot be set for an `end' event");
    }
    else
      event->start = atof (start);
    g_free (start);
  }

  if (fp->type == GTS_ERROR)
    return;

  if (var[2].set && var[5].set) {
    gts_file_variable_error (fp, var, "istep", 
			     "step and istep cannot be set simultaneously");
    return;
  }

  if (var[2].set && event->step <= 0.) {
    gts_file_variable_error (fp, var, "step",
			     "step `%g' must be strictly positive", 
			     event->step);
    return;
  }
  if (!var[2].set && !var[5].set && var[1].set) {
    gts_file_error (fp, "expecting a number (step or istep)");
    return;
  }
  if (var[1].set && event->end <= event->start) {
    gts_file_variable_error (fp, var, "end",
			     "end `%g' must be larger than start `%g'", 
			     event->end, event->start);
    return;
  }
  if (event->start < 0. && var[1].set) {
    gts_file_variable_error (fp, var, "end",
			     "end cannot be specified for an `init' event");
    return;
  }
  if (event->start < 0. && var[2].set) {
    gts_file_variable_error (fp, var, "step",
			     "step cannot be specified for an `init' event");
    return;
  }
  if (var[0].set || !var[3].set)
    event->t = event->start;
  else
    event->t = event->start = G_MAXDOUBLE/2.;

  if (!var[5].set && !var[2].set && var[4].set) {
    gts_file_error (fp, "expecting a number (istep or step)");
    return;
  }
  if (var[3].set && event->iend <= event->istart) {
    gts_file_variable_error (fp, var, "iend",
			     "iend `%u' must be larger than istart `%u'", 
			     event->iend, event->istart);
    return;
  }
  if (var[3].set || !var[0].set)
    event->i = event->istart;
  else
    event->i = event->istart = G_MAXINT/2;

  if (class_changed && fp->type != '\n' && klass->read)
    (* klass->read) (o, fp);
}

static void gfs_event_class_init (GfsEventClass * klass)
{
  klass->event = gfs_event_event;

  GTS_OBJECT_CLASS (klass)->write = gfs_event_write;
  GTS_OBJECT_CLASS (klass)->read  = gfs_event_read;
}

GfsEventClass * gfs_event_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_event_info = {
      "GfsEvent",
      sizeof (GfsEvent),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_event_class_init,
      (GtsObjectInitFunc) gfs_event_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = 
      gts_object_class_new (GTS_OBJECT_CLASS (gts_slist_containee_class ()),
			    &gfs_event_info);
  }

  return klass;
}

GfsEvent * gfs_event_new (GfsEventClass * klass)
{
  GfsEvent * object;

  g_return_val_if_fail (klass != NULL, NULL);

  object = GFS_EVENT (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

/**
 * gfs_event_set:
 * @e: a #GfsEvent.
 * @start: start time.
 * @end: end time.
 * @step: time step.
 * @istart: start iteration.
 * @iend: end iteration.
 * @istep: iteration step.
 *
 * Sets the properties of event @e.
 *
 * If any of the arguments is negative, the corresponding value in @e
 * is unchanged.
 */
void gfs_event_set (GfsEvent * e, 
		    gdouble start, gdouble end, gdouble step,
		    gint istart, gint iend, gint istep)
{
  g_return_if_fail (e != NULL);
  g_return_if_fail (step < 0. || istep < 0.);
  g_return_if_fail (end < 0. || start < 0. || start <= end);
  g_return_if_fail (istep >= 0 || step >= 0. || iend < 0);
  g_return_if_fail (istart < 0 || iend < 0 || istart <= iend);
  
  if (start >= 0.) e->start = start;
  if (end >= 0.)   e->end = end;
  if (step >= 0.)  e->step = step;
  if (istart >= 0) e->istart = istart;
  if (iend >= 0)   e->iend = iend;
  if (istep >= 0)  e->istep = istep;
  if (start >= 0. || istart < 0)
    e->t = e->start;
  else
    e->t = e->start = G_MAXDOUBLE/2.;
  if (istart >= 0 || start < 0.)
    e->i = e->istart;
  else
    e->i = e->istart = G_MAXINT/2;
}

/* GfsGenericInit: Object */

static void gfs_generic_init_init (GfsEvent * event)
{
  event->start = -1;
}

GfsEventClass * gfs_generic_init_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_generic_init_info = {
      "GfsGenericInit",
      sizeof (GfsEvent),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_generic_init_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_generic_init_info);
  }

  return klass;
}

/* GfsInit: Object */

static void gfs_init_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_init_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  
  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return;
  }
  fp->scope_max++;
  gts_file_next_token (fp);
  while (fp->type != GTS_ERROR && fp->type != '}') {
    if (fp->type == '\n') {
      gts_file_next_token (fp);
      continue;
    }
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a variable name");
      return;
    }
    else {
      GfsVariable * v = gfs_variable_from_name (GFS_DOMAIN (gfs_object_simulation (*o))->variables, fp->token->str);
      GfsFunction * f;

      if (!v) {
	gts_file_error (fp, "unknown variable `%s'\n", fp->token->str);
	return;
      }
      gts_file_next_token (fp);

      if (fp->type != '=') {
	 gts_file_error (fp, "expecting `=`");
	 return;
      }
      gts_file_next_token (fp);

      f = gfs_function_new (gfs_function_class (), 0.);
      gfs_function_read (f, fp);
      if (fp->type == GTS_ERROR) {
	gts_object_destroy (GTS_OBJECT (f));
	return;
      }
      g_hash_table_insert (GFS_INIT (*o)->f, v, f);
    }
  }
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);
}

static void write_f (GfsVariable * v, GfsFunction * f, FILE * fp)
{
  fprintf (fp, "  %s =", v->name);
  gfs_function_write (f, fp);
  fputc ('\n', fp);
}

static void gfs_init_write (GtsObject * o, FILE * fp)
{
  GfsInit * i = GFS_INIT (o);

  if (GTS_OBJECT_CLASS (gfs_init_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_init_class ())->parent_class->write) 
      (o, fp);

  fputs (" {\n", fp);
  g_hash_table_foreach (i->f, (GHFunc) write_f, fp);
  fputc ('}', fp);
}

static void destroy_f (GfsVariable * v, GtsObject * o)
{
  gts_object_destroy (o);
}

static void gfs_init_destroy (GtsObject * object)
{
  GfsInit * i = GFS_INIT (object);

  g_hash_table_foreach (i->f, (GHFunc) destroy_f, NULL);
  g_hash_table_destroy (i->f);

  (* GTS_OBJECT_CLASS (gfs_init_class ())->parent_class->destroy) 
    (object);
}

static void init_vf (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];
  GfsFunction * f = data[1];
  GfsSimulation * sim = data[2];
  FttVector p;

  if (v->centered)
    ftt_cell_pos (cell, &p);
  else
    gfs_cell_cm (cell, &p);
  GFS_VARIABLE (cell, v->i) = gfs_function_value (f, &p, sim->time.t);
}

static void init_f (GfsVariable * v, GfsFunction * f, GfsDomain * domain)
{
  gpointer data[3];
  
  data[0] = v;
  data[1] = f;
  data[2] = domain;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) init_vf, data);
}

static gboolean gfs_init_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_class ())->parent_class)->event) (event, sim)) {
    g_hash_table_foreach (GFS_INIT (event)->f, (GHFunc) init_f, sim);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_class_init (GfsGenericInitClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_init_destroy;
}

static void gfs_init_init (GfsInit * object)
{
  object->f = g_hash_table_new (NULL, NULL);
}

GfsGenericInitClass * gfs_init_class (void)
{
  static GfsGenericInitClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_info = {
      "GfsInit",
      sizeof (GfsInit),
      sizeof (GfsGenericInitClass),
      (GtsObjectClassInitFunc) gfs_init_class_init,
      (GtsObjectInitFunc) gfs_init_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_info);
  }

  return klass;
}

/* GfsInitFlowConstant: Object: fixme: deprecated */

static void gfs_init_flow_constant_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_init_flow_constant_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_flow_constant_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  g_warning ("GfsInitFlowConstant is deprecated you should use GfsInit instead");
}

static void gfs_init_flow_constant_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_init_flow_constant_read;
}

GfsEventClass * gfs_init_flow_constant_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_flow_constant_info = {
      "GfsInitFlowConstant",
      sizeof (GfsInit),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_init_flow_constant_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_init_class ()),
				  &gfs_init_flow_constant_info);
  }

  return klass;
}

#if FTT_2D

/* GfsInitVorticity: Object */

static void sum_volume (FttCell * cell, GtsRange * vol) 
{
  gdouble size = ftt_cell_size (cell);
  
  if (GFS_IS_MIXED (cell))
    gts_range_add_value (vol, size*size*GFS_STATE (cell)->solid->a);
  else
    gts_range_add_value (vol, size*size);
}

static void add_ddiv (FttCell * cell, gdouble * ddiv) 
{
  gdouble size = ftt_cell_size (cell);
  
  if (GFS_IS_MIXED (cell))
    GFS_STATE (cell)->div += size*size*GFS_STATE (cell)->solid->a*(*ddiv);
  else
    GFS_STATE (cell)->div += size*size*(*ddiv);
}

static void correct_div (GfsDomain * domain)
{
  GtsRange div, vol;
  gdouble ddiv;

  div = gfs_domain_stats_variable (domain, gfs_div, FTT_TRAVERSE_LEAFS, -1);
  gts_range_init (&vol);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) sum_volume, &vol);
  gts_range_update (&vol);
  ddiv = - div.mean/vol.mean;

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) add_ddiv, &ddiv);
}

static void multiply (FttCell * cell, GfsVariable * v)
{
  gdouble size = ftt_cell_size (cell);
  
  GFS_STATE (cell)->div *= size*size;
}

static void stream_from_vorticity (GfsDomain * domain,
				   GfsVariable * stream,
				   GfsVariable * vorticity,
				   gdouble tolerance)
{
  GfsNorm norm;
  guint maxlevel, maxit = 100;

  g_return_if_fail (domain != NULL);

  gfs_poisson_coefficients (domain, NULL, 1.);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) multiply, vorticity);
  correct_div (domain); /* enforce solvability condition */
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, 
			    stream);
  gfs_residual (domain, FTT_DIMENSION, FTT_TRAVERSE_LEAFS, -1, stream, vorticity, gfs_res);
  norm = gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, 1.);
  maxlevel = gfs_domain_depth (domain);
  while (norm.infty > tolerance && maxit) {
    gfs_poisson_cycle (domain, FTT_DIMENSION, 0, maxlevel, 4, stream, vorticity);
    norm = gfs_domain_norm_residual (domain, FTT_TRAVERSE_LEAFS, -1, 1.);
    maxit--;
  }
  if (maxit == 0)
    g_warning ("GfsInitVorticity: cannot solve streamfunction from vorticity\n"
	       "  (residual: %g)", norm.infty);
}

static void init_from_streamfunction (FttCell * cell)
{
  gdouble size = ftt_cell_size (cell);

  GFS_STATE (cell)->u = - gfs_center_gradient (cell, FTT_Y, GFS_GX)/size;
  GFS_STATE (cell)->v = gfs_center_gradient (cell, FTT_X, GFS_GX)/size;
}

static gboolean gfs_init_vorticity_event (GfsEvent * event, 
					  GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_vorticity_class ())->parent_class)->event) (event, sim)) {
    stream_from_vorticity (GFS_DOMAIN (sim), gfs_gx, gfs_div, 1e-9);
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_from_streamfunction,
			      NULL);
    return TRUE;
  }
  return FALSE;
}


static void gfs_init_vorticity_class_init (GfsInitVorticityClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_vorticity_event;
}

GfsInitVorticityClass * gfs_init_vorticity_class (void)
{
  static GfsInitVorticityClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_vorticity_info = {
      "GfsInitVorticity",
      sizeof (GfsInitVorticity),
      sizeof (GfsInitVorticityClass),
      (GtsObjectClassInitFunc) gfs_init_vorticity_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_vorticity_info);
  }

  return klass;
}

#endif /* FTT_2D */

/* GfsEventSum: Object */

static void gfs_event_sum_write (GtsObject * o, FILE * fp)
{
  GfsEventSum * s = GFS_EVENT_SUM (o);

  if (GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class->write)
      (o, fp);

  fprintf (fp, " %s %s", s->v->name, s->sv->name);
}

static void gfs_event_sum_read (GtsObject ** o, GtsFile * fp)
{
  GfsEventSum * s = GFS_EVENT_SUM (*o);
  GfsDomain * domain =  GFS_DOMAIN (gfs_object_simulation (s));

  if (GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  if (!(s->v = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (sv)");
    return;
  }
  if (!(s->sv = gfs_variable_from_name (domain->variables, fp->token->str))) {
    s->sv = gfs_domain_add_variable (domain, fp->token->str);
    g_assert (s->sv);
  }
  s->sv->fine_coarse = s->v->fine_coarse;
  gts_file_next_token (fp);
}

static gboolean gfs_event_sum_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_sum_class ())->parent_class)->event) (event, sim)) {
    GfsEventSum * s = GFS_EVENT_SUM (event);

    if (s->last < 0.)
      gfs_domain_cell_traverse (GFS_DOMAIN (sim),
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) gfs_cell_reset, s->sv);
    else {
      s->dt = sim->time.t - s->last;
      gfs_domain_cell_traverse (GFS_DOMAIN (sim),
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				s->sum, s);
    }
    s->last = sim->time.t;
    return TRUE;
  }
  return FALSE;
}

static void gfs_event_sum_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_event_sum_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_event_sum_write;
  GFS_EVENT_CLASS (klass)->event = gfs_event_sum_event;
}

static void sum (FttCell * cell, GfsEventSum * s)
{
  GFS_VARIABLE (cell, s->sv->i) += s->dt*GFS_VARIABLE (cell, s->v->i);
}

static void gfs_event_sum_init (GfsEventSum * object)
{
  object->last = -1.;
  object->sum = (FttCellTraverseFunc) sum;
}

GfsEventClass * gfs_event_sum_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_event_sum_info = {
      "GfsEventSum",
      sizeof (GfsEventSum),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_event_sum_class_init,
      (GtsObjectInitFunc) gfs_event_sum_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_event_sum_info);
  }

  return klass;
}

/* GfsEventSum2: Object */

static void sum2 (FttCell * cell, GfsEventSum * s)
{
  gdouble val = GFS_VARIABLE (cell, s->v->i);

  GFS_VARIABLE (cell, s->sv->i) += s->dt*val*val;
}

static void gfs_event_sum2_init (GfsEventSum * object)
{
  object->sum = (FttCellTraverseFunc) sum2;
}

GfsEventClass * gfs_event_sum2_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_event_sum2_info = {
      "GfsEventSum2",
      sizeof (GfsEventSum),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_event_sum2_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_sum_class ()),
				  &gfs_event_sum2_info);
  }

  return klass;
}

/* GfsEventStop: Object */

static void gfs_event_stop_write (GtsObject * o, FILE * fp)
{
  GfsEventStop * s = GFS_EVENT_STOP (o);

  if (GTS_OBJECT_CLASS (gfs_event_stop_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_event_stop_class ())->parent_class->write)
      (o, fp);

  fprintf (fp, " %s %g", s->v->name, s->max);
}

static void gfs_event_stop_read (GtsObject ** o, GtsFile * fp)
{
  GfsEventStop * s = GFS_EVENT_STOP (*o);
  GfsDomain * domain =  GFS_DOMAIN (gfs_object_simulation (s));

  if (GTS_OBJECT_CLASS (gfs_event_stop_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_event_stop_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  if (!(s->v = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (max)");
    return;
  }
  s->max = atof (fp->token->str);
  gts_file_next_token (fp);

  s->oldv = gfs_domain_add_variable (domain, NULL);
  s->oldv->fine_coarse = s->v->fine_coarse;
  gts_file_next_token (fp);
}

static void diff (FttCell * cell, GfsEventStop * s)
{
  GFS_VARIABLE (cell, s->oldv->i) -= GFS_VARIABLE (cell, s->v->i);
}

static void copy (FttCell * cell, GfsEventStop * s)
{
  GFS_VARIABLE (cell, s->oldv->i) = GFS_VARIABLE (cell, s->v->i);
}

static gboolean gfs_event_stop_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_stop_class ())->parent_class)->event) (event, sim)) {
    GfsEventStop * s = GFS_EVENT_STOP (event);

    if (s->last >= 0.) {
      GfsNorm n;

      gfs_domain_cell_traverse (GFS_DOMAIN (sim),
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) diff, s);
      n = gfs_domain_norm_variable (GFS_DOMAIN (sim), s->oldv,
				    FTT_TRAVERSE_LEAFS, -1);
      if (n.infty <= s->max)
	sim->time.end = sim->time.t;
    }
    gfs_domain_cell_traverse (GFS_DOMAIN (sim),
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) copy, s);
    s->last = sim->time.t;
    return TRUE;
  }
  return FALSE;
}

static void gfs_event_stop_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_event_stop_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_event_stop_write;
  GFS_EVENT_CLASS (klass)->event = gfs_event_stop_event;
}

static void gfs_event_stop_init (GfsEventStop * object)
{
  object->last = -1.;
}

GfsEventClass * gfs_event_stop_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_event_stop_info = {
      "GfsEventStop",
      sizeof (GfsEventStop),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_event_stop_class_init,
      (GtsObjectInitFunc) gfs_event_stop_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_event_stop_info);
  }

  return klass;
}

/* GfsEventScript: Object */

static void gfs_event_script_destroy (GtsObject * o)
{
  GfsEventScript * s = GFS_EVENT_SCRIPT (o);

  if (s->script) g_string_free (s->script, TRUE);

  (* GTS_OBJECT_CLASS (gfs_event_script_class ())->parent_class->destroy) (o);
} 

static void gfs_event_script_write (GtsObject * o, FILE * fp)
{
  GfsEventScript * s = GFS_EVENT_SCRIPT (o);

  if (GTS_OBJECT_CLASS (gfs_event_script_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_event_script_class ())->parent_class->write)
      (o, fp);

  fputs (" {", fp);
  if (s->script)
    fputs (s->script->str, fp);
  fputc ('}', fp);
}

static void gfs_event_script_read (GtsObject ** o, GtsFile * fp)
{
  GfsEventScript * s = GFS_EVENT_SCRIPT (*o);
  guint scope;
  gint c;

  if (GTS_OBJECT_CLASS (gfs_event_script_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_event_script_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return;
  }
  if (s->script)
    g_string_free (s->script, TRUE);
  s->script = g_string_new ("");
  scope = fp->scope_max;
  c = gts_file_getc (fp);
  while (c != EOF && fp->scope > scope) {
    g_string_append_c (s->script, c);
    c = gts_file_getc (fp);
  }
  if (fp->scope != scope) {
    gts_file_error (fp, "parse error");
    return;
  }
  gts_file_next_token (fp);
}

static gboolean gfs_event_script_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_event_script_class ())->parent_class)->event) 
      (event, sim)) {
    GfsEventScript * s = GFS_EVENT_SCRIPT (event);
    if (s->script) {
      gchar * scommand;
      gchar sname[] = "/tmp/gfsXXXXXX";
      gchar ename[] = "/tmp/gfsXXXXXX";
      gint sf = mkstemp (sname);
      gint ef = mkstemp (ename);
      gint status;
      FILE * f;

      if (sf < 0 || ef < 0) {
	g_warning ("GfsEventScript cannot create temporary files");
	return TRUE;
      }
      f = fdopen (sf, "w");
      fputs (s->script->str, f);
      fclose (f);
      scommand = g_strdup_printf ("GfsTime=%g GfsIter=%d GfsPid=%d "
				  "GFS_STOP=%d sh %s",
				  sim->time.t, sim->time.i, 
				  GFS_DOMAIN (sim)->pid,
				  GFS_EVENT_SCRIPT_STOP,
				  sname);
      fflush (stdout);
      fflush (stderr);
      status = system (scommand);
      if (status != -1)
	status = WEXITSTATUS (status);
      g_free (scommand);
      remove (sname);
      if (status == GFS_EVENT_SCRIPT_STOP)
	sim->time.end = sim->time.t;
      else if (status == -1 || status != 0) {
	FILE * ferr = fdopen (ef, "r");
	gint c;

	fputs ("Error while executing GfsEventScript:\n", stderr);
	while ((c = fgetc (ferr)) != EOF)
	  fputc (c, stderr);
	fclose (ferr);
      }
      remove (ename);
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_event_script_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gfs_event_script_destroy;
  GTS_OBJECT_CLASS (klass)->read = gfs_event_script_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_event_script_write;
  GFS_EVENT_CLASS (klass)->event = gfs_event_script_event;
}

GfsEventClass * gfs_event_script_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_event_script_info = {
      "GfsEventScript",
      sizeof (GfsEventScript),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_event_script_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_event_script_info);
  }

  return klass;
}

/* GfsInitFraction: Object */

static void gfs_init_fraction_destroy (GtsObject * object)
{
  GfsInitFraction * init = GFS_INIT_FRACTION (object);

  if (init->surface)
    gts_object_destroy (GTS_OBJECT (init->surface));

  (* GTS_OBJECT_CLASS (gfs_init_fraction_class ())->parent_class->destroy) 
    (object);
}

static void gfs_init_fraction_read (GtsObject ** o, GtsFile * fp)
{
  GfsInitFraction * init;
  GfsDomain * domain;

  if (GTS_OBJECT_CLASS (gfs_init_fraction_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_fraction_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  init = GFS_INIT_FRACTION (*o);
  domain =  GFS_DOMAIN (gfs_object_simulation (init));
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (variable)");
    return;
  }
  if ((init->c = gfs_variable_from_name (domain->variables, fp->token->str))
      == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != '{') {
    FILE * f;
    GtsFile * gf;

    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a string (filename)\n");
      return;
    }
    f = fopen (fp->token->str, "rt");
    if (f == NULL) {
      gts_file_error (fp, "cannot open file `%s'\n", fp->token->str);
      return;
    }
    gf = gts_file_new (f);
    if (gts_surface_read (init->surface, gf)) {
      gts_file_error (fp, 
		      "file `%s' is not a valid GTS file\n"
		      "%s:%d:%d: %s",
		      fp->token->str, fp->token->str,
		      gf->line, gf->pos, gf->error);
      gts_file_destroy (gf);
      fclose (f);
      return;
    }
    gts_file_destroy (gf);
    fclose (f);
  }
  else { /* embedded GTS file */
    fp->scope_max++;
    gts_file_next_token (fp);
    if (gts_surface_read (init->surface, fp))
      return;
    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    fp->scope_max--;
  }
  
  if (!gts_surface_is_orientable (init->surface)) {
    gts_file_error (fp, "surface is not orientable");
    return;
  }
  
  gts_file_next_token (fp);
}

static void gfs_init_fraction_write (GtsObject * o, FILE * fp)
{
  GfsInitFraction * init = GFS_INIT_FRACTION (o);

  if (GTS_OBJECT_CLASS (gfs_init_fraction_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_init_fraction_class ())->parent_class->write) 
      (o, fp);
  fprintf (fp, " %s { ", init->c->name);
  gts_surface_write (init->surface, fp);
  fputs ("}\n", fp);
}

static gboolean gfs_init_fraction_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_fraction_class ())->parent_class)->event) 
      (event, sim)) {
    gfs_domain_init_fraction (GFS_DOMAIN (sim), 
			      GFS_INIT_FRACTION (event)->surface,
			      GFS_INIT_FRACTION (event)->c);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_fraction_class_init (GfsInitFractionClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_fraction_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_fraction_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_fraction_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_init_fraction_destroy;
}

static void gfs_init_fraction_init (GfsInitFraction * object)
{
  object->surface = gts_surface_new (gts_surface_class (),
				     gts_face_class (),
				     gts_edge_class (),
				     gts_vertex_class ());
}

GfsInitFractionClass * gfs_init_fraction_class (void)
{
  static GfsInitFractionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_fraction_info = {
      "GfsInitFraction",
      sizeof (GfsInitFraction),
      sizeof (GfsInitFractionClass),
      (GtsObjectClassInitFunc) gfs_init_fraction_class_init,
      (GtsObjectInitFunc) gfs_init_fraction_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_fraction_info);
  }

  return klass;
}

/* GfsRemoveDroplets: Object */

static gboolean gfs_remove_droplets_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_remove_droplets_class ())->parent_class)->event) (event, sim)) {
    gfs_domain_remove_droplets (GFS_DOMAIN (sim), 
				GFS_REMOVE_DROPLETS (event)->c, 
				GFS_REMOVE_DROPLETS (event)->min);
    return TRUE;
  }
  return FALSE;
}

static void gfs_remove_droplets_read (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain;

  if (GTS_OBJECT_CLASS (gfs_remove_droplets_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_remove_droplets_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (variable)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if ((GFS_REMOVE_DROPLETS (*o)->c = gfs_variable_from_name (domain->variables, fp->token->str))
      == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (min)");
    return;
  }
  GFS_REMOVE_DROPLETS (*o)->min = atoi (fp->token->str);
  gts_file_next_token (fp);
}

static void gfs_remove_droplets_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_remove_droplets_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_remove_droplets_class ())->parent_class->write) 
      (o, fp);
  fprintf (fp, " %s %d", GFS_REMOVE_DROPLETS (o)->c->name, GFS_REMOVE_DROPLETS (o)->min);
}

static void gfs_remove_droplets_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_remove_droplets_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_remove_droplets_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_remove_droplets_write;
}

GfsEventClass * gfs_remove_droplets_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_remove_droplets_info = {
      "GfsRemoveDroplets",
      sizeof (GfsRemoveDroplets),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_remove_droplets_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_remove_droplets_info);
  }

  return klass;
}
