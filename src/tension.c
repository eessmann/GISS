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

#include <math.h>
#include <stdlib.h>

#include "tension.h"
#include "vof.h"
#include "levelset.h"

/* GfsSourceTensionGeneric: Object */

static void gfs_source_tension_generic_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceTensionGeneric * s = GFS_SOURCE_TENSION_GENERIC (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  (* GTS_OBJECT_CLASS (gfs_source_tension_generic_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable (C)");
    return;
  }
  if ((s->c = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  s->sigma = gfs_read_constant (fp, domain)/pow (GFS_SIMULATION (domain)->physical_params.L, 3.);
}

static void gfs_source_tension_generic_write (GtsObject * o, FILE * fp)
{
  GfsSourceTensionGeneric * t = GFS_SOURCE_TENSION_GENERIC (o);
  (* GTS_OBJECT_CLASS (gfs_source_tension_generic_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s %g", t->c->name, 
	   t->sigma*pow (gfs_object_simulation (o)->physical_params.L, 3.));
}

typedef struct {
  gdouble amin, amax;
  guint depth;
  GfsFunction * alpha;
  GfsVariable * c;
} StabilityParams;

static void interface_level (FttCell * cell, StabilityParams * p)
{
  guint level = ftt_cell_level (cell);
  if (level > p->depth &&
      GFS_VARIABLE (cell, p->c->i) > 1e-3 && 
      GFS_VARIABLE (cell, p->c->i) < 1. - 1.e-3)
    p->depth = level;
}

static void min_max_alpha (FttCell * cell, StabilityParams * p)
{
  interface_level (cell, p);
  if (p->alpha) {
    gdouble a = gfs_function_value (p->alpha, cell);
    if (a < p->amin) p->amin = a;
    if (a > p->amax) p->amax = a;
  }
}

static gdouble gfs_source_tension_generic_stability (GfsSourceGeneric * s,
						     GfsSimulation * sim)
{
  GfsSourceTensionGeneric * t = GFS_SOURCE_TENSION_GENERIC (s);
  gdouble h;
  StabilityParams p = { G_MAXDOUBLE, -G_MAXDOUBLE, 0 };

  p.alpha = sim->physical_params.alpha;
  p.c = t->c;
  gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) min_max_alpha, &p);
  h = ftt_level_size (p.depth);
  if (p.alpha) {
    gdouble rhom = (1./p.amin + 1./p.amax)/2.;
    return sqrt (rhom*h*h*h/(M_PI*t->sigma));
  }
  else
    return sqrt (h*h*h/(M_PI*t->sigma));
}

static void gfs_source_tension_generic_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read =       gfs_source_tension_generic_read;
  GTS_OBJECT_CLASS (klass)->write =      gfs_source_tension_generic_write;
  klass->stability =                     gfs_source_tension_generic_stability;
}

GfsSourceGenericClass * gfs_source_tension_generic_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_tension_generic_info = {
      "GfsSourceTensionGeneric",
      sizeof (GfsSourceTensionGeneric),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_tension_generic_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = 
      gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
			    &gfs_source_tension_generic_info);
  }

  return klass;
}

/* GfsSourceTensionCSS: Object */

static void gfs_source_tension_css_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceTensionCSS * s = GFS_SOURCE_TENSION_CSS (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  FttComponent c;

  (* GTS_OBJECT_CLASS (gfs_source_tension_css_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  for (c = 0; c < FTT_DIMENSION; c++) {
    static gchar * name[3] = {"_Tx", "_Ty", "_Tz"};
    if ((s->t[c] = gfs_variable_from_name (domain->variables, name[c])) == NULL)
      s->t[c] = gfs_domain_add_variable (domain, name[c], NULL);
  }
}

static void foreach_cell_normal (FttCell * cell, GfsSourceTensionCSS * s)
{
  FttVector n;
  gdouble nn = 0.;
  gdouble sigh = GFS_SOURCE_TENSION_GENERIC (s)->sigma/ftt_cell_size (cell);
  FttComponent c;

  gfs_youngs_gradient (cell, GFS_SOURCE_TENSION_GENERIC (s)->c, &n);
  for (c = 0; c < FTT_DIMENSION; c++)
    nn += (&n.x)[c]*(&n.x)[c];
  nn = sqrt (nn + 1e-50);
  GFS_VARIABLE (cell, s->g[0]->i) = sigh*n.x*n.x/nn;
  GFS_VARIABLE (cell, s->g[1]->i) = sigh*n.y*n.y/nn;
  GFS_VARIABLE (cell, s->g[2]->i) = sigh*n.x*n.y/nn;
}

static void foreach_cell_tension_css (FttCell * cell, GfsSourceTensionCSS * s)
{
  gdouble h = ftt_cell_size (cell);
  FttVector nx, ny, nxy;
  GfsSimulation * sim = gfs_object_simulation (s);
  gdouble alpha = sim->physical_params.alpha ? 
    gfs_function_value (sim->physical_params.alpha, cell) : 1.;

  gfs_youngs_gradient (cell, s->g[0], &nx);
  gfs_youngs_gradient (cell, s->g[1], &ny);
  gfs_youngs_gradient (cell, s->g[2], &nxy);

  GFS_VARIABLE (cell, s->t[0]->i) = alpha*(ny.x - nxy.y)/h;
  GFS_VARIABLE (cell, s->t[1]->i) = alpha*(nx.y - nxy.x)/h;
}

static gboolean gfs_source_tension_css_event (GfsEvent * event, 
					      GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_source_tension_css_class ())->parent_class)->event)
      (event, sim)) {
    GfsSourceTensionCSS * s = GFS_SOURCE_TENSION_CSS (event);
    guint i;

#if (!FTT_2D)
    g_assert_not_implemented ();
#endif

    for (i = 0; i < 3; i++)
      s->g[i] = gfs_temporary_variable (GFS_DOMAIN (sim));

    gfs_domain_cell_traverse (GFS_DOMAIN (sim),
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) foreach_cell_normal, event);
    /* fixme: boundary conditions for normal */
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) foreach_cell_tension_css, event);
    for (i = 0; i < 3; i++)
      gts_object_destroy (GTS_OBJECT (s->g[i]));
    return TRUE;
  }
  return FALSE;
}

static gdouble gfs_source_tension_css_value (GfsSourceGeneric * s, 
					     FttCell * cell,
					     GfsVariable * v)
{
  return GFS_VARIABLE (cell, GFS_SOURCE_TENSION_CSS (s)->t[v->component]->i);
}

static void gfs_source_tension_css_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_source_tension_css_read;
  GFS_EVENT_CLASS (klass)->event = gfs_source_tension_css_event;
}

static void gfs_source_tension_css_init (GfsSourceGeneric * s)
{
  s->centered_value = gfs_source_tension_css_value;
}

GfsSourceGenericClass * gfs_source_tension_css_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_tension_css_info = {
      "GfsSourceTensionCSS",
      sizeof (GfsSourceTensionCSS),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_tension_css_class_init,
      (GtsObjectInitFunc) gfs_source_tension_css_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = 
      gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_tension_generic_class ()),
			    &gfs_source_tension_css_info);
  }

  return klass;
}

/* GfsSourceTension: Object */

static void gfs_source_tension_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceTension * s = GFS_SOURCE_TENSION (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  (* GTS_OBJECT_CLASS (gfs_source_tension_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable (Kappa)");
    return;
  }
  if ((s->k = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
}

static void gfs_source_tension_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_tension_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s", GFS_SOURCE_TENSION (o)->k->name);
}

static gdouble gfs_source_tension_stability (GfsSourceGeneric * s,
					     GfsSimulation * sim)
{
  if (GFS_IS_VARIABLE_POSITION (GFS_SOURCE_TENSION (s)->k)) {
    /* reduced gravity */
    StabilityParams p = { G_MAXDOUBLE, -G_MAXDOUBLE, 0 };
    p.c = GFS_SOURCE_TENSION_GENERIC (s)->c;
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) interface_level, &p);
    return sqrt (ftt_level_size (p.depth)/fabs (GFS_SOURCE_TENSION_GENERIC (s)->sigma));
  }
  else 
    /* surface tension */
    return gfs_source_tension_generic_stability (s, sim);
}

static void gfs_source_tension_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read =       gfs_source_tension_read;
  GTS_OBJECT_CLASS (klass)->write =      gfs_source_tension_write;
  klass->stability =                     gfs_source_tension_stability;
}

GfsSourceGenericClass * gfs_source_tension_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_tension_info = {
      "GfsSourceTension",
      sizeof (GfsSourceTension),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_tension_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = 
      gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_tension_generic_class ()),
			    &gfs_source_tension_info);
  }

  return klass;
}

/* GfsVariableCurvature: object */

static void variable_curvature_destroy (GtsObject * o)
{
  if (GFS_VARIABLE_CURVATURE (o)->kmax)
    gts_object_destroy (GTS_OBJECT (GFS_VARIABLE_CURVATURE (o)->kmax));

  (* GTS_OBJECT_CLASS (gfs_variable_curvature_class ())->parent_class->destroy) (o);
}

static void curvature_coarse_fine (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  guint n;

  ftt_cell_children (parent, &child);
  for (n = 0; n < FTT_CELLS; n++)
    if (child.c[n])
      GFS_VARIABLE (child.c[n], v->i) = GFS_VARIABLE (parent, v->i);
}

static void curvature_fine_coarse (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  gdouble val = 0., sa = 0.;
  guint i;

  ftt_cell_children (parent, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i] && GFS_VARIABLE (child.c[i], v->i) < G_MAXDOUBLE) {
      val += GFS_VARIABLE (child.c[i], v->i);
      sa += 1.;
    }
  if (sa > 0.)
    GFS_VARIABLE (parent, v->i) = val/sa;
  else
    GFS_VARIABLE (parent, v->i) = G_MAXDOUBLE;
}

static void variable_curvature_read (GtsObject ** o, GtsFile * fp)
{
  GfsVariableCurvature * v = GFS_VARIABLE_CURVATURE (*o);
  GfsDomain * domain;

  (* GTS_OBJECT_CLASS (gfs_variable_curvature_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (fraction or distance)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(v->f = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  g_free (GFS_VARIABLE1 (v)->description);
  GFS_VARIABLE1 (v)->description = NULL;
  if (GFS_IS_VARIABLE_TRACER (v->f)) {
    if (!GFS_IS_VARIABLE_TRACER_VOF (v->f)) {
       gts_file_error (fp, "variable `%s' is not a VOF tracer", fp->token->str);
       return;
    }
    GFS_VARIABLE1 (v)->description = g_strjoin (" ", 
						"Curvature of the interface defined by tracer",
						v->f->name, NULL);
    gts_file_next_token (fp);
    if (fp->type == GTS_STRING) {
      v->kmax = gfs_domain_get_or_add_variable (domain, fp->token->str, "Maximum curvature");
      if (v->kmax) {
	v->kmax->coarse_fine = curvature_coarse_fine;
	v->kmax->fine_coarse = curvature_fine_coarse;
	gts_file_next_token (fp);
      }
      else if (!GFS_IS_VARIABLE_POSITION (v)) {
	gts_file_error (fp, "`%s' is a reserved variable name", fp->token->str);
	return;
      }
    }
  }
  else if (GFS_IS_VARIABLE_DISTANCE (v->f)) {
    GFS_VARIABLE1 (v)->description = g_strjoin (" ", 
						"Curvature of the interface defined by distance",
						v->f->name, NULL); 
    gts_file_next_token (fp);
  }
  else {
    gts_file_error (fp, "variable `%s' is neither a tracer nor a distance", fp->token->str);
    return;
  }
}

static void variable_curvature_write (GtsObject * o, FILE * fp)
{
  GfsVariableCurvature * v = GFS_VARIABLE_CURVATURE (o);

  (* GTS_OBJECT_CLASS (gfs_variable_curvature_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s", v->f->name);
  if (v->kmax)
    fprintf (fp, " %s", v->kmax->name);
}

static void height_curvature (FttCell * cell, GfsVariable * v)
{
  GfsVariable * t = GFS_VARIABLE_CURVATURE (v)->f;
  GfsVariable * kmax = GFS_VARIABLE_CURVATURE (v)->kmax;
  gdouble f = GFS_VALUE (cell, t);

  if (GFS_IS_FULL (f)) {
    GFS_VALUE (cell, v) = G_MAXDOUBLE;
    if (kmax)
      GFS_VALUE (cell, kmax) = G_MAXDOUBLE;
  }
  else {
    if (kmax) {
      gdouble k;
      GFS_VALUE (cell, v) = gfs_height_curvature (cell, GFS_VARIABLE_TRACER_VOF (t), &k);
      GFS_VALUE (cell, kmax) = k;
    }
    else
      GFS_VALUE (cell, v) = gfs_height_curvature (cell, GFS_VARIABLE_TRACER_VOF (t), NULL);
  }
}

static void fit_curvature (FttCell * cell, GfsVariable * v)
{
  GfsVariable * t = GFS_VARIABLE_CURVATURE (v)->f;
  gdouble f = GFS_VALUE (cell, t);

  if (!GFS_IS_FULL (f) && GFS_VALUE (cell, v) == G_MAXDOUBLE) {
    GfsVariable * kmax = GFS_VARIABLE_CURVATURE (v)->kmax;
    if (kmax) {
      gdouble k;
      GFS_VALUE (cell, v) = gfs_fit_curvature (cell, GFS_VARIABLE_TRACER_VOF (t), &k);
      GFS_VALUE (cell, kmax) = k;
    }
    else
      GFS_VALUE (cell, v) = gfs_fit_curvature (cell, GFS_VARIABLE_TRACER_VOF (t), NULL);    
  }
}

typedef struct {
  GfsVariable * v, * f, * tmp;
} DiffuseParms;

#define FMIN 0.01

static void diffuse_kmax (FttCell * cell, DiffuseParms * p)
{
  gdouble f = GFS_VARIABLE (cell, p->f->i);
  if (GFS_VARIABLE (cell, p->v->i) < G_MAXDOUBLE && f*(1. - f) > FMIN*(1. - FMIN))
    GFS_VARIABLE (cell, p->tmp->i) = GFS_VARIABLE (cell, p->v->i);
  else {
    FttCellNeighbors neighbor;
    gdouble sa = 0., s = 0.;
    FttDirection d;

    ftt_cell_neighbors (cell, &neighbor);
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (neighbor.c[d] && GFS_VARIABLE (neighbor.c[d], p->v->i) < G_MAXDOUBLE) {
	gdouble f = GFS_VARIABLE (neighbor.c[d], p->f->i);
	if (f*(1. - f) > FMIN*(1. - FMIN)) {
	  f *= 1. - f;
	  s += f*GFS_VARIABLE (neighbor.c[d], p->v->i);
	  sa += f;
	}
      }
    if (sa > 0.)
      GFS_VARIABLE (cell, p->tmp->i) = s/sa;
    else
      GFS_VARIABLE (cell, p->tmp->i) = GFS_VARIABLE (cell, p->v->i);
  }
}

static void diffuse (FttCell * cell, DiffuseParms * p)
{
  if (GFS_VARIABLE (cell, p->v->i) < G_MAXDOUBLE)
    GFS_VARIABLE (cell, p->tmp->i) = GFS_VARIABLE (cell, p->v->i);
  else {
    FttCellNeighbors neighbor;
    gdouble sa = 0., s = 0.;
    FttDirection d;

    ftt_cell_neighbors (cell, &neighbor);
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (neighbor.c[d] && GFS_VARIABLE (neighbor.c[d], p->v->i) < G_MAXDOUBLE) {
	s += GFS_VARIABLE (neighbor.c[d], p->v->i);
	sa += 1.;
      }
    if (sa > 0.)
      GFS_VARIABLE (cell, p->tmp->i) = s/sa;
    else
      GFS_VARIABLE (cell, p->tmp->i) = G_MAXDOUBLE;
  }
}

static void variable_curvature_diffuse (GfsVariable * v, GfsVariable * f, 
					GfsSimulation * sim, guint n)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  FttCellTraverseFunc diff = f ? (FttCellTraverseFunc) diffuse_kmax : (FttCellTraverseFunc) diffuse;
  DiffuseParms p;
  p.v = v;
  p.f = f;
  p.tmp = gfs_temporary_variable (domain);

  while (n--) {
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, diff, &p);
    gfs_variables_swap (p.v, p.tmp);
    gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) p.v->fine_coarse, p.v);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, p.v);
  }

  gts_object_destroy (GTS_OBJECT (p.tmp));  
}

static void variable_curvature_from_fraction (GfsEvent * event, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsVariable * kmax = GFS_VARIABLE_CURVATURE (event)->kmax;

  gfs_domain_timer_start (domain, "variable_curvature");

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) height_curvature, event);
  gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) GFS_VARIABLE1 (event)->fine_coarse, event);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, GFS_VARIABLE1 (event));
  if (kmax) {
    gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) kmax->fine_coarse, kmax);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, kmax);
    variable_curvature_diffuse (kmax, GFS_VARIABLE_CURVATURE (event)->f, sim, 1);
  }
  variable_curvature_diffuse (GFS_VARIABLE1 (event), NULL, sim, 1);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) fit_curvature, event);
  gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) GFS_VARIABLE1 (event)->fine_coarse, event);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, GFS_VARIABLE1 (event));
  if (kmax) {
    gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) kmax->fine_coarse, kmax);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, kmax);
    variable_curvature_diffuse (kmax, GFS_VARIABLE_CURVATURE (event)->f, sim, 1);
  }
  variable_curvature_diffuse (GFS_VARIABLE1 (event), NULL, sim, 1);

  gfs_domain_timer_stop (domain, "variable_curvature");
}

static void normal (FttCell * cell, gpointer * data)
{
  GfsVariable ** nv = data[0];
  GfsVariable * d = GFS_VARIABLE_CURVATURE (data[1])->f;
  GtsVector n = { 0., 0., 0. };
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++)
    n[c] = gfs_center_gradient (cell, c, d->i);
  gts_vector_normalize (n);
  for (c = 0; c < FTT_DIMENSION; c++)
    GFS_VARIABLE (cell, nv[c]->i) = n[c];
}

static void distance_curvature (FttCell * cell, gpointer * data)
{
  GfsVariable ** nv = data[0];
  gdouble kappa = 0.;
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++)
    kappa += gfs_center_gradient (cell, c, nv[c]->i);
  GFS_VALUE (cell, nv[FTT_DIMENSION]) = kappa/ftt_cell_size (cell);
}

static void interface_curvature (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[1];
  GfsVariableCurvature * k = GFS_VARIABLE_CURVATURE (v);
  gdouble f = GFS_VARIABLE (cell, GFS_VARIABLE_DISTANCE (k->f)->v->i);

  if (GFS_IS_FULL (f))
    GFS_VARIABLE (cell, v->i) = G_MAXDOUBLE;
  else {
    GfsVariable ** nv = data[0];
    gdouble h = ftt_cell_size (cell)/2.;
    FttCell * target = cell;
    FttComponent c;
    FttVector p;

    ftt_cell_pos (cell, &p);
    for (c = 0; c < FTT_DIMENSION; c++) {
      gdouble delta = GFS_VARIABLE (cell, k->f->i)*GFS_VARIABLE (cell, nv[c]->i);
      (&p.x)[c] -= delta;
      if (fabs (delta) > h)
	target = NULL;
    }
    if (!target)
      target = gfs_domain_locate (v->domain, p, -1);
    GFS_VARIABLE (cell, v->i) = gfs_interpolate (target, p, nv[FTT_DIMENSION]);
  }
}

static void variable_curvature_from_distance (GfsEvent * event, GfsSimulation * sim)
{
  GfsVariable * n[FTT_DIMENSION + 1];
  GfsDomain * domain = GFS_DOMAIN (sim);
  gpointer data[2];
  FttComponent c;

  gfs_domain_timer_start (domain, "variable_curvature");

  for (c = 0; c < FTT_DIMENSION + 1; c++) {
    n[c] = gfs_temporary_variable (domain);
    gfs_variable_set_vector (n[c], c);
  }
  data[0] = n;
  data[1] = event;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) normal, data);
  for (c = 0; c < FTT_DIMENSION; c++)
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, n[c]);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) distance_curvature, data);
  gfs_domain_copy_bc (domain, FTT_TRAVERSE_LEAFS, -1, 
		      GFS_VARIABLE1 (event), n[FTT_DIMENSION]);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) interface_curvature, data);
  gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) GFS_VARIABLE1 (event)->fine_coarse, event);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, GFS_VARIABLE1 (event));
  for (c = 0; c < FTT_DIMENSION + 1; c++)
    gts_object_destroy (GTS_OBJECT (n[c]));

  variable_curvature_diffuse (GFS_VARIABLE1 (event), NULL, sim, 2);

  gfs_domain_timer_stop (domain, "variable_curvature");
}

static gboolean variable_curvature_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_curvature_class ())->parent_class)->event)
      (event, sim)) {
    if (GFS_IS_VARIABLE_TRACER (GFS_VARIABLE_CURVATURE (event)->f))
      variable_curvature_from_fraction (event, sim);
    else /* distance */
      variable_curvature_from_distance (event, sim);
    return TRUE;
  }
  return FALSE;
}

static void variable_curvature_class_init (GtsObjectClass * klass)
{
  klass->destroy = variable_curvature_destroy;
  klass->read = variable_curvature_read;
  klass->write = variable_curvature_write;
  GFS_EVENT_CLASS (klass)->event = variable_curvature_event;
}

static void variable_curvature_init (GfsVariable * v)
{
  v->coarse_fine = curvature_coarse_fine;
  v->fine_coarse = curvature_fine_coarse;
  v->units = -1.;
}

GfsVariableClass * gfs_variable_curvature_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_curvature_info = {
      "GfsVariableCurvature",
      sizeof (GfsVariableCurvature),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_curvature_class_init,
      (GtsObjectInitFunc) variable_curvature_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_curvature_info);
  }

  return klass;
}

/* GfsVariablePosition: object */

static void variable_position_read (GtsObject ** o, GtsFile * fp)
{
  GfsVariablePosition * v = GFS_VARIABLE_POSITION (*o);

  (* GTS_OBJECT_CLASS (gfs_variable_position_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (component)");
    return;
  }
  if (!strcmp (fp->token->str, "x"))
    v->c = FTT_X;
  else if (!strcmp (fp->token->str, "y"))
    v->c = FTT_Y;
#if !FTT_2D
  else if (!strcmp (fp->token->str, "z"))
    v->c = FTT_Z;
#endif /* 3D */
  else {
    gts_file_error (fp, "`%s' is not a valid component", fp->token->str);
    return;
  }
  if (GFS_VARIABLE1 (v)->description)
    g_free (GFS_VARIABLE1 (v)->description);
  GFS_VARIABLE1 (v)->description = g_strjoin (" ", fp->token->str,
					      "coordinate of the interface defined by tracer",
					      GFS_VARIABLE_CURVATURE (v)->f->name, NULL);
  gts_file_next_token (fp);
  if (fp->type != '\n')
    /* fixme: mapping? */
    v->ref = gfs_read_constant (fp, gfs_object_simulation (*o));
}

static void variable_position_write (GtsObject * o, FILE * fp)
{
  GfsVariablePosition * v = GFS_VARIABLE_POSITION (o);

  (* GTS_OBJECT_CLASS (gfs_variable_position_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s", v->c == FTT_X ? "x" : v->c == FTT_Y ? "y" : "z");
  if (v->ref != 0.)
    fprintf (fp, " %g", v->ref);
}

static void position (FttCell * cell, GfsVariable * v)
{
  FttVector p;

  if (gfs_vof_center (cell, GFS_VARIABLE_TRACER_VOF (GFS_VARIABLE_CURVATURE (v)->f), &p))
    GFS_VARIABLE (cell, v->i) = (&p.x)[GFS_VARIABLE_POSITION (v)->c] - 
      GFS_VARIABLE_POSITION (v)->ref;
  else
    GFS_VARIABLE (cell, v->i) = G_MAXDOUBLE;
}

static gboolean variable_position_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    
    gfs_domain_timer_start (domain, "variable_position");
    
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) position, event);
    gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) GFS_VARIABLE1 (event)->fine_coarse, event);
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, GFS_VARIABLE1 (event));
    
    variable_curvature_diffuse (GFS_VARIABLE1 (event), NULL, sim, 2);
    
    gfs_domain_timer_stop (domain, "variable_position");
    return TRUE;
  }
  return FALSE;
}

static void variable_position_class_init (GtsObjectClass * klass)
{
  klass->read = variable_position_read;
  klass->write = variable_position_write;
  GFS_EVENT_CLASS (klass)->event = variable_position_event;
}

static void variable_position_init (GfsVariable * v)
{
  v->units = 1.;
}

GfsVariableClass * gfs_variable_position_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_position_info = {
      "GfsVariablePosition",
      sizeof (GfsVariablePosition),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_position_class_init,
      (GtsObjectInitFunc) variable_position_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_curvature_class ()), 
				  &gfs_variable_position_info);
  }

  return klass;
}
