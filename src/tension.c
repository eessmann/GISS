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

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (sigma)");
    return;
  }
  s->sigma = atof (fp->token->str);

  gts_file_next_token (fp);
}

static void gfs_source_tension_generic_write (GtsObject * o, FILE * fp)
{
  GfsSourceTensionGeneric * t = GFS_SOURCE_TENSION_GENERIC (o);
  (* GTS_OBJECT_CLASS (gfs_source_tension_generic_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s %g", t->c->name, t->sigma);
}

typedef struct {
  gdouble amin, amax;
  guint depth;
  GfsFunction * alpha;
  GfsVariable * c;
} StabilityParams;

static void min_max_alpha (FttCell * cell, StabilityParams * p)
{
  guint level = ftt_cell_level (cell);
  if (level > p->depth && 
      GFS_VARIABLE (cell, p->c->i) > 1e-3 && 
      GFS_VARIABLE (cell, p->c->i) < 1. - 1.e-3)
    p->depth = level;
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

  gfs_youngs_normal (cell, GFS_SOURCE_TENSION_GENERIC (s)->c, &n);
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

static void gfs_source_tension_css_event (GfsEvent * event, 
					  GfsSimulation * sim)
{
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
}

static gdouble gfs_source_tension_css_value (GfsSourceGeneric * s, 
					     FttCell * cell,
					     GfsVariable * v)
{
  return GFS_VARIABLE (cell, GFS_SOURCE_TENSION_CSS (s)->t[v->component]->i);
}

static void gfs_source_tension_css_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read =       gfs_source_tension_css_read;
  GFS_EVENT_CLASS (klass)->event_half =  gfs_source_tension_css_event;
  klass->centered_value =                gfs_source_tension_css_value;
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
      (GtsObjectInitFunc) NULL,
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
  if (GFS_IS_VARIABLE_CURVATURE (GFS_SOURCE_TENSION (s)->k))
    return gfs_source_tension_generic_stability (s, sim);
  else
    return G_MAXDOUBLE;
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

static void variable_curvature_read (GtsObject ** o, GtsFile * fp)
{
  GfsVariableCurvature * v = GFS_VARIABLE_CURVATURE (*o);
  GfsDomain * domain;

  (* GTS_OBJECT_CLASS (gfs_variable_curvature_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (f)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(v->f = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  
  if (GFS_VARIABLE1 (v)->description)
    g_free (GFS_VARIABLE1 (v)->description);
  GFS_VARIABLE1 (v)->description = g_strjoin (" ", "Curvature of the interface defined by tracer",
					      v->f->name, NULL);
}

static void variable_curvature_write (GtsObject * o, FILE * fp)
{
  GfsVariableCurvature * v = GFS_VARIABLE_CURVATURE (o);

  (* GTS_OBJECT_CLASS (gfs_variable_curvature_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s", v->f->name);
}

static void curvature (FttCell * cell, GfsVariable * v)
{
  GfsVariable * t = GFS_VARIABLE_CURVATURE (v)->f;
  gdouble f = GFS_VARIABLE (cell, t->i);

  if (GFS_IS_FULL (f))
    GFS_VARIABLE (cell, v->i) = G_MAXDOUBLE;
  else
    GFS_VARIABLE (cell, v->i) = gfs_height_curvature (cell, t);
}

typedef struct {
  GfsVariable * v, * tmp;
} DiffuseParms;

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

static void variable_curvature_event_half (GfsEvent * event, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);

  gfs_domain_timer_start (domain, "variable_curvature");

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) curvature, event);
  gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) GFS_VARIABLE1 (event)->fine_coarse, event);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, GFS_VARIABLE1 (event));

  DiffuseParms p;
  p.v = GFS_VARIABLE1 (event);
  p.tmp = gfs_temporary_variable (domain);
  guint n = 2;

  while (n--) {
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) diffuse, &p);
    gfs_variables_swap (p.v, p.tmp);
    gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) p.v->fine_coarse, p.v);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, p.v);
  }

  gts_object_destroy (GTS_OBJECT (p.tmp));

  gfs_domain_timer_stop (domain, "variable_curvature");
}

static gboolean variable_curvature_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_curvature_class ())->parent_class)->event)
      (event, sim)) {
    if (!GFS_VARIABLE_CURVATURE (event)->first_done) {
      variable_curvature_event_half (event, sim);
      GFS_VARIABLE_CURVATURE (event)->first_done = TRUE;
    }
    return TRUE;
  }
  return FALSE;
}

static void variable_curvature_class_init (GtsObjectClass * klass)
{
  klass->read = variable_curvature_read;
  klass->write = variable_curvature_write;
  GFS_EVENT_CLASS (klass)->event = variable_curvature_event;
}

static void curvature_coarse_fine (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  guint n;

  ftt_cell_children (parent, &child);
  for (n = 0; n < FTT_CELLS; n++) {
    GfsVariableCurvature * k = GFS_VARIABLE_CURVATURE (v);
    gdouble f = GFS_VARIABLE (child.c[n], k->f->i);

    GFS_VARIABLE (child.c[n], v->i) = GFS_VARIABLE (parent, v->i);
    if (!GFS_IS_FULL (f))
      g_assert (fabs (GFS_VARIABLE (child.c[n], v->i)) < G_MAXDOUBLE);
  }
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

static void variable_curvature_init (GfsVariable * v)
{
  v->coarse_fine = curvature_coarse_fine;
  v->fine_coarse = curvature_fine_coarse;
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
