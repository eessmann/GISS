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

/* GfsSourceTensionCSS: Object */

static void gfs_source_tension_css_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceTensionCSS * s = GFS_SOURCE_TENSION_CSS (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  FttComponent c;

  (* GTS_OBJECT_CLASS (gfs_source_tension_css_class ())->parent_class->read) (o, fp);
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

  for (c = 0; c < FTT_DIMENSION; c++) {
    static gchar * name[3] = {"_Tx", "_Ty", "_Tz"};
    if ((s->t[c] = gfs_variable_from_name (domain->variables, name[c])) == NULL)
      s->t[c] = gfs_domain_add_variable (domain, name[c], NULL);
  }

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (sigma)");
    return;
  }
  s->sigma = atof (fp->token->str);

  gts_file_next_token (fp);
}

static void gfs_source_tension_css_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_tension_css_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s %g", GFS_SOURCE_TENSION_CSS (o)->c->name, GFS_SOURCE_TENSION_CSS (o)->sigma);
}

static void foreach_cell_normal (FttCell * cell, GfsSourceTensionCSS * s)
{
  FttVector n;
  gdouble nn = 0.;
  gdouble sigh = s->sigma/ftt_cell_size (cell);
  FttComponent c;

  gfs_youngs_normal (cell, s->c, &n);
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

  gfs_youngs_normal (cell, s->g[0], &nx);
  gfs_youngs_normal (cell, s->g[1], &ny);
  gfs_youngs_normal (cell, s->g[2], &nxy);

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
  GTS_OBJECT_CLASS (klass)->write =      gfs_source_tension_css_write;
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
      gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
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
    gts_file_error (fp, "expecting a variable (C)");
    return;
  }
  if ((s->c = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable (Kappa)");
    return;
  }
  if ((s->k = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  if (!GFS_IS_VARIABLE_CURVATURE (s->k)) {
    gts_file_error (fp, "variable `%s' is not a VariableCurvature", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
}

static void gfs_source_tension_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_tension_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s %s", GFS_SOURCE_TENSION (o)->c->name, GFS_SOURCE_TENSION (o)->k->name);
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

static gdouble gfs_source_tension_stability (GfsSourceGeneric * s,
					     GfsSimulation * sim)
{
  GfsSourceTension * t = GFS_SOURCE_TENSION (s);
  gdouble h, sigma = 1.;
  StabilityParams p = { G_MAXDOUBLE, -G_MAXDOUBLE, 0 };

  p.alpha = sim->physical_params.alpha;
  p.c = t->c;
  gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) min_max_alpha, &p);
  h = ftt_level_size (p.depth);
  if (p.alpha) {
    gdouble rhom = (1./p.amin + 1./p.amax)/2.;
    return sqrt (rhom*h*h*h/(2.*M_PI*sigma));
  }
  else
    return sqrt (h*h*h/(2.*M_PI*sigma));
}

static void gfs_source_tension_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read =       gfs_source_tension_read;
  GTS_OBJECT_CLASS (klass)->write =      gfs_source_tension_write;
  klass->stability = gfs_source_tension_stability;
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
      gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
			    &gfs_source_tension_info);
  }

  return klass;
}
