/* Gerris - The GNU Flow Solver                           (-*-C-*-)
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

#include "event.h"
#include "domain.h"

typedef struct _Spectrum Spectrum;

struct _Spectrum {
  gdouble ** r, ** i;
  guint n;
};

static Spectrum * spectrum_new (guint n)
{
  Spectrum * s = g_malloc (sizeof (Spectrum));
  guint i;

  s->r = g_malloc (sizeof (gdouble *)*n);
  s->i = g_malloc (sizeof (gdouble *)*n);
  for (i = 0; i < n; i++) {
    guint j;

    s->r[i] = g_malloc (sizeof (gdouble)*n);
    s->i[i] = g_malloc (sizeof (gdouble)*n);
    for (j = 0; j < n; j++) {
      gdouble k = sqrt ((i + 1)*(i + 1) + (j + 1)*(j + 1));

      s->r[i][j] = 2.*(rand()/(gdouble) RAND_MAX - 0.5)
	/(k*(1 + exp (4.*log (k/6.))));
      s->i[i][j] = 2.*(rand()/(gdouble) RAND_MAX - 0.5)
	/(k*(1 + exp (4.*log (k/6.))));
    }
  }
  s->n = n;
    
  return  s;
}

static void spectrum_destroy (Spectrum * s)
{
  guint i;

  g_return_if_fail (s != NULL);

  for (i = 0; i < s->n; i++) {
    g_free (s->r[i]);
    g_free (s->i[i]);
  }
  g_free (s->r);
  g_free (s->i);
  g_free (s);
}

static gdouble spectrum_value (Spectrum * s,
			       FttVector * pos)
{
  guint k1, k2;
  gdouble val = 0.;

  g_return_val_if_fail (s != NULL, 0.);
  g_return_val_if_fail (pos != NULL, 0.);

  for (k1 = 0; k1 < s->n; k1++)
    for (k2 = 0; k2 < s->n; k2++)
      val += s->r[k1][k2]*cos (2.*M_PI*((k1 + 1)*pos->x + (k2 + 1)*pos->y)) -
	     s->i[k1][k2]*sin (2.*M_PI*((k1 + 1)*pos->x + (k2 + 1)*pos->y));
  return val;
}

static void init_spectrum_streamfunction (FttCell * cell,
					  Spectrum * s)
{
  FttVector pos;

  ftt_cell_pos (cell, &pos);
  GFS_STATE (cell)->div = spectrum_value (s, &pos);
}

static void init_velocity_from_streamfunction (FttCell * cell,
					       GfsVariable * stream)
{
  gdouble size = ftt_cell_size (cell);

  GFS_STATE (cell)->u = - gfs_center_gradient (cell, FTT_Y, stream->i)/size;
  GFS_STATE (cell)->v =   gfs_center_gradient (cell, FTT_X, stream->i)/size;
}

/* GfsInitVorticitySpectrum: Header */

GfsEventClass * gfs_init_vorticity_spectrum_class  (void);

/* GfsInitVorticitySpectrum: Object */

static gboolean gfs_init_vorticity_spectrum_event (GfsEvent * event, 
						   GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_vorticity_spectrum_class ())->parent_class)->event) (event, sim)) {
    Spectrum * s = spectrum_new (128);

    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
        (FttCellTraverseFunc) init_spectrum_streamfunction, s);
    spectrum_destroy (s);
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
        (FttCellTraverseFunc) init_velocity_from_streamfunction, gfs_div);    
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_vorticity_spectrum_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_vorticity_spectrum_event;
}

GfsEventClass * gfs_init_vorticity_spectrum_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_vorticity_spectrum_info = {
      "GfsInitVorticitySpectrum",
      sizeof (GfsEvent),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_init_vorticity_spectrum_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_vorticity_spectrum_info);
  }

  return klass;
}

/* GfsInitVorticityRandom: Header */

GfsEventClass * gfs_init_vorticity_random_class  (void);

/* GfsInitVorticityRandom: Object */

static void init_random_streamfunction (FttCell * cell)
{
  GFS_STATE (cell)->div = rand()/(gdouble) RAND_MAX;
}

static gboolean gfs_init_vorticity_random_event (GfsEvent * event, 
						 GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_vorticity_random_class ())->parent_class)->event) (event, sim)) {
    srand (GFS_DOMAIN (sim)->pid);
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
        (FttCellTraverseFunc) init_random_streamfunction, NULL);
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
        (FttCellTraverseFunc) init_velocity_from_streamfunction, gfs_div);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_vorticity_random_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_vorticity_random_event;
}

GfsEventClass * gfs_init_vorticity_random_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_vorticity_random_info = {
      "GfsInitVorticityRandom",
      sizeof (GfsEvent),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_init_vorticity_random_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_vorticity_random_info);
  }

  return klass;
}

/* Initialize module */

const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_init_vorticity_spectrum_class ();
  gfs_init_vorticity_random_class ();
  return NULL;
}
