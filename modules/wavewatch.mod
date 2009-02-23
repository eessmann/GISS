/* Gerris - The GNU Flow Solver                       (-*-C-*-)
 * Copyright (C) 2001-2008 National Institute of Water and Atmospheric Research
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
#include "wave.h"
#include "wavewatch/wavewatch.h"

#define DEBUG 0

/* fixme: needs to be identical to the same function in wave.c */
static double frequency (int ik)
{
  double gamma = 1.1;
  double f0 = 0.04;
  return f0*pow(gamma, ik);
}
      
static double theta (guint ith, guint ntheta)
{
  return 2.*M_PI*ith/ntheta;
}

typedef struct {
  GfsWave * wave;
  GfsVariable * ustar, * fpi, * u10, * v10;
  REAL * A;     /* Actions (NK*NTHETA) */
  REAL * CG;    /* Group velocities (NK) */
  REAL * WN;    /* Wavenumbers (NK) */
  REAL * ALPHA; /* Nondimensional 1-D spectrum (NK) */  
} SourceParams;

static void energy_to_action (FttCell * cell, SourceParams * p)
{
  guint i, j;
  for (i = 0; i < p->wave->nk; i++)
    for (j = 0; j < p->wave->ntheta; j++)
      p->A[j + i*p->wave->ntheta] =
	GFS_VALUE (cell, p->wave->F[i][j])*p->CG[i]/(2.*M_PI*frequency (i));
}

static void action_to_energy (FttCell * cell, SourceParams * p)
{
  guint i, j;
  for (i = 0; i < p->wave->nk; i++)
    for (j = 0; j < p->wave->ntheta; j++)
      GFS_VALUE (cell, p->wave->F[i][j]) = 
	p->A[j + i*p->wave->ntheta]*(2.*M_PI*frequency (i))/p->CG[i];
}

static void source (FttCell * cell, SourceParams * p)
{
  energy_to_action (cell, p);

  INTEGER IX, IY, IMOD = 1;
  REAL DEPTH = 1000.; /* fixme: depth is fixed at 1000 m for now */

  double u10 = p->u10 ? GFS_VALUE (cell, p->u10) : 0.;
  double v10 = p->v10 ? GFS_VALUE (cell, p->v10) : 0.;
  REAL U10ABS = sqrt (u10*u10 + v10*v10);
  REAL U10DIR = atan2 (v10, u10);

  REAL USTAR = GFS_VALUE (cell, p->ustar);
  REAL FPI = GFS_VALUE (cell, p->fpi);
  REAL DTG = GFS_SIMULATION (p->wave)->advection_params.dt*3600.;

  REAL USTDIR;
  REAL EMEAN, FMEAN, WMEAN, AMAX;
  REAL CD, Z0;
  REAL DTDYN, FCUT;
  
  __w3srcemd__w3srce (&IX, &IY, &IMOD, p->A, p->ALPHA, p->WN, p->CG, &DEPTH, 
		      &U10ABS, &U10DIR, &USTAR, &USTDIR,
		      &EMEAN, &FMEAN, &WMEAN, &AMAX, 
		      &FPI, &CD, &Z0, 
		      &DTDYN, &FCUT, &DTG);

#if DEBUG
  guint i, j;
  for (i = 0; i < p->wave->nk; i++) {
    for (j = 0; j < p->wave->ntheta; j++)      
      fprintf (stderr, "%g %g %g\n", frequency (i), theta (j, p->wave->ntheta), 
	       p->A[j + i*p->wave->ntheta]);
    fprintf (stderr, "\n");
  }
#endif

  action_to_energy (cell, p);
  GFS_VALUE (cell, p->ustar) = USTAR;
  GFS_VALUE (cell, p->fpi) = FPI;
}

static void wavewatch_source (GfsWave * wave)
{
  GfsDomain * domain = GFS_DOMAIN (wave);

  if (wave->nk != 25 || wave->ntheta != 24)
    g_assert_not_implemented ();

  static gboolean initialized = FALSE;
  if (!initialized) {
    __gfsw3init__gfsw3_init ();
    initialized = TRUE;
  }

  SourceParams p;
  p.wave = wave;
  p.A = g_malloc (wave->nk*wave->ntheta*sizeof (REAL));
  p.CG = g_malloc (wave->nk*sizeof (REAL));
  p.WN = g_malloc (wave->nk*sizeof (REAL));
  p.ALPHA = g_malloc (wave->nk*sizeof (REAL));
  p.ustar = gfs_variable_from_name (domain->variables, "Ustar");
  p.fpi = gfs_variable_from_name (domain->variables, "Fpi");
  p.u10 = gfs_variable_from_name (domain->variables, "U10");
  p.v10 = gfs_variable_from_name (domain->variables, "V10");

  guint i;
  for (i = 0; i < wave->nk; i++) {
    REAL omega = 2.*M_PI*frequency (i);
    p.WN[i] = omega*omega/9.81;
    p.CG[i] = 9.81/omega/2.;
  }

  gfs_domain_traverse_leaves (domain, (FttCellTraverseFunc) source, &p);

  g_free (p.A);
  g_free (p.CG);
  g_free (p.WN);
  g_free (p.ALPHA);
}

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "wavewatch";
const gchar * g_module_check_init (void);
void          gfs_module_read     (GtsFile * fp, GfsSimulation * sim);

const gchar * g_module_check_init (void)
{
  return NULL;
}

void gfs_module_read (GtsFile * fp, GfsSimulation * sim)
{
  g_return_if_fail (fp != NULL);
  g_return_if_fail (sim != NULL);

  if (!GFS_IS_WAVE (sim)) {
    gts_file_error (fp, "wavewatch module can only be used with GfsWave");
    return;
  }

  GFS_WAVE (sim)->source = wavewatch_source;
  gfs_domain_get_or_add_variable (GFS_DOMAIN (sim), "Ustar", "Friction velocity");
  gfs_domain_get_or_add_variable (GFS_DOMAIN (sim), "Fpi",   "Peak-input frequency");
}
