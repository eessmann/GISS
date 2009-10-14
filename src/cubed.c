/* Gerris - The GNU Flow Solver
 * Copyright (C) 2009 National Institute of Water and Atmospheric Research
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

#include "cubed.h"
#include <complex.h>
#include "map.h"

#define N 30

#if 0
/* Conformal mapping Taylor coefficients: from Rancic et al, 1996, Table B.1 */

static double A[N] = {
   1.47713062600964, -0.38183510510174, -0.05573058001191, -0.01895883606818, -0.00791315785221,
  -0.00486625437708, -0.00329251751279, -0.00235481488325, -0.00175870527475, -0.00135681133278,
  -0.00107459847699, -0.00086944475948, -0.00071607115121, -0.00059867100093, -0.00050699063239,
  -0.00043415191279, -0.00037541003286, -0.00032741060100, -0.00028773091482, -0.00025458777519,
  -0.00022664642371, -0.00020289261022, -0.00018254510830, -0.00016499474461, -0.00014976117168,
  -0.00013646173946, -0.00012478875823, -0.00011449267279, -0.00010536946150, -0.00009725109376
};

static double B[N] = {
  0.67698819751739, 0.11847293456554, 0.05317178134668, 0.02965810434052, 0.01912447304028,
  0.01342565621117, 0.00998873323180, 0.00774868996406, 0.00620346979888, 0.00509010874883,
  0.00425981184328, 0.00362308956077, 0.00312341468940, 0.00272360948942, 0.00239838086555,
  0.00213001905118, 0.00190581316131, 0.00171644156404, 0.00155493768255, 0.00141600715207,
  0.00129556597754, 0.00119042140226, 0.00109804711790, 0.00101642216628, 0.00094391366522,
  0.00087919021224, 0.00082115710311, 0.00076890728775, 0.00072168382969, 0.00067885087750
};

#else
/* Conformal mapping Taylor coefficients: from map_xy2xyz.m from mitgcm */

static double A[N] = {
  1.47713057321600, -0.38183513110512, -0.05573055466344, -0.01895884801823, -0.00791314396396,
  -0.00486626515498, -0.00329250387158, -0.00235482619663, -0.00175869000970, -0.00135682443774,
  -0.00107458043205, -0.00086946107050, -0.00071604933286, -0.00059869243613, -0.00050696402446,
  -0.00043418115349, -0.00037537743098, -0.00032745130951, -0.00028769063795, -0.00025464473946,
  -0.00022659577923, -0.00020297175587, -0.00018247947703, -0.00016510295548, -0.00014967258633,
  -0.00013660647356, -0.00012466390509, -0.00011468147908, -0.00010518717478, -0.00009749136078
};

static double B[N] = {
  0.67698822171341, 0.11847295533659, 0.05317179075349, 0.02965811274764, 0.01912447871071,
  0.01342566129383, 0.00998873721022, 0.00774869352561, 0.00620347278164, 0.00509011141874,
  0.00425981415542, 0.00362309163280, 0.00312341651697, 0.00272361113245, 0.00239838233411,
  0.00213002038153, 0.00190581436893, 0.00171644267546, 0.00155493871562, 0.00141600812949,
  0.00129556691848, 0.00119042232809, 0.00109804804853, 0.00101642312253, 0.00094391466713,
  0.00087919127990, 0.00082115825576, 0.00076890854394, 0.00072168520663, 0.00067885239089
};
#endif

static complex WofZ (complex Z)
{
  complex W = 0.;
  int n = N;
  while (n-- > 0)
    W = (W + A[n])*Z;
  return W;
}

static complex ZofW (complex W)
{
  complex Z = 0.;
  int n = N;
  while (n-- > 0)
    Z = (Z + B[n])*W;
  return Z;
}

/* I^(1/3) */
#define I3 (0.86602540378444 + I/2.)
/* sqrt (3.) - 1. */
#define RA 0.73205080756888

/* Conformal mapping of a cube onto a sphere. Maps (x,y) on the
 * north-pole face of a cube to (X,Y,Z) coordinates in physical space.
 *
 * Based on f77 code from Jim Purser & Misha Rancic.
 *
 * Face is oriented normal to Z-axis with X and Y increasing with x
 * and y.
 *
 * valid ranges:  -1 < x < 1   -1 < y < 1.
 */
static void cmap_xy2XYZ (double x, double y, double * X, double * Y, double * Z)
{
  int kx = x < 0., ky = y < 0.;
  x = fabs (x); y = fabs (y);
  int kxy = y > x;

  if (kxy) {
    double tmp = x;
    x = 1. - y;
    y = 1. - tmp;
  }
  else {
    x = 1. - x;
    y = 1. - y;
  }

  complex z = (x + I*y)/2.;
  complex W;
  if (cabs (z) > 0.) {
    z = z*z*z*z;
    W = WofZ (z);
    W = I3*cpow (W*I, 1./3.);
  }
  else
    W = 0.;
  complex cb = I - 1.;
  complex cc = RA*cb/2.;
  W = (W - RA)/(cb + cc*W);
  *X = creal (W);
  *Y = cimag (W);
  double H = 2./(1. + (*X)*(*X) + (*Y)*(*Y));
  *X *= H;
  *Y *= H;
  *Z = H - 1.;
  
  if (kxy) {
    double tmp = *X;
    *X = *Y;
    *Y = tmp;
  }
  if (kx)
    *X = - *X;
  if (ky)
    *Y = - *Y;
}

/* Conformal mapping of a sphere onto a cube. Maps (X,Y,Z) coordinates
 * in physical space to (x,y) on the north-pole face of a cube.
 *
 * This is the inverse transform of cmap_xy2XYZ().
 */
static void cmap_XYZ2xy (double X, double Y, double Z, double * x, double * y)
{
  int kx = X < 0., ky = Y < 0.;
  X = fabs (X); Y = fabs (Y);
  int kxy = Y > X;

  if (kxy) {
    double tmp = X;
    X = Y;
    Y = tmp;
  }

  double H = Z + 1.;
  X /= H; Y /= H;
  complex W = X + Y*I;
  complex cb = I - 1.;
  complex cc = RA*cb/2.;
  W = (W*cb + RA)/(1. - W*cc);
  W = W/I3;
  W = W*W*W;
  W /= I;
  complex z = ZofW (W);
  z = cpow (z, 1./4.)*2.;
  *x = fabs (creal (z));
  *y = fabs (cimag (z));

  if (kxy) {
    *x = 1. - *x;
    *y = 1. - *y;
  }
  else {
    double tmp = *x;
    *x = 1. - *y;
    *y = 1. - tmp;
  }
  if (kx)
    *x = - *x;
  if (ky)
    *y = - *y;
}

static double angle_between_vectors (const double v1[3], const double v2[3])
{
  return acos (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

static void plane_normal (const double p1[3], const double p2[3], double plane[3])
{
  plane[0] = p1[1]*p2[2] - p1[2]*p2[1];
  plane[1] = p1[2]*p2[0] - p1[0]*p2[2];
  plane[2] = p1[0]*p2[1] - p1[1]*p2[0];
  double mag = sqrt (plane[0]*plane[0] + plane[1]*plane[1] + plane[2]*plane[2]);
  plane[0] /= mag;
  plane[1] /= mag;
  plane[2] /= mag;
}

static double excess_of_quad (const double v1[3], const double v2[3],
			      const double v3[3], const double v4[3])
{
  double plane1[3], plane2[3], plane3[3], plane4[3];

  plane_normal (v1, v2, plane1);
  plane_normal (v2, v3, plane2);
  plane_normal (v3, v4, plane3);
  plane_normal (v4, v1, plane4);
  
  return 2.*M_PI -
    angle_between_vectors (plane2, plane1) -
    angle_between_vectors (plane3, plane2) -
    angle_between_vectors (plane4, plane3) -
    angle_between_vectors (plane1, plane4);
}

/* GfsMapCubed: Object */

static void gfs_map_cubed_read (GtsObject ** o, GtsFile * fp)
{
  /* this mapping cannot be used independently from GfsMetricCubed */
}

static void gfs_map_cubed_write (GtsObject * o, FILE * fp)
{
  /* this mapping cannot be used independently from GfsMetricCubed */
}

static void gfs_map_cubed_class_init (GfsMapClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_map_cubed_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_map_cubed_write;
}

static void map_cubed_transform (GfsMap * map, const FttVector * src, FttVector * dest)
{
  GfsSimulation * sim = gfs_object_simulation (map);
  double lon = src->x*M_PI/180., lat = src->y*M_PI/180.;
  double X = cos (lat)*sin (lon), Y = sin (lat), Z = sqrt (1. - X*X - Y*Y);
  double x, y;
  cmap_XYZ2xy (X, Y, Z, &x, &y);
  dest->x = x/2.*sim->physical_params.L;
  dest->y = y/2.*sim->physical_params.L;
  dest->z = src->z;
}

static void map_cubed_inverse (GfsMap * map, const FttVector * src, FttVector * dest)
{
  GfsSimulation * sim = gfs_object_simulation (map);
  double X, Y, Z;
  cmap_xy2XYZ (2.*src->x/sim->physical_params.L, 2.*src->y/sim->physical_params.L, &X, &Y, &Z);
  double lat = asin (Y);
  double lon = asin (X/cos (lat));
  dest->x = lon*180./M_PI;
  dest->y = lat*180./M_PI;
  dest->z = src->z;
}

static void gfs_map_cubed_init (GfsMap * map)
{
  map->transform = map_cubed_transform;
  map->inverse =   map_cubed_inverse;
}

static GfsMapClass * gfs_map_cubed_class (void)
{
  static GfsMapClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_map_cubed_info = {
      "GfsMapCubed",
      sizeof (GfsMap),
      sizeof (GfsMapClass),
      (GtsObjectClassInitFunc) gfs_map_cubed_class_init,
      (GtsObjectInitFunc) gfs_map_cubed_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_map_class ()), &gfs_map_cubed_info);
  }

  return klass;
}

/* GfsMetricCubed: Object */

static gdouble metric_cubed_face_metric (const GfsDomain * domain, const FttCellFace * face, 
					 const gpointer data)
{
  if (face->d/2 > 1)
    return 1.;
  FttVector p;
  ftt_face_pos (face, &p);
  double h = ftt_cell_size (face->cell);
  double v1[3], v2[3];
  v1[0] = v2[0] = 2.*p.x; v1[1] = v2[1] = 2.*p.y;
  FttComponent c = FTT_ORTHOGONAL_COMPONENT (face->d/2);
  v1[c] -= h; v2[c] += h;
  cmap_xy2XYZ (v1[0], v1[1], &v1[0], &v1[1], &v1[2]);
  cmap_xy2XYZ (v2[0], v2[1], &v2[0], &v2[1], &v2[2]);
  return 2.*angle_between_vectors (v1, v2)/(M_PI*h);
}

static gdouble metric_cubed_cell_metric (const GfsDomain * domain, const FttCell * cell, 
					 const gpointer data)
{
  FttVector p;
  ftt_cell_pos (cell, &p);
  double h = ftt_cell_size (cell);
  double v1[3], v2[3], v3[3], v4[3];
  v1[0] = 2.*p.x - h; v1[1] = 2.*p.y - h;
  cmap_xy2XYZ (v1[0], v1[1], &v1[0], &v1[1], &v1[2]);
  v2[0] = 2.*p.x - h; v2[1] = 2.*p.y + h;
  cmap_xy2XYZ (v2[0], v2[1], &v2[0], &v2[1], &v2[2]);
  v3[0] = 2.*p.x + h; v3[1] = 2.*p.y + h;
  cmap_xy2XYZ (v3[0], v3[1], &v3[0], &v3[1], &v3[2]);
  v4[0] = 2.*p.x + h; v4[1] = 2.*p.y - h;
  cmap_xy2XYZ (v4[0], v4[1], &v4[0], &v4[1], &v4[2]);
  return 4.*excess_of_quad (v1, v2, v3, v4)/(M_PI*M_PI*h*h);
}

static gdouble metric_cubed_solid_metric (const GfsDomain * domain, const FttCell * cell, 
					  const gpointer data)
{
  g_assert (GFS_IS_MIXED (cell));
  g_assert_not_implemented ();
  return 1.;
}

static void metric_cubed_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_metric_cubed_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (domain->metric_data || domain->face_metric || domain->cell_metric || domain->solid_metric) {
    gts_file_error (fp, "cannot use multiple metrics (yet)");
    return;
  }

  GtsObject * map = gts_object_new (GTS_OBJECT_CLASS (gfs_map_cubed_class ()));
  gfs_object_simulation_set (map, domain);
  gts_container_add (GTS_CONTAINER (GFS_SIMULATION (domain)->maps), GTS_CONTAINEE (map));

  domain->face_metric  = metric_cubed_face_metric;
  domain->cell_metric  = metric_cubed_cell_metric;
  domain->solid_metric = metric_cubed_solid_metric;
}

static void metric_cubed_class_init (GtsObjectClass * klass)
{
  klass->read = metric_cubed_read;
}

static void metric_cubed_init (GfsEvent * m)
{
  m->istep = G_MAXINT/2;
}

GfsEventClass * gfs_metric_cubed_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_metric_cubed_info = {
      "GfsMetricCubed",
      sizeof (GfsEvent),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) metric_cubed_class_init,
      (GtsObjectInitFunc) metric_cubed_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()), 
				  &gfs_metric_cubed_info);
  }

  return klass;
}
