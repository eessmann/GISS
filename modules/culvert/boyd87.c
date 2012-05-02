/* Implementation of culvert model based on
 *
 * "Generalised Head-Discharge Equations for Culverts", M. J. Boyd,
 * Fourth national local government engineering conference, Perth, 17-20
 * August 1987.  
 *
 */

#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "boyd87.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/**** Section 3: critical depth ****/

static double dc_box (double B, double Q)
{
  return 0.4672*pow (Q/B, 0.667); /* equation (1) */
}

static double dc_pipe (double D, double Q, double g)
{
  double dc;
  dc = D*pow ((Q/sqrt(g)*pow (D, 2.5))/1.26, 1./3.75); /* equation (4a) */
  if (dc/D < 0.85)
    dc = D*pow ((Q/sqrt(g)*pow (D, 2.5))/0.95, 1./1.95); /* equation (4b) */
  return dc;
}

/**** Section 4.3: Generalised equations for culverts with inlet control ****/

/* Box culvert */

static double Q_inlet_box_1 (double B, double HW, double D, double g)
{
  if (HW/D < 1.35)
    /* Inlet not submerged */
    return 0.544*sqrt(g)*B*pow(HW, 1.50); /* equation (9a) */
  else
    /* Inlet submerged */
    return 0.702*sqrt(g)*B*pow(D, 0.89)*pow(HW, 0.61); /* equation (9b) */
}

/**
 * Flow rate for box culvert, inlet control.
 * HW: headwater depth
 * B: box width.
 * D: box height.
 * type: entrance type.
 * g: acceleration of gravity.
 */
double Q_inlet_box (double HW, double B, double D, int type, double g)
{
  double HW1 = HW;
  switch (type) {
  case 1: /* entrance type 1: Wingwall Flare 30 degrees to 75 degrees */
    HW1 = HW; 
    break;
  case 2: /* entrance type 2: Wingwall Flare 90 degrees and 15 degrees */
    HW1 = D*pow (HW/D/1.09, 1./0.99); /* equation (9c) */
    break;
  case 3: /* entrance type 3: Wingwall Flare 0 degrees */
    HW1 = D*pow (HW/D/1.07, 1./1.08); /* equation (9d) */
    break;
  }
  return Q_inlet_box_1 (B, HW1, D, g);
}

/* Circular pipe culvert */

static double Q_inlet_pipe_1 (double HW, double D, double g)
{
  if (HW/D < 1.2)
    /* Inlet not submerged */
    return 0.421*sqrt(g)*pow(D, 0.87)*pow(HW, 1.63); /* equation (10a) */
  else
    /* Inlet submerged */
    return 0.530*sqrt(g)*pow(D, 1.87)*pow(HW, 0.63); /* equation (10b) */
}

/**
 * Flow rate for circular pipe culvert, inlet control.
 * HW: headwater depth
 * D: diameter.
 * type: entrance type.
 * g: the acceleration of gravity.
 */
double Q_inlet_pipe (double HW, double D, int type, double g)
{
  double HW1 = HW;
  switch (type) {
  case 1: /* entrance type 1: Square edge with headwall */
    HW1 = HW; 
    break;
  case 2: /* entrance type 2: Groove end with headwall */
    HW1 = D*pow (HW/D/0.92, 1./0.90); /* equation (10c) */
    break;
  case 3: /* entrance type 3: Groove end projecting */
    HW1 = D*pow (HW/D/0.91, 1./0.94); /* equation (10d) */
    break;
  }
  return Q_inlet_pipe_1 (HW1, D, g);
}

/**** Section 5: outlet control ****/

static double Q_outlet (double HW, double TW,
			double area, double perimeter,
			double S0, double L,
			double n, double ke,
			double g)
{
  /* equations (11) and (12a) */
  double Rh = area/perimeter;
  assert (HW + S0*L >= TW);
  double v2 = (HW + S0*L - TW)/(0.5/g + n*n/pow(Rh, 1.333)*L + ke/2./g);
  return area*sqrt (v2);
}

static double Q_outlet_box_effective (double HW, double TW,
				      double B, double D,
				      double S0, double L,
				      double n, double ke,
				      double g)
{
  if (B*D > 0.)
    return Q_outlet (HW, TW, B*D, 2.*(B + D), S0, L, n, ke, g);
  return 0.;
}

/**
 * Flow rate for box culvert, outlet control.
 * HW: headwater depth
 * TW: tailwater depth.
 * B: box width.
 * D: box height.
 * S0: slope.
 * L: length.
 * n: Manning friction coefficient.
 * ke: entrance loss coefficient.
 * g: acceleration of gravity.
 */
double Q_outlet_box (double HW, double TW,
		     double B, double D,
		     double S0, double L,
		     double n, double ke,
		     double g)
{
  double Q = Q_outlet_box_effective (HW, TW, B, D, S0, L, n, ke, g);
  if (TW > D)
    return Q;
  /* outlet not submerged */
  int nmax = 10;
  double Q0;
  do {
    double dc = dc_box (B, Q);
    double h0 = MAX ((dc + D)/2., TW);
    if (h0 > D) h0 = D;
    if (h0 > HW + S0*L) h0 = HW + S0*L;
    Q0 = Q;
    /* assume that the effective culvert cross-section is B*h0 */
    Q = Q_outlet_box_effective (HW, h0, B, h0, S0, L, n, ke, g);
  } while (nmax-- && fabs (Q - Q0) > 0.0001);
  assert (nmax > 0);
  return Q;
}

static double Q_outlet_pipe_effective (double HW, double TW,
				       double D,
				       double S0, double L,
				       double n, double ke,
				       double g)
{
  if (TW > D)
    /* outlet submerged */  
    return Q_outlet (HW, TW, M_PI*D*D/4., M_PI*D, S0, L, n, ke, g);
  else if (TW > 0.) {
    /* outlet not submerged */
    double theta = acos (1. - 2.*TW/D);               /* equation (3a) */
    double B = D*sin (theta);                         /* equation (3b) */
    double area = D*D*(theta - sin (2.*theta)/2.)/4.; /* equation (3c) */
    double perimeter = B + theta*D;
    return Q_outlet (HW, TW, area, perimeter, S0, L, n, ke, g);
  }
  else
    return 0.;
}

/**
 * Flow rate for circular pipe culvert, outlet control.
 * HW: headwater depth
 * TW: tailwater depth.
 * D: diameter.
 * S0: slope.
 * L: length.
 * n: Manning friction coefficient.
 * ke: entrance loss coefficient.
 * g: acceleration of gravity.
 */
double Q_outlet_pipe (double HW, double TW,
		      double D,
		      double S0, double L,
		      double n, double ke,
		      double g)
{
  double Q = Q_outlet_pipe_effective (HW, TW, D, S0, L, n, ke, g);
  if (TW > D)
    /* outlet submerged */
    return Q;
  /* outlet not submerged */
  int nmax = 10;
  double Q0;
  do {
    double dc = dc_pipe (D, Q, g);
    double h0 = MAX ((dc + D)/2., TW);
    if (h0 > D) h0 = D;
    if (h0 > HW + S0*L) h0 = HW + S0*L;
    Q0 = Q;
    /* assume that the effective area/perimeter is controled by h0 */
    Q = Q_outlet_pipe_effective (HW, h0, D, S0, L, n, ke, g);
  } while (nmax-- && fabs (Q - Q0) > 0.0001);
  assert (nmax > 0);
  return Q;
}

/**
 * Flow rate for circular pipe culvert.
 * HW: headwater depth
 * TW: tailwater depth.
 * D: diameter.
 * type: entrance type.
 * S0: slope.
 * L: length.
 * n: Manning friction coefficient.
 * ke: entrance loss coefficient.
 * g: acceleration of gravity.
 */
double Q_pipe (double HW, double TW,
	       double D,
	       int type,
	       double S0, double L,
	       double n, double ke,
	       double g)
{
  double Q_outlet_control = Q_outlet_pipe (HW, TW, D, S0, L, n, ke, g);
  double Q_inlet_control = Q_inlet_pipe (HW, D, type, g);
  return MIN (Q_outlet_control, Q_inlet_control);
}

/**
 * Flow rate for circular pipe culvert, outlet control.
 * HW: headwater depth
 * TW: tailwater depth.
 * B: box width.
 * D: box height.
 * type: entrance type.
 * S0: slope.
 * L: length.
 * n: Manning friction coefficient.
 * ke: entrance loss coefficient.
 * g: acceleration of gravity.
 */
double Q_box (double HW, double TW,
	      double B, double D,
	      int type,
	      double S0, double L,
	      double n, double ke,
	      double g)
{
  double Q_outlet_control = Q_outlet_box (HW, TW, B, D, S0, L, n, ke, g);
  double Q_inlet_control = Q_inlet_box (HW, B, D, type, g);
  return MIN (Q_outlet_control, Q_inlet_control);
}
