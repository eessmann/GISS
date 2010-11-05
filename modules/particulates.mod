/* Gerris - The GNU Flow Solver			(-*-C-*-)
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

#include <stdlib.h>
#include "particle.h"
#include "source.h"

/* Adaptive 4/5 Runge-Kutta method */
/* Header */

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

typedef void      (* RKFunc)            ( gdouble t, gdouble *y, gdouble *dydt, gpointer data);

typedef struct _GfsAdaptiveRK    GfsAdaptiveRK;

struct _GfsAdaptiveRK {
  gdouble t, dtgoal, dtmax; /* independent variable t, interval to be integrated, max dt */
  gdouble * y;  /* goal function */
  gdouble * yc; /* characteristic values for nondimensionalization */
  guint n;      /* number of elements in y */ 
  gpointer data;/* user data required for f */
};

/* Object */

static void rkck (gdouble *dttry, gdouble *errmax, gdouble *yerr, 
                  gdouble *ytmp, gdouble *ak1, GfsAdaptiveRK *RKdata,
		  RKFunc derivatives)
{
  gint i;
  gdouble A2 = 0.2, A3 = 0.3, A4 = 0.6, A5 = 1.0, A6 = 0.875;
  gdouble b21 = 0.2;
  gdouble b31 = 3.0/40.0, b32 = 9.0/40.0;
  gdouble b41 = 0.3, b42 = -0.9, b43 = 1.2;
  gdouble b51 = -11.0/54.0, b52 = 2.5, b53 = -70.0/27.0, b54 = 35./27.;
  gdouble b61 = 1631.0/55296.0, b62 = 175.0/512.0, b63 = 575./13824.;
  gdouble b64 = 44275./110592., b65 = 253./4096.;
  gdouble c1 = 37./378., c3 = 250./621., c4 = 125.0/594.0, c6 = 512.0/1771.;
  gdouble dc1 = c1-2825.0/27648.0, dc3 = c3-18575.0/48384.;
  gdouble dc4 = c4 - 13525.0/55296.0, dc5 = -277.0/14336.0, dc6 = c6 - 0.25;

  gdouble *ak2, *ak3, *ak4, *ak5, *ak6;
     
  ak2 = g_malloc (RKdata->n*sizeof (gdouble));
  ak3 = g_malloc (RKdata->n*sizeof (gdouble));
  ak4 = g_malloc (RKdata->n*sizeof (gdouble));
  ak5 = g_malloc (RKdata->n*sizeof (gdouble));
  ak6 = g_malloc (RKdata->n*sizeof (gdouble));

  /* First Step */
  for (i = 0; i < RKdata->n; i++)   
    ytmp[i] = RKdata->y[i] + (*dttry)*b21*ak1[i];

  derivatives (RKdata->t + A2*(*dttry), ytmp, ak2, RKdata->data);
    
  for (i = 0; i < RKdata->n; i++) 
    ytmp[i] = RKdata->y[i] + (*dttry)*(b31*ak1[i] + b32*ak2[i]);
        
  derivatives (RKdata->t + A3*(*dttry), ytmp, ak3, RKdata->data);

  for (i = 0; i < RKdata->n; i++) 
    ytmp[i] = RKdata->y[i] + (*dttry)*(b41*ak1[i] + b42*ak2[i] + b43*ak3[i]);

  derivatives (RKdata->t + A4*(*dttry), ytmp, ak4, RKdata->data);
    
  for (i = 0; i < RKdata->n; i++) 
    ytmp[i] = RKdata->y[i] + (*dttry)*(b51*ak1[i] + b52*ak2[i] + b53*ak3[i]
				       + b54*ak4[i]);
    
  derivatives (RKdata->t + A5*(*dttry), ytmp, ak5, RKdata->data);
    
  for (i = 0; i < RKdata->n; i++) 
    ytmp[i] = RKdata->y[i] + (*dttry)*(b61*ak1[i] + b62*ak2[i] + b63*ak3[i]
				       + b64*ak4[i] + b65*ak5[i]);

  derivatives (RKdata->t + A6*(*dttry), ytmp, ak6, RKdata->data);
    
  for (i = 0; i < RKdata->n; i++) 
    ytmp[i] = RKdata->y[i] + (*dttry)*(c1*ak1[i] + c3*ak3[i]
				       + c4*ak4[i] + c6*ak6[i]);
    
  for (i = 0; i < RKdata->n; i++) 
    yerr[i] = (*dttry)*(dc1*ak1[i] + dc3*ak3[i]
			+ dc4*ak4[i] + dc5*ak5[i] + dc6*ak6[i]);
    
  g_free (ak2); g_free (ak3); g_free (ak4); g_free (ak5); g_free (ak6);
}

/* Adaptive RK method
   It integrates a function df/dt=f(t,...)
   f(t,y,dydt,data) is RKFunc whose arguments are:
   t: independent variable
   y: function
   dydt: derivative of the function
   data: extra data
   RKdata contains the info of the temporal integration 
*/
static void rkqs (GfsAdaptiveRK * RKdata, RKFunc derivatives)
{
  guint i;
  gdouble eps = 1.e-4, dtnext, errmax;
  gdouble dt, dttmp, tnew, *yerr, *ytmp;
  gdouble *ak1;
     
  ak1 = g_malloc (RKdata->n*sizeof (gdouble));
  yerr = g_malloc (RKdata->n*sizeof (gdouble));     
  ytmp = g_malloc (RKdata->n*sizeof (gdouble));     
    
  dt = MIN (RKdata->dtgoal, RKdata->dtmax); //step size to the trial value
    
  while (dt > 0.) {
    derivatives (RKdata->t, RKdata->y, ak1, RKdata->data);
    for (;;) {
      errmax = 0.;
      rkck (&dt, &errmax, yerr, ytmp, ak1, RKdata, derivatives);
      for (i = 0; i < RKdata->n; i++) 
	errmax = MAX (errmax, fabs (yerr[i]/RKdata->yc[i]));
      errmax /= eps;
      if (errmax <= 1.)
	break;
      dttmp = SAFETY*dt*pow(errmax,PSHRNK);
      dt = MAX (dttmp, 0.1*dt);
      tnew = RKdata->t + dt;
      if (tnew == RKdata->t)
	g_error ("Time step in RK equal 0");
    }

    if (errmax > ERRCON) 
      dtnext = SAFETY*dt*pow (errmax, PGROW);
    else 
      dtnext = 5.0*dt;

    RKdata->t += dt;
    for (i = 0; i < RKdata->n; i++) 
      RKdata->y[i] = ytmp[i];
    RKdata->dtgoal -= dt;
   
    dt = MIN (MIN (dtnext, RKdata->dtgoal), RKdata->dtmax);
  }

  g_free (ak1); g_free (yerr); g_free (ytmp);
}


/* GfsParticulate: Header */

typedef struct _GfsParticulate GfsParticulate;

struct _GfsParticulate {
  GfsParticle parent;
  FttVector vel;
  gdouble mass, volume;
  FttVector force;
  GtsSListContainer * forces;
};

#define GFS_PARTICULATE(obj)            GTS_OBJECT_CAST (obj,		\
							 GfsParticulate, gfs_particulate_class ())
#define GFS_IS_PARTICULATE(obj)         (gts_object_is_from_class (obj, gfs_particulate_class ()))

static GfsEventClass * gfs_particulate_class  (void);

/* GfsParticleList: Header */

typedef struct _GfsParticleList GfsParticleList;

struct _GfsParticleList {
  GfsEventList parent;
  gint idlast;
  GtsSListContainer * forces;
};

#define GFS_PARTICLE_LIST(obj)            GTS_OBJECT_CAST (obj,		\
							   GfsParticleList, \
							   gfs_particle_list_class ())

#define GFS_IS_PARTICLE_LIST(obj)         (gts_object_is_from_class (obj, \
								     gfs_particle_list_class ()))

static GfsEventClass * gfs_particle_list_class  (void);

/* GfsParticleForce: header */

typedef struct _GfsParticleForce GfsParticleForce;

struct _GfsParticleForce{
  GtsSListContainee parent;
  FttVector (* force) (GfsParticle *p, GfsParticleForce *force);
};

#define GFS_PARTICLE_FORCE(obj)            GTS_OBJECT_CAST (obj,		\
							GfsParticleForce, \
							gfs_particle_force_class ())
#define GFS_IS_PARTICLE_FORCE(obj)         (gts_object_is_from_class (obj, \
								      gfs_particle_force_class ()))

static GtsSListContaineeClass * gfs_particle_force_class  (void);

/* GfsForceCoeff: header */

typedef struct _GfsForceCoeff GfsForceCoeff;

struct _GfsForceCoeff{
  GfsParticleForce parent;
  GfsFunction * coefficient;
  GfsVariable *re_p, *u_rel, *v_rel, *w_rel, *pdia;
  GfsParticulate *p;
};

#define FORCE_COEFF(obj)            GTS_OBJECT_CAST (obj,		\
						    GfsForceCoeff,		\
						    gfs_force_coeff_class ())
#define GFS_IS_FORCE_COEFF(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_coeff_class ()))
static GtsSListContaineeClass * gfs_force_coeff_class  (void);

/* GfsForceLift: header */

#define GFS_IS_FORCE_LIFT(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_lift_class ()))
static GtsSListContaineeClass * gfs_force_lift_class  (void);

/* GfsForceDrag: header */

#define GFS_IS_FORCE_DRAG(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_drag_class ()))
static GtsSListContaineeClass * gfs_force_drag_class  (void);

/* GfsForceBuoy: header */

#define GFS_IS_FORCE_BUOY(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_buoy_class ()))
static GtsSListContaineeClass * gfs_force_buoy_class  (void);

/* Forces on the Particle */

static FttVector subs_fttvectors (FttVector *a, FttVector *b)
{
  FttVector result;
  FttComponent c;
  for(c = 0; c< FTT_DIMENSION; c++)    
    (&result.x)[c]  = (&a->x)[c] - (&b->x)[c];  
  return result;
}

/* Same as in source.c used here to obtained viscosity */
static GfsSourceDiffusion * source_diffusion_viscosity (GfsVariable * v)
{
  if (v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
    
    while (i) {
      GtsObject * o = i->data;
      
      if (GFS_IS_SOURCE_DIFFUSION (o))
	return GFS_SOURCE_DIFFUSION (o);
      i = i->next;
    }
  }
  return NULL;
}

/* Similar to gfs_vorticity which returns norm of the vorticity */
static void vorticity_vector (FttCell *cell, GfsVariable **v, 
			      FttVector *vort) 
{
  gdouble size;

  if (cell == NULL) return;
  if (v == NULL) return;

  size = ftt_cell_size (cell);
#if FTT_2D
  vort->x = 0.;
  vort->y = 0.;
  vort->z = (gfs_center_gradient (cell, FTT_X, v[1]->i) -
	     gfs_center_gradient (cell, FTT_Y, v[0]->i))/size;
#else  /* FTT_3D */
  vort->x = (gfs_center_gradient (cell, FTT_Y, v[2]->i) -
	     gfs_center_gradient (cell, FTT_Z, v[1]->i))/size;
  vort->y = (gfs_center_gradient (cell, FTT_Z, v[0]->i) -
	     gfs_center_gradient (cell, FTT_X, v[2]->i))/size;
  vort->z = (gfs_center_gradient (cell, FTT_X, v[1]->i) -
	     gfs_center_gradient (cell, FTT_Y, v[0]->i))/size;
#endif
}

/* GfsForceCoeff: object */

static void gfs_force_coeff_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_force_coeff_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_force_coeff_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '\n') {
    GfsForceCoeff * force = FORCE_COEFF (*o);
    force->coefficient = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (force->coefficient, gfs_object_simulation (*o), fp);
    GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
    
    /* fixme: "Rep", "Urelp" etc... should be derived variables not
       straight variables (i.e. there is no need to allocate memory
       for these as they are only used temporarily to compute the
       coefficient) */
    force->re_p = gfs_domain_get_or_add_variable (domain, "Rep", 
						  "Particle Reynolds number");  
    force->u_rel = gfs_domain_get_or_add_variable (domain, "Urelp", 
						   "Particle x - relative velocity");
    force->v_rel = gfs_domain_get_or_add_variable (domain, "Vrelp", 
						   "Particle y - relative velocity");
#if !FTT_2D
    force->w_rel = gfs_domain_get_or_add_variable (domain, "Wrelp", 
						   "Particle z - relative velocity");
#endif
    force->pdia = gfs_domain_get_or_add_variable (domain, "Pdia", 
						  "Particle radii");
  }
}

static void gfs_force_coeff_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_force_coeff_class ())->parent_class->write) (o, fp);
  GfsForceCoeff * force = FORCE_COEFF (o);
  if (force->coefficient)
    gfs_function_write (force->coefficient, fp);
}

static void gfs_force_coeff_destroy (GtsObject * o)
{
  if (FORCE_COEFF (o)->coefficient)
    gts_object_destroy (GTS_OBJECT (FORCE_COEFF (o)->coefficient));

  (* GTS_OBJECT_CLASS (gfs_force_coeff_class ())->parent_class->destroy) (o);
}

static void gfs_force_coeff_class_init (GtsObjectClass * klass)
{
  klass->read = gfs_force_coeff_read;
  klass->write = gfs_force_coeff_write;
  klass->destroy = gfs_force_coeff_destroy;
}
 
GtsSListContaineeClass * gfs_force_coeff_class (void)
{
  static GtsSListContaineeClass * klass = NULL;
  
  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_coeff_info = {
      "GfsForceCoeff",
      sizeof (GfsForceCoeff),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) gfs_force_coeff_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_force_class ()),
				  &gfs_force_coeff_info);
  }
  return klass;
}

/* GfsForceLift: object */

static FttVector compute_lift_force (GfsParticle * p, GfsParticleForce * liftforce)
{
  GfsParticulate * particulate = GFS_PARTICULATE (p);
  GfsForceCoeff * coeff = FORCE_COEFF (liftforce);

  GfsSimulation * sim = gfs_object_simulation (particulate);
  GfsDomain * domain = GFS_DOMAIN (sim);
  
  FttVector force;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
  if (cell == NULL) return force;

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha, cell) : 1.;
  GfsVariable ** u = gfs_domain_velocity (domain);
 
  gdouble viscosity = 0.;
  GfsSourceDiffusion * d = source_diffusion_viscosity (u[0]); 
  if (d) viscosity = gfs_diffusion_cell (d->D, cell);
  
  FttVector fluid_vel;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&fluid_vel.x)[c] = gfs_interpolate (cell, p->pos, u[c]);

  FttVector relative_vel = subs_fttvectors (&fluid_vel, &particulate->vel);
  FttVector vorticity;
  vorticity_vector (cell, u, &vorticity);

  gdouble cl = 0.5;
  if (coeff->coefficient) {
    gdouble norm_relative_vel = sqrt (relative_vel.x*relative_vel.x + 
				      relative_vel.y*relative_vel.y +
				      relative_vel.z*relative_vel.z);
    gdouble dia =  2.*pow(3.0*(particulate->volume)/4.0/M_PI, 1./3.);
    if (viscosity == 0) {
      g_warning ("Viscosity is 0. cannot compute lift force on particulate\n");
      g_assert_not_reached ();
    }
    gdouble Re = norm_relative_vel*dia*fluid_rho/viscosity;

    GFS_VALUE (cell, coeff->re_p) = Re;
    GFS_VALUE (cell, coeff->pdia) = dia;
    GFS_VALUE (cell, coeff->u_rel) = relative_vel.x;
    GFS_VALUE (cell, coeff->v_rel) = relative_vel.y;
#if !FTT_2D
    GFS_VALUE (cell, coeff->w_rel) = relative_vel.z;
#endif
    cl = gfs_function_value (coeff->coefficient, cell); 
  }
 
#if FTT_2D
  force.x = fluid_rho*cl*relative_vel.y*vorticity.z;
  force.y = -fluid_rho*cl*relative_vel.x*vorticity.z;
#else
  force.x = fluid_rho*cl*(relative_vel.y*vorticity.z
			  -relative_vel.z*vorticity.y);
  force.y = fluid_rho*cl*(relative_vel.z*vorticity.x
			  -relative_vel.x*vorticity.z);
  force.z = fluid_rho*cl*(relative_vel.x*vorticity.y
			  -relative_vel.y*vorticity.x);
#endif

  return force; 
}

static void gfs_force_lift_init (GfsParticleForce * force)
{
  force->force = compute_lift_force;
}

GtsSListContaineeClass * gfs_force_lift_class (void)
{
  static GtsSListContaineeClass * klass = NULL;
  
  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_lift_info = {
      "GfsForceLift",
      sizeof (GfsForceCoeff),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_force_lift_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_force_coeff_class ()),
				  &gfs_force_lift_info);
  }
  return klass;
}

/* GfsForceDrag: object */

static FttVector compute_drag_force (GfsParticle * p, GfsParticleForce * dragforce)
{
  GfsParticulate * particulate = GFS_PARTICULATE (p);
  GfsForceCoeff * coeff = FORCE_COEFF (dragforce);
  GfsSimulation * sim = gfs_object_simulation (particulate);
  GfsDomain * domain = GFS_DOMAIN (sim);

  FttVector force;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
  if (cell == NULL) return force;

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha,cell) : 1.;
  GfsVariable ** u = gfs_domain_velocity (domain);  

  gdouble viscosity = 0.;  

  GfsSourceDiffusion * d = source_diffusion_viscosity (u[0]); 
  if (d) viscosity = gfs_diffusion_cell (d->D, cell);
  
  FttVector fluid_vel;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&fluid_vel.x)[c] = gfs_interpolate (cell, p->pos, u[c]);

  FttVector relative_vel = subs_fttvectors (&fluid_vel, &particulate->vel);

  gdouble dia = 2.*pow(3.0*(particulate->volume)/4.0/M_PI, 1./3.);
#if !FTT_2D
  gdouble norm_relative_vel = sqrt (relative_vel.x*relative_vel.x + 
				    relative_vel.y*relative_vel.y +
				    relative_vel.z*relative_vel.z);
#else
  gdouble norm_relative_vel = sqrt (relative_vel.x*relative_vel.x + 
				    relative_vel.y*relative_vel.y);
#endif

  gdouble cd = 0.;
  gdouble Re;
  if (viscosity == 0)    
    return force;
  else
    Re = norm_relative_vel*dia*fluid_rho/viscosity;
  
  if (coeff->coefficient) {
    GFS_VALUE (cell, coeff->re_p) = Re;
    GFS_VALUE (cell, coeff->u_rel) = relative_vel.x;
    GFS_VALUE (cell, coeff->v_rel) = relative_vel.y;
#if !FTT_2D
    GFS_VALUE (cell, coeff->w_rel) = relative_vel.z;
#endif
    GFS_VALUE (cell, coeff->pdia) = dia;
    cd = gfs_function_value (coeff->coefficient, cell); 
  }
  else {
    if (Re < 1e-8)
      return force;
    else if (Re < 50.0)
      cd = 16.*(1. + 0.15*pow(Re,0.5))/Re;
    else
      cd = 48.*(1. - 2.21/pow(Re,0.5))/Re;
  }
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] += 3./(4.*dia)*cd*norm_relative_vel*(&relative_vel.x)[c]*fluid_rho;
  
  return force;
}

static void gfs_force_drag_init (GfsParticleForce * force)
{
  force->force = compute_drag_force;
}

static GtsSListContaineeClass * gfs_force_drag_class (void)
{
  static GtsSListContaineeClass * klass = NULL;
  
  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_drag_info = {
      "GfsForceDrag",
      sizeof (GfsForceCoeff),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_force_drag_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_force_coeff_class ()),
				  &gfs_force_drag_info);
  }
  return klass;
}

/* GfsForceBuoy: object */

static FttVector compute_buoyancy_force (GfsParticle * p, GfsParticleForce * buoyforce)
{
  GfsParticulate * particulate = GFS_PARTICULATE (p); 
  GfsSimulation * sim = gfs_object_simulation (particulate);
  GfsDomain * domain = GFS_DOMAIN (sim);

  FttVector force;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
  if (cell == NULL) return force;

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha, cell) : 1.;
  GfsVariable ** u = gfs_domain_velocity (domain);
 
  gdouble g[3];
  for (c = 0; c < FTT_DIMENSION; c++) {
    g[c] = 0.;
    if (u[c]->sources) {
      GSList * i = GTS_SLIST_CONTAINER (u[c]->sources)->items;
      
      while (i) {
	if (GFS_IS_SOURCE (i->data)) {
	  g[c] += gfs_function_value (GFS_SOURCE ((GfsSourceGeneric *) i->data)->intensity, 
				      cell);
	}
	i = i->next;
      }
    }
  }

  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] += (particulate->mass/particulate->volume-fluid_rho)*g[c];

  return force;
}

static void gfs_force_buoy_init (GfsParticleForce * force)
{
  force->force = compute_buoyancy_force;
}

static GtsSListContaineeClass * gfs_force_buoy_class (void)
{
  static GtsSListContaineeClass * klass = NULL;
  
  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_buoy_info = {
      "GfsForceBuoy",
      sizeof (GfsParticleForce),
      sizeof (GtsSListContaineeClass), 
      (GtsObjectClassInitFunc) NULL, 
      (GtsObjectInitFunc) gfs_force_buoy_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_force_class ()),
				  &gfs_force_buoy_info);
  }
  return klass;
}

/* GfsParticleForce: object */

static void gfs_particle_force_read (GtsObject ** o, GtsFile * fp)
{ 
  GtsObjectClass *klass;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsParticleClass)");
    return;
  }

  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_particle_force_class ())) {
    gts_file_error (fp, "`%s' is not a GfsParticleForce", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
}

static void gfs_particle_force_write (GtsObject * o, FILE * fp)
{
  fprintf (fp, "%s", o->klass->info.name);
}

static void gfs_particle_force_class_init (GtsObjectClass * klass)
{
  GTS_OBJECT_CLASS(klass)->read = gfs_particle_force_read;
  GTS_OBJECT_CLASS(klass)->write = gfs_particle_force_write;
}

GtsSListContaineeClass * gfs_particle_force_class (void)
{
  static GtsSListContaineeClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_particle_force_info = {
      "GfsParticleForce",
      sizeof (GfsParticleForce),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) gfs_particle_force_class_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_slist_containee_class ()),
				  &gfs_particle_force_info);
  }
  return klass;
}

/* GfsParticulate: Object */

static void compute_forces (GfsParticleForce * event, GfsParticulate * p)
{ 
  FttComponent c;
  FttVector new_force = (event->force) (GFS_PARTICLE (p), event);
  FttVector total_force;
     
  for ( c = 0 ; c < FTT_DIMENSION; c++)
    (&total_force.x)[c] = (&new_force.x)[c]*p->volume + (&p->force.x)[c];
    
  p->force = total_force;
}

static gboolean gfs_particulate_event (GfsEvent * event, 
				       GfsSimulation * sim)
{
  GfsParticle * p = GFS_PARTICLE (event);
  GfsParticulate * particulate = GFS_PARTICULATE (event);
  FttVector pos = p->pos;
  gfs_simulation_map (sim, &pos);
  
  FttComponent c;
  /* Velocity Verlet Algorithm */
  for (c = 0; c < FTT_DIMENSION; c++) {
    (&pos.x)[c] += (&particulate->force.x)[c]*sim->advection_params.dt*sim->advection_params.dt
      /particulate->mass/2.+ (&particulate->vel.x)[c]*sim->advection_params.dt;
    (&particulate->vel.x)[c] += (&particulate->force.x)[c]*sim->advection_params.dt
      /(2.*particulate->mass);
  }
  
  /* Compute forces */
  if (particulate->forces != NULL) {
    for (c = 0; c < FTT_DIMENSION; c++)
      (&particulate->force.x)[c] = 0.;      
    gts_container_foreach (GTS_CONTAINER (particulate->forces), 
			   (GtsFunc) compute_forces, particulate);
  }
  
  for (c = 0; c < FTT_DIMENSION; c++)
    (&particulate->vel.x)[c] += 
      (&particulate->force.x)[c]*sim->advection_params.dt/(2.*particulate->mass);
  
  gfs_simulation_map_inverse (sim, &pos);
  p->pos = pos;   
  return TRUE;
} 

static void gfs_particulate_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  GfsParticulate * p = GFS_PARTICULATE (*o);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (mass)");
    return;
  }
  p->mass = atof (fp->token->str);
  gts_file_next_token (fp);
  
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (volume)");
    return;
  }
  p->volume = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (v.x)");
    return;
  }
  p->vel.x = atof (fp->token->str);
  gts_file_next_token (fp);
  
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (v.y)");
    return;
  }
  p->vel.y = atof (fp->token->str);
  gts_file_next_token (fp);
  
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (v.z)");
    return;
  }
  p->vel.z = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type == GTS_INT || fp->type == GTS_FLOAT) {
    p->force.x = atof (fp->token->str);
    gts_file_next_token (fp);
  }
  
  if (fp->type == GTS_INT || fp->type == GTS_FLOAT) {
    p->force.y = atof (fp->token->str);
    gts_file_next_token (fp);
  }
  
  if (fp->type == GTS_INT || fp->type == GTS_FLOAT) {
    p->force.z = atof (fp->token->str);
    gts_file_next_token (fp);
  }
}

static void gfs_particulate_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class->write) (o, fp);
 
 GfsParticulate * p = GFS_PARTICULATE (o);
  fprintf (fp, " %g %g %g %g %g", p->mass, p->volume, p->vel.x, p->vel.y, p->vel.z);
  fprintf (fp, " %g %g %g", p->force.x, p->force.y, p->force.z);
}

static void gfs_particulate_class_init (GfsEventClass * klass)
{
  klass->event = gfs_particulate_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_particulate_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_particulate_write;
}

GfsEventClass * gfs_particulate_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_particulate_info = {
      "GfsParticulate",
      sizeof (GfsParticulate),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_particulate_class_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_class ()),
				  &gfs_particulate_info);
  }
  return klass;
}

/* GfsParticleList: Object */

static void assign_forces(GfsParticulate *particulate, GtsSListContainer *forces)
{
  particulate->forces = forces;
}

static void gfs_particle_list_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsParticleList * p = GFS_PARTICLE_LIST (*o);  
  GfsEventList * l = GFS_EVENT_LIST (p);  
  if (fp->type == '{') {
    fp->scope_max++;
    gts_file_next_token (fp);
    
    while (fp->type == '\n')
      gts_file_next_token (fp);
  
    GfsSimulation * sim = gfs_object_simulation (*o);
    GtsObjectClass * klass;
    while (fp->type != '}') {      

      if (fp->type != GTS_STRING) {
	gts_file_error (fp, "expecting a keyword (GfsParticleForce)");
	break;
      }
      klass = gfs_object_class_from_name (fp->token->str);
 
      if (klass == NULL) {
	gts_file_error (fp, "unknown class `%s'", fp->token->str);
	break;
      }
      if (!gts_object_class_is_from_class (klass, gfs_particle_force_class ())) {
	gts_file_error (fp, "'%s' is not a GfsParticleForce", fp->token->str);
	break;
      }
  
      GtsObject * object = gts_object_new (klass);
      gfs_object_simulation_set (object, sim);
  
      (* klass->read) (&object, fp);
    
      if (fp->type == GTS_ERROR) {
	gts_object_destroy (object);
	break;
      }
  
      while (fp->type == '\n') 
	gts_file_next_token (fp);
      
      gts_container_add (GTS_CONTAINER (p->forces), GTS_CONTAINEE (object));   
    }
    
    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    fp->scope_max--;
    gts_file_next_token (fp);
 
  }

  p->forces->items = g_slist_reverse (p->forces->items);
 
  gts_container_foreach (GTS_CONTAINER (l->list), (GtsFunc) assign_forces, p->forces);

  if(fp->type == GTS_INT){
    p->idlast = atoi (fp->token->str);
    gts_file_next_token (fp);
  }    
}

static void gfs_particle_list_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->write) (o, fp);

  GfsParticleList * p = GFS_PARTICLE_LIST (o);
  fputs (" {\n", fp);
  GSList * i = p->forces->items;
  while (i) {
    fputs ("    ", fp);
    (* GTS_OBJECT (i->data)->klass->write) (i->data, fp);
    fputc ('\n', fp);
    i = i->next; 
  }
  fputc ('}', fp);

  fprintf (fp, " %d", p->idlast);
}

static void gfs_particle_list_init (GtsObject *o){

  GfsParticleList * plist = GFS_PARTICLE_LIST(o);

  plist->forces = 
    GTS_SLIST_CONTAINER (gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ())));
 
}

static void gfs_particle_list_destroy (GtsObject * o)
{
  GfsParticleList * plist = GFS_PARTICLE_LIST(o);
 
  gts_container_foreach (GTS_CONTAINER (plist->forces), (GtsFunc) gts_object_destroy, NULL);
  gts_object_destroy (GTS_OBJECT (plist->forces));

  (* GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->destroy) (o);
}

static void gfs_particle_list_class_init (GtsObjectClass * klass)
{
  klass->read = gfs_particle_list_read;
  klass->write = gfs_particle_list_write;  
  klass->destroy = gfs_particle_list_destroy;  
}

GfsEventClass * gfs_particle_list_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_particle_list_info = {
      "GfsParticleList",
      sizeof (GfsParticleList),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_particle_list_class_init,
      (GtsObjectInitFunc) gfs_particle_list_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_list_class ()),
				  &gfs_particle_list_info);
  }
  return klass;
}

/* GfsDropletToParticle: header */

typedef struct _GfsDropletToParticle GfsDropletToParticle;

struct _GfsDropletToParticle{
  /*< private >*/
  GfsParticleList parent;
  GfsVariable * v;
  
  /*< public >*/
  GfsFunction * fc;
  GfsVariable * c;
  gint min;
  gdouble resetwith;
  gdouble density;
};

#define DROPLET_TO_PARTICLE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsDropletToParticle,\
					         gfs_droplet_to_particle_class ())

#define IS_DROPLET_TO_PARTICLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_droplet_to_particle_class ()))

static GfsEventClass * gfs_droplet_to_particle_class  (void);

typedef struct {
  FttVector pos, vel;
  gdouble volume;
} Droplets;

typedef struct {
  GfsVariable * tag, * c, *t;
  Droplets * drops;
  GfsVariable **u;
  guint * sizes;
  guint n, min;
  gdouble resetval;
  gdouble density;
  GfsFunction *fc;
} DropletsPar;

static int greater (const void * a, const void * b)
{
  return *((guint *)a) > *((guint *)b) ? -1 : 1;
}

static void reset_small_fraction (FttCell * cell, DropletsPar * p)
{
  gint i = GFS_VALUE (cell, p->tag);
  if (i > 0 && p->sizes[i - 1] < p->min)
    GFS_VALUE (cell, p->c) = p->resetval;
}

static void compute_droplet_properties (FttCell * cell, DropletsPar * p)
{
  gint i = GFS_VALUE (cell, p->tag);
  gdouble h = ftt_cell_size (cell), vol;
  FttVector pos; 
  ftt_cell_pos (cell, &pos);
  GfsVariable ** u = p->u;

  if (i > 0) {
    p->sizes[i - 1]++;
    vol = pow (h, FTT_DIMENSION);
    p->drops[i-1].volume += vol*GFS_VALUE (cell, p->c);
    FttComponent c;
    for(c = 0; c < FTT_DIMENSION; c++){
      (&(p->drops[i-1].pos.x))[c] +=  (&pos.x)[c];
      (&(p->drops[i-1].vel.x))[c] += GFS_VALUE (cell,u[c]);
    }
  }  
}

static void convert_droplets (GfsDomain * domain, 
			      DropletsPar * pars, GfsParticleList * plist)
{
  GfsSimulation * sim = gfs_object_simulation (plist); 
  guint i;
  
  GfsDropletToParticle * d = DROPLET_TO_PARTICLE (plist);
  GfsEventList * l = GFS_EVENT_LIST (plist); 

  pars->sizes = g_malloc0 (pars->n*sizeof (guint));  
  pars->drops = g_malloc0 (pars->n*sizeof (Droplets));

  FttComponent c;
  /* Initialize drops */
  for (i = 0; i < pars->n; i++){
    pars->drops[i].volume = 0.;
    pars->sizes[i] = 0;
    for(c = 0; c < FTT_DIMENSION; c++) {
      (&(pars->drops[i].pos.x))[c] = 0.;
      (&(pars->drops[i].vel.x))[c] = 0.;
    }
  }
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) compute_droplet_properties, pars);

#ifdef HAVE_MPI
  if (domain->pid >= 0) {
    guint * sizes = g_malloc0 (pars->n*sizeof (guint));
    MPI_Allreduce (pars->sizes, sizes, pars->n, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    g_free (pars->sizes);
    pars->sizes = sizes;
  }
#endif
  if (d->min >= 0)
    pars->min = d->min;
  else {
    guint * tmp = g_malloc (pars->n*sizeof (guint));
    memcpy (tmp, pars->sizes, pars->n*sizeof (guint));
    qsort (tmp, pars->n, sizeof (guint), greater);
    g_assert (-1 - d->min < pars->n);
    /* fixme: this won't work for parallel jobs */
    pars->min = tmp[-1 - d->min];
    g_free (tmp);
  }
  
  for (i = 0; i < pars->n; i++) {
    if (pars->sizes[i] < pars->min){
      for (c = 0; c < FTT_DIMENSION; c++) {
      	(&pars->drops[i].pos.x)[c] = (&pars->drops[i].pos.x)[c]/pars->sizes[i];
	(&pars->drops[i].vel.x)[c] = (&pars->drops[i].vel.x)[c]/pars->sizes[i];
      }
      FttCell * cell = gfs_domain_locate (domain, pars->drops[i].pos, -1, NULL);    
      if (cell) {
	/* Construct an Object */
	GtsObjectClass * klass = l->klass;
	if (klass == NULL) {
	  gfs_error (0, "Unknown particle class\n");
	  return;
	}
	GtsObject * object = gts_object_new (klass);
	gfs_object_simulation_set (object, sim);
	l->list->items = g_slist_reverse (l->list->items);	
	gts_container_add (GTS_CONTAINER (l->list), GTS_CONTAINEE (object));
	l->list->items = g_slist_reverse (l->list->items);
	GfsEvent * list = GFS_EVENT (l);	
	gfs_event_set (GFS_EVENT (object), 
		       list->start, list->end, list->step, list->istart, list->iend, list->istep);
	GfsParticulate * drop = GFS_PARTICULATE (object);
	GfsParticle * p = GFS_PARTICLE (drop);
	
	drop->vel = pars->drops[i].vel;
	p->pos = pars->drops[i].pos;
	drop->volume = pars->drops[i].volume;
	p->id = ++plist->idlast;
	drop->mass = sim->physical_params.alpha ? 1./
	  gfs_function_value (sim->physical_params.alpha, cell) : 1.;
	drop->mass *= drop->volume;
	for (c = 0; c < FTT_DIMENSION; c++)
	  (&drop->force.x)[c] = 0.;
      }       
    }   
  }  
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) reset_small_fraction, pars); 
  g_free (pars->drops);
  g_free (pars->sizes);
 
}

/* GfsDropletToParticle: object */

static void compute_v (FttCell * cell, GfsRemoveDroplets * d)
{
  GFS_VALUE (cell, d->v) = gfs_function_value (d->fc, cell);
}

static gboolean gfs_droplet_to_particle_event (GfsEvent * event, GfsSimulation * sim)
{ 
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class)->event)
      (event, sim)) {  
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsParticleList * plist = GFS_PARTICLE_LIST (event);
    GfsDropletToParticle *d = DROPLET_TO_PARTICLE (event);
    d->v = d->fc ? gfs_function_get_variable (d->fc) : d->c;
    DropletsPar p ;   
  
    p.resetval = d->resetwith;
    p.tag = gfs_temporary_variable (domain);
    p.u = gfs_domain_velocity (domain);
    p.density = d->density;
    p.t = d->c;
  
    if (d->v){
      p.c = d->v;
      p.n = gfs_domain_tag_droplets (domain, p.c, p.tag);
      if (p.n > 0 && -d->min < (gint) p.n){
	p.c = d->c;
	convert_droplets (domain, &p, plist);
      }
    }
    else {      
      d->v = gfs_temporary_variable (domain);      
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
				(FttCellTraverseFunc) compute_v, d);
      p.c = d->v;      
      p.n = gfs_domain_tag_droplets (domain, p.c, p.tag);
      if (p.n > 0 && -d->min < (gint) p.n){
	p.c = d->c;	
	convert_droplets (domain, &p, plist);	
      }              
      gts_object_destroy (GTS_OBJECT (d->v));      
    } 

    gts_object_destroy (GTS_OBJECT (p.tag));
    return TRUE;
  }
  return FALSE;
}

static void gfs_droplet_to_particle_read (GtsObject ** o, GtsFile * fp)
{  
  if (GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsDropletToParticle * r = DROPLET_TO_PARTICLE (*o);  
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (r));

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (variable)");
    return;
  }

  if ((r->c = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_INT, "min",  TRUE},
      {GTS_DOUBLE, "reset",    TRUE},
      {GTS_DOUBLE, "density",   TRUE},
      {GTS_NONE}
    };

    var[0].data = &r->min;
    var[1].data = &r->resetwith;
    var[2].data = &r->density;

    gts_file_assign_variables (fp, var);
  }
 
  if (fp->type != '\n') {
    r->fc = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (r->fc, gfs_object_simulation (r), fp);
  }
}

static void gfs_droplet_to_particle_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->write) (o, fp);

  GfsDropletToParticle * r = DROPLET_TO_PARTICLE(o);

  fprintf (fp, " %s { min = %d reset = %g density = %g } ",
	   r->c->name, r->min, r->resetwith, r->density);
  if (r->fc)
    gfs_function_write (r->fc, fp);
}

static void gfs_droplet_to_particle_destroy (GtsObject * o)
{
  GfsDropletToParticle * drops = DROPLET_TO_PARTICLE (o);
  if (drops->fc)
    gts_object_destroy (GTS_OBJECT (drops->fc));
  
  (* GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->destroy) (o);
}

static void gfs_droplet_to_particle_init (GfsDropletToParticle * r)
{
  r->resetwith = 0.;
  r->min = 20;
  r->density = 1.;
}

static void gfs_droplet_to_particle_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_droplet_to_particle_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_droplet_to_particle_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_droplet_to_particle_write;  
  GTS_OBJECT_CLASS (klass)->destroy = gfs_droplet_to_particle_destroy;  
}

GfsEventClass * gfs_droplet_to_particle_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_droplet_to_particle_info = {
      "GfsDropletToParticle",
      sizeof (GfsDropletToParticle),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_droplet_to_particle_class_init,
      (GtsObjectInitFunc) gfs_droplet_to_particle_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_list_class ()),
				  &gfs_droplet_to_particle_info);
  }
  return klass;
}

/* GfsBubble: Header */

typedef struct _GfsBubble GfsBubble;

struct _GfsBubble {
  /*< private >*/
  GfsParticulate parent;
  
  /*< public >*/
  gdouble velR, p0, R0;
};

#define GFS_BUBBLE(obj)            GTS_OBJECT_CAST (obj, GfsBubble, gfs_bubble_class ())
#define GFS_IS_BUBBLE(obj)         (gts_object_is_from_class (obj, gfs_bubble_class ()))

static GfsEventClass * gfs_bubble_class  (void);

/* GfsBubble: Object */
/* The radius of each bubble varies according to the Rayleigh-Plesset equation */

typedef struct {
  gdouble liqpres, liqdens;
  GfsBubble * bubble;
} RPData;

void static bubble_derivs (gdouble t, gdouble * y, gdouble * dydt, RPData * rp)
{
  /* dVbdt */
  gdouble Rb = pow (y[0]*3./(4.*M_PI), 1./3.);
  dydt[0] = 4.*M_PI*pow (Rb, 2)*y[1];

  /* interface acceleration */
  gdouble pbubble = rp->bubble->p0*pow (rp->bubble->R0/Rb, 3.*1.4);
  
  /* incompressible RP equation */
  dydt[1] = ((pbubble - rp->liqpres)/rp->liqdens - 3./2.*pow (y[1], 2))/Rb;
}

static gboolean gfs_bubble_event (GfsEvent * event, 
				  GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class)->event) 
      (event, sim)) {
    GfsParticle * p = GFS_PARTICLE (event);
    GfsParticulate * particulate = GFS_PARTICULATE (event);
    GfsBubble * bubble = GFS_BUBBLE (event);
    GfsDomain * domain = GFS_DOMAIN (sim);

    GfsVariable * liqpres = gfs_variable_from_name (domain->variables, "P");
  
    FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
    if (cell == NULL) 
      return TRUE;
    gdouble liq_rho = sim->physical_params.alpha ? 1./
      gfs_function_value (sim->physical_params.alpha, cell) : 1.;

    FttVector pos = p->pos;
    gfs_simulation_map (sim, &pos);

    gdouble point_pres = gfs_interpolate (cell, p->pos, liqpres);

    gdouble y[2];
    y[0] = particulate->volume;
    y[1] = bubble->velR;
    gdouble yc[2];
    yc[0] = particulate->volume;           /* characteristic distance */
    yc[1] = 0.01*sqrt(bubble->p0/liq_rho); /* characteristic velocity */

    RPData rp = { point_pres, liq_rho, bubble };
    GfsAdaptiveRK RKdata;
    RKdata.t = sim->time.t;
    RKdata.dtgoal = sim->advection_params.dt;
    RKdata.dtmax = yc[1]*0.001;
    RKdata.y = y;
    RKdata.yc = yc;
    RKdata.n = 2;
    RKdata.data = &rp;

    rkqs (&RKdata, (RKFunc) bubble_derivs);

    bubble->velR = y[1];
    particulate->volume = y[0];

    return TRUE;
  }
  return FALSE;
} 

static void gfs_bubble_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  GfsBubble * p = GFS_BUBBLE (*o);
  GfsParticulate * part = GFS_PARTICULATE (*o);

  p->R0 = pow (part->volume*3./(4.*M_PI), 1./3.);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (radial velocity)");
    return;
  }
  p->velR = atof (fp->token->str);
  gts_file_next_token (fp);
  
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (reference pressure)");
    return;
  }
  p->p0 = atof (fp->token->str);
  gts_file_next_token (fp);
}

static void gfs_bubble_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class->write) (o, fp); 
  GfsBubble * p = GFS_BUBBLE (o);
  fprintf (fp, " %g %g", p->velR, p->p0);
}

static void gfs_bubble_class_init (GfsEventClass * klass)
{
  klass->event = gfs_bubble_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_bubble_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_bubble_write;
}

GfsEventClass * gfs_bubble_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bubble_info = {
      "GfsBubble",
      sizeof (GfsBubble),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_bubble_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particulate_class ()),
				  &gfs_bubble_info);
  }
  return klass;
}

/* Initialize module */

const gchar gfs_module_name[] = "particulates";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
  gfs_particulate_class ();
    gfs_bubble_class ();
  gfs_particle_list_class ();
  gfs_force_lift_class ();
  gfs_force_drag_class ();
  gfs_force_buoy_class ();
  gfs_particle_force_class ();

  gfs_droplet_to_particle_class ();
  return NULL; 
}
