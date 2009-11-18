/* Gerris - The GNU Flow Solver			(-*-C-*-)
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
 * Author of the object: Gaurav Tomar
 */

#include "particulates.h"
#include "source.h"


static FttVector  compute_buoyancy_force(GfsParticle *p, GfsParticleForce *force);
static FttVector  compute_lift_force(GfsParticle *p, GfsParticleForce *force);
static FttVector  compute_drag_force(GfsParticle *p, GfsParticleForce *force);

/* !Forces on the Particle! */

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

/**Lift**/
static FttVector compute_lift_force(GfsParticle *p, GfsParticleForce *liftforce)
{

  GfsParticulate *particulate = GFS_PARTICULATE(p);
  ForceLift * lift = FORCE_LIFT(liftforce);

  GfsSimulation *sim = gfs_object_simulation(particulate);
  GfsDomain * domain = GFS_DOMAIN(sim);
  
  FttVector force;
  FttComponent c;
  for(c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate(domain, p->pos, -1, NULL);
  if(cell==NULL) return force;

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value(sim->physical_params.alpha,cell) : 1.;
  GfsVariable ** u = gfs_domain_velocity (domain);
 
  gdouble viscosity = 0.;
   
  GfsSourceDiffusion *d = source_diffusion_viscosity(u[0]); 
  if(d) viscosity = gfs_diffusion_cell(d->D, cell);
  
  FttVector fluid_vel;
  for(c = 0; c < FTT_DIMENSION; c++)
    (&fluid_vel.x)[c] = gfs_interpolate(cell, p->pos, u[c]);

  FttVector relative_vel = subs_fttvectors(&fluid_vel, &particulate->vel);

  FttVector vorticity;
  vorticity_vector (cell, u, &vorticity);

  gdouble cl = 0.5;

  if(lift->coefficient){
    gdouble norm_relative_vel = sqrt(relative_vel.x*relative_vel.x + 
				     relative_vel.y*relative_vel.y +
				     relative_vel.z*relative_vel.z);

    gdouble dia =  2.*pow(3.0*(particulate->volume)/4.0/M_PI, 1./3.);
    
    gdouble Re = norm_relative_vel*dia*fluid_rho/viscosity;

    GFS_VARIABLE(cell, lift->re_p->i) = Re;
    GFS_VARIABLE(cell, lift->pdia->i) = dia;
    GFS_VARIABLE(cell, lift->u_rel->i) = relative_vel.x;
    GFS_VARIABLE(cell, lift->v_rel->i) = relative_vel.y;
#if !FTT_2D
    GFS_VARIABLE(cell, lift->w_rel->i) = relative_vel.z;
#endif
    cl = gfs_function_value (lift->coefficient, cell); 
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

static void gfs_force_lift_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_force_lift_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_force_lift_class ())->parent_class->read) 
      (o, fp);

  if (fp->type == GTS_ERROR)
    return;

  ForceLift *force = FORCE_LIFT(*o);
  
  if (fp->type == '{') {

    fp->scope_max++;
    gts_file_next_token (fp);
    
    force->coefficient = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (force->coefficient, gfs_object_simulation (*o), fp);
    gts_file_next_token (fp);
    
    GfsDomain *domain = GFS_DOMAIN(gfs_object_simulation(*o));
    if(force->coefficient){
      
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
    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    fp->scope_max--;
  }
  gts_file_next_token (fp);
}

static void gfs_force_lift_write (GtsObject * o, FILE * fp)
{
  /* call write method of parent */
  (* GTS_OBJECT_CLASS (gfs_force_lift_class ())->parent_class->write) (o, fp);
  ForceLift *force = FORCE_LIFT(o);
  if(force->coefficient)
    gfs_function_write (force->coefficient, fp);
}

static void gfs_force_lift_init (GtsObject * o)
{
  GfsParticleForce *force = PARTICLE_FORCE(o);
  force->force = &compute_lift_force;
}

static void gfs_force_lift_class_init (GfsEventClass * klass)
{
  /* define new methods and overload inherited methods here */
  GTS_OBJECT_CLASS (klass)->read = gfs_force_lift_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_force_lift_write;
}
 
GtsSListContaineeClass * gfs_force_lift_class (void)
{
  static GtsSListContaineeClass * klass = NULL;
  
  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_lift_info = {
      "ForceLift",
      sizeof (ForceLift),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) gfs_force_lift_class_init,
      (GtsObjectInitFunc) gfs_force_lift_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_force_class ()),
				  &gfs_force_lift_info);
  }
  return klass;
}

/* Drag Force */
static FttVector compute_drag_force (GfsParticle *p, GfsParticleForce *dragforce)
{

  GfsParticulate *particulate = GFS_PARTICULATE(p);
  ForceDrag * drag = FORCE_DRAG(dragforce);

  GfsSimulation *sim = gfs_object_simulation(particulate);
  GfsDomain * domain = GFS_DOMAIN(sim);

  FttVector force;
  FttComponent c;
  for(c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate(domain, p->pos, -1, NULL);
  if(cell==NULL) return force;

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value(sim->physical_params.alpha,cell) : 1.;
  GfsVariable **u = gfs_domain_velocity (domain);  

  gdouble viscosity = 0.;  

  GfsSourceDiffusion *d = source_diffusion_viscosity(u[0]); 
  if(d) viscosity = gfs_diffusion_cell(d->D, cell);
  
  FttVector fluid_vel;
  for(c = 0; c < FTT_DIMENSION; c++)
    (&fluid_vel.x)[c] = gfs_interpolate(cell, p->pos, u[c]);

  FttVector relative_vel = subs_fttvectors(&fluid_vel, &particulate->vel);

  gdouble dia = 2.*pow(3.0*(particulate->volume)/4.0/M_PI, 1./3.);
#if !FTT_2D
  gdouble norm_relative_vel = sqrt(relative_vel.x*relative_vel.x + 
				   relative_vel.y*relative_vel.y +
				   relative_vel.z*relative_vel.z);
#else
  gdouble norm_relative_vel = sqrt(relative_vel.x*relative_vel.x + 
				   relative_vel.y*relative_vel.y);
#endif

  gdouble cd = 0.;
  gdouble Re;
  if(viscosity == 0)    
    return force;
  else
    Re = norm_relative_vel*dia*fluid_rho/viscosity;
  
  if(drag->coefficient){

    GFS_VARIABLE(cell, drag->re_p->i) = Re;

    GFS_VARIABLE(cell, drag->u_rel->i) = relative_vel.x;

    GFS_VARIABLE(cell, drag->v_rel->i) = relative_vel.y;

#if !FTT_2D
    GFS_VARIABLE(cell, drag->w_rel->i) = relative_vel.z;
#endif

    GFS_VARIABLE(cell, drag->pdia->i) = dia;

    cd = gfs_function_value (drag->coefficient, cell); 
  }
  else
    if(Re < 1e-8)
      return force;
    else if(Re < 50.0)
      cd = 16.*(1. + 0.15*pow(Re,0.5))/Re;
    else
      cd = 48.*(1. - 2.21/pow(Re,0.5))/Re;
  
  for(c = 0; c< FTT_DIMENSION; c++)
    (&force.x)[c] += 3./(4.*dia)*cd*norm_relative_vel*(&relative_vel.x)[c]*fluid_rho;
  
  return force;
}

static void gfs_force_drag_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_force_lift_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_force_lift_class ())->parent_class->read) 
      (o, fp);
  ForceDrag *force = FORCE_DRAG(*o);
  
  if (fp->type == '{') {

    fp->scope_max++;
    gts_file_next_token (fp);
    
    force->coefficient = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (force->coefficient, gfs_object_simulation (*o), fp);
    gts_file_next_token (fp);
    
    GfsDomain *domain = GFS_DOMAIN(gfs_object_simulation(*o));
    if(force->coefficient){
      
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
    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    fp->scope_max--;
  }
  gts_file_next_token (fp);
}


static void gfs_force_drag_init (GtsObject * o)
{
  GfsParticleForce *force = PARTICLE_FORCE(o);
  force->force = &compute_drag_force;
}


static void gfs_force_drag_write (GtsObject * o, FILE * fp)
{ 
  /* call write method of parent */
  (* GTS_OBJECT_CLASS (gfs_force_drag_class ())->parent_class->write) (o, fp);
  ForceDrag *force = FORCE_DRAG(o);
  if(force->coefficient)
    gfs_function_write (force->coefficient, fp);
}
static void gfs_force_drag_class_init (GfsEventClass * klass)
{
  /* define new methods and overload inherited methods here */
  GTS_OBJECT_CLASS (klass)->read = gfs_force_drag_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_force_drag_write;
}
 
GtsSListContaineeClass * gfs_force_drag_class (void)
{
  static GtsSListContaineeClass * klass = NULL;
  
  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_drag_info = {
      "ForceDrag",
      sizeof (ForceDrag),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) gfs_force_drag_class_init,
      (GtsObjectInitFunc) gfs_force_drag_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_force_class ()),
				  &gfs_force_drag_info);
  }
  return klass;
}

/* Buoyancy Force */
static FttVector compute_buoyancy_force(GfsParticle *p, GfsParticleForce *buoyforce)
{
 
  GfsParticulate *particulate = GFS_PARTICULATE(p);
 
  GfsSimulation *sim = gfs_object_simulation(particulate);
  GfsDomain * domain = GFS_DOMAIN(sim);

  FttVector force;
  FttComponent c;
  for(c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate(domain, p->pos, -1, NULL);
  if(cell==NULL) return force;

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value(sim->physical_params.alpha,cell) : 1.;
  GfsVariable ** u = gfs_domain_velocity (domain);
 
  gdouble g[3];

  for(c = 0; c < FTT_DIMENSION; c++){
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

  for(c = 0; c< FTT_DIMENSION; c++)
    (&force.x)[c] += (particulate->mass/particulate->volume-fluid_rho)*g[c];

  return force;
}

static void gfs_force_buoy_init (GtsObject * o)
{
  GfsParticleForce *force = PARTICLE_FORCE(o);
  force->force = &compute_buoyancy_force;
}

GtsSListContaineeClass * gfs_force_buoy_class (void)
{
  static GtsSListContaineeClass * klass = NULL;
  
  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_buoy_info = {
      "ForceBuoy",
      sizeof (ForceBuoy),
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

/* Particle Force */
static void gfs_particle_force_read (GtsObject ** o, GtsFile * fp)
{

  printf("Here\n");
  GtsObjectClass *klass;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsEventClass)");
    return;
  }

  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_particle_force_class ())) {
    gts_file_error (fp, "`%s' is not a GfsEvent", fp->token->str);
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
  /* define new methods and overload inherited methods here */
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

static void compute_forces(GfsParticleForce *event, GfsParticulate *p){
 
  FttVector force;
  force = (event->force)(GFS_PARTICLE(p), event);
 
  FttComponent c;
  for(c = 0; c < FTT_DIMENSION; c++)
    (&p->force.x)[c] += (&force.x)[c];

}
/* GfsParticulate: Object */
static gboolean gfs_particulate_event (GfsEvent * event, 
				       GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class)->event)
      (event, sim)) {

    GfsParticle * p = GFS_PARTICLE(event);
    GfsParticulate * particulate = GFS_PARTICULATE(event);
    FttVector pos = p->pos;
    gfs_simulation_map (sim, &pos);
 
    FttComponent c;
    /*Velocity Verlet Algorithm*/
    for (c = 0; c < FTT_DIMENSION; c++){
      (&pos.x)[c] += (&particulate->force.x)[c]*sim->advection_params.dt*sim->advection_params.dt
	/particulate->mass/2.;
      (&particulate->vel.x)[c] += (&particulate->force.x)[c]*sim->advection_params.dt
	/(2.*particulate->mass);
    }
      
    /*Compute forces*/
    if(particulate->forces != NULL){
      for(c = 0; c<FTT_DIMENSION;c++)
	(&particulate->force.x)[c] = 0.;
      
      gts_container_foreach (GTS_CONTAINER (particulate->forces), 
			     (GtsFunc) compute_forces, particulate);

    }
    
    for(c = 0; c<FTT_DIMENSION;c++)
      (&particulate->vel.x)[c] += (&particulate->force.x)[c]*sim->advection_params.dt/(2.*particulate->mass);
    
    gfs_simulation_map_inverse (sim, &pos);
    p->pos = pos;   
    
    return TRUE;
  }
  return FALSE;
} 

static void gfs_particulate_read (GtsObject ** o, GtsFile * fp)
{

  if (GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  GfsParticulate * p = GFS_PARTICULATE(*o);

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
  /* call write method of parent */
  (* GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class->write) (o, fp);
  GfsParticulate * p = GFS_PARTICULATE(o);
  fprintf(fp," %g %g %g %g %g", p->mass, p->volume, p->vel.x, p->vel.y, p->vel.z);
  fprintf(fp," %g %g %g", p->force.x, p->force.y, p->force.z);

}

static void gfs_particulate_class_init (GfsEventClass * klass)
{
  /* define new methods and overload inherited methods here */
  GFS_EVENT_CLASS (klass)->event = gfs_particulate_event;
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
  GfsParticleList * p = GFS_PARTICLE_LIST(*o);
  
  GfsEventList *l = GFS_EVENT_LIST(p);
  
  if (fp->type == '{') {

    fp->scope_max++;
    gts_file_next_token (fp);
    
    while (fp->type == '\n')
      gts_file_next_token (fp);
  
    GfsSimulation * sim = gfs_object_simulation (*o);
    GtsObjectClass * klass;
    while (fp->type != '}') {
      
      if (fp->type != GTS_STRING) {
	gts_file_error (fp, "expecting a keyword");
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
}



static void gfs_particle_list_write (GtsObject * o, FILE * fp)
{

  /* call write method of parent */
  (* GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->write) (o, fp);

  GfsParticleList * p = GFS_PARTICLE_LIST(o);
  fputs (" {\n", fp);
  GSList * i = p->forces->items;
  while (i) {
    fputs ("    ", fp);
    (* GTS_OBJECT (i->data)->klass->write) (i->data, fp);
    fputc ('\n', fp);
    i = i->next; 
  }
  fputc ('}', fp);
 
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

  /* call destroy method of parent */
  (* GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->destroy) (o);
}

static void gfs_particle_list_class_init (GfsEventClass * klass)
{
  /* define new methods and overload inherited methods here */
  GTS_OBJECT_CLASS (klass)->read = gfs_particle_list_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_particle_list_write;  
  GTS_OBJECT_CLASS (klass)->destroy = gfs_particle_list_destroy;  
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

/* Initialize modules */
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
  gfs_particulate_class ();
  gfs_particle_list_class ();
  gfs_force_lift_class ();
  gfs_force_drag_class ();
  gfs_force_buoy_class ();
  gfs_particle_force_class ();
  return NULL; 
}
