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

#include "droplet2particle.h"

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
  g_return_if_fail (cell != NULL);

  gint i = GFS_VALUE (cell, p->tag);
  gdouble h = ftt_cell_size(cell),vol;
  FttVector pos; 
  ftt_cell_pos(cell, &pos);
  GfsVariable **u = p->u;

  if (i > 0){
    p->sizes[i - 1]++;
    vol = pow(h, FTT_DIMENSION);
    p->drops[i-1].volume += vol*GFS_VALUE(cell, p->c);
    FttComponent c;
    for(c = 0; c < FTT_DIMENSION; c++){
      (&(p->drops[i-1].pos.x))[c] +=  (&pos.x)[c];
      (&(p->drops[i-1].vel.x))[c] += GFS_VALUE(cell,u[c]);
    }
  }  
}

static void convert_droplets (GfsDomain * domain, 
			      DropletsPar *pars, GfsParticleList * plist)
{
  
  GfsSimulation *sim = gfs_object_simulation(plist); 
  
  g_return_if_fail (pars->c != NULL);
  
  guint i;
  
  GfsDropletToParticle *d = DROPLET_TO_PARTICLE(plist);
  GfsEventList *l = GFS_EVENT_LIST(plist); 

  pars->sizes = g_malloc0 (pars->n*sizeof (guint));
  
  pars->drops = g_malloc0 (pars->n*sizeof (Droplets));

  FttComponent c;
  /* Initialize drops */
  for(i = 0; i < pars->n; i++){
    pars->drops[i].volume = 0.;
    pars->sizes[i] = 0;
    for(c = 0; c < FTT_DIMENSION; c++){
      (&(pars->drops[i].pos.x))[c] = 0.;
      (&(pars->drops[i].vel.x))[c] = 0.;
    }
  }
  /*Compute Droplet Properties here*/
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
  
  for(i = 0; i < pars->n; i++){    
    printf("1)size: %d, pos.x: %g, min:%d, index: %d\n", pars->sizes[i],(&pars->drops[i].pos.x)[1]/pars->sizes[i],pars->min, i);
    if(pars->sizes[i] < pars->min){
      printf("2)size: %d, min:%d, index: %d\n", pars->sizes[i],pars->min, i);
      for(c = 0; c < FTT_DIMENSION; c++){
      	(&pars->drops[i].pos.x)[c] = (&pars->drops[i].pos.x)[c]/pars->sizes[i];
	(&pars->drops[i].vel.x)[c] = (&pars->drops[i].vel.x)[c]/pars->sizes[i];
      }
      FttCell * cell = gfs_domain_locate(domain, pars->drops[i].pos, -1, NULL);    
      if(cell){

	/* Construct an Object */
	GtsObjectClass * klass = l->klass;
	if (klass == NULL) {
	  gfs_error (0,"Unknown particle class\n");
	  return;
	}
	GtsObject * object = gts_object_new (klass);
	gfs_object_simulation_set (object, sim);
	l->list->items = g_slist_reverse (l->list->items);	
	gts_container_add (GTS_CONTAINER (l->list), GTS_CONTAINEE (object));
	l->list->items = g_slist_reverse (l->list->items);
	GfsEvent *list = GFS_EVENT(l);	
	gfs_event_set (GFS_EVENT(object), list->start, list->end, list->step, list->istart, list->iend, list->istep);
	GfsParticulate * drop = GFS_PARTICULATE(object);
	GfsParticle *p = GFS_PARTICLE(drop);
	
	drop->vel =  pars->drops[i].vel;
	p->pos =  pars->drops[i].pos;
	drop->volume = pars->drops[i].volume;
	p->id = ++plist->idlast;
	drop->mass = sim->physical_params.alpha ? 1./
	  gfs_function_value(sim->physical_params.alpha, cell) : 1.;
	drop->mass *= drop->volume;
	for(c = 0; c < FTT_DIMENSION; c++){
	  (&drop->force.x)[c] = 0.;
	}
      }       
    }   
  }  
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) reset_small_fraction, pars); 
  g_free(pars->drops);
  g_free(pars->sizes);
 
}
/* Droplet to Particle */
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
    p.u = gfs_domain_velocity(domain);
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

    g_free(p.tag);
    return TRUE;
  }
  return FALSE;
}

static void gfs_droplet_to_particle_read (GtsObject ** o, GtsFile * fp)
{
  
  if (GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsDropletToParticle * r = DROPLET_TO_PARTICLE(*o);
  
  GfsDomain * domain = GFS_DOMAIN(gfs_object_simulation(r));

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (variable)");
    return;
  }

  if ((r->c = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if(fp->type == '{') {
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

  /* call write method of parent */
  (* GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->write) (o, fp);

  GfsDropletToParticle * r = DROPLET_TO_PARTICLE(o);

  fprintf (fp, " %s ", r->c->name);
  fprintf (fp, " { ");
  fprintf (fp, " min = %d reset = %g density = %g } ", r->min, r->resetwith, r->density);
  
  if (r->fc){
    gfs_function_write (r->fc, fp);
  }
  
}

static void gfs_droplet_to_particle_destroy (GtsObject * o)
{
  GfsDropletToParticle * drops = DROPLET_TO_PARTICLE(o);
  if (drops->fc)
    gts_object_destroy (GTS_OBJECT (drops->fc));
    /* call destroy method of parent */
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
  /* define new methods and overload inherited methods here */
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

/* Initialize modules */
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
  gfs_droplet_to_particle_class ();
  return NULL; 
}
