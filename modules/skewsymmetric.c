/* Gerris - The GNU Flow Solver                       (-*-C-*-)
 * Copyright (C) 2011 Daniel Fuster
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

#include "simulation.h"
#include "source.h"
#include "adaptive.h"
#include "output.h"
#include "solid.h"

/* GfsSkewSymmetric: Header */

typedef struct _GfsSkewSymmetric              GfsSkewSymmetric;

struct _GfsSkewSymmetric {
  /*< private >*/
  GfsSimulation parent;

  /*< public >*/
  gdouble beta;         /*parameter to define the position of the intermediate step*/
  GfsVariable * velfaces[FTT_NEIGHBORS],  *velold[FTT_NEIGHBORS];
};

#define GFS_SKEW_SYMMETRIC(obj)            GTS_OBJECT_CAST (obj,		\
							 GfsSkewSymmetric,	\
							 gfs_skew_symmetric_class ())
#define GFS_IS_SKEW_SYMMETRIC(obj)         (gts_object_is_from_class (obj,	\
								   gfs_skew_symmetric_class ()))

GfsSimulationClass * gfs_skew_symmetric_class  (void);

typedef struct {
  GfsVariable **velfaces , **velold , **u; 
  GfsVariable *p;
  gdouble * dt, beta; 
} FaceData;

typedef struct {
  GfsSourceDiffusion * d; 
  GfsFunction * alpha;
  FaceData * fd;
} DataDif;

/* GfsSkewSymmetric: Object */

static void gfs_skew_symmetric_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_skew_symmetric_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '{') return;

  fp->scope_max++;
  gts_file_next_token (fp);

  while (fp->type != GTS_ERROR && fp->type != '}') {
    if (fp->type == '\n') {
      gts_file_next_token (fp);
      continue;
    }
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a keyword");
      return;
    }
    else if (!strcmp (fp->token->str, "beta")) {
      gts_file_next_token (fp);
      if (fp->type != '=')
        gts_file_error (fp, "expecting `='");
      else {
        gts_file_next_token (fp);
        GFS_SKEW_SYMMETRIC(*o)->beta = atof (fp->token->str);
        gts_file_next_token (fp);
      }
    }
    else
      gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
  }

  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);

}

static void gfs_skew_symmetric_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_skew_symmetric_class ())->parent_class->write) (o, fp);

  fprintf (fp, " { beta = %g } \n", GFS_SKEW_SYMMETRIC(o)->beta);
}

static void gfs_skew_symmetric_run (GfsSimulation * sim);

static void gfs_skew_symmetric_class_init (GfsSimulationClass * klass) 
{
  GTS_OBJECT_CLASS (klass)->read  = gfs_skew_symmetric_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_skew_symmetric_write;
  klass->run  = gfs_skew_symmetric_run;
}

static void gfs_skew_symmetric_init (GfsSkewSymmetric * object)
{
  object->beta = 0.20;

  GfsDomain * domain = GFS_DOMAIN (object);
  GfsVariable  * velfaces[FTT_NEIGHBORS], * velold[FTT_NEIGHBORS];
  FttDirection d;
  gchar sufix[2];
  gchar *name, * descr;


  for (d = 0; d < FTT_NEIGHBORS; d++) {
    sprintf(sufix, "%d", (guint)d);
   
    name = g_strconcat ("Uface", sufix, NULL);
    descr = g_strconcat (sufix, "-component of face velocity", NULL);
    velfaces[d] = gfs_domain_add_variable (domain, name, descr);
    velfaces[d]->units = 1.;
    g_free(name); g_free(descr);

    name = g_strconcat ("Ufaceold", sufix, NULL);
    descr = g_strconcat (sufix, "-component of old face velocity", NULL);
    velold[d]   = gfs_domain_add_variable (domain, name, descr);
    velold[d]->units = 1.;
    g_free(name); g_free(descr);
  }

  gfs_variable_set_vector (velfaces, FTT_DIMENSION);
  gfs_variable_set_vector (velold, FTT_DIMENSION);
}

GfsSimulationClass * gfs_skew_symmetric_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_skew_symmetric_info = {
      "GfsSkewSymmetric",
      sizeof (GfsSkewSymmetric),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) gfs_skew_symmetric_class_init,
      (GtsObjectInitFunc) gfs_skew_symmetric_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()),
        &gfs_skew_symmetric_info);
  }

  return klass;
}

static void get_face_values (FttCell * cell, FaceData * fd)
{
  GfsStateVector * s = GFS_STATE (cell);

  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    FttComponent c = d/2;
    s->f[d].un = GFS_VALUE (cell, fd->u[c])/2.;
    if (ftt_cell_neighbor (cell, d))
      s->f[d].un += GFS_VALUE (ftt_cell_neighbor (cell, d), fd->u[c])/2.;
    else
      s->f[d].un  = 0;
  }
}

static void initialize_unold (FttCell * cell, FaceData * fd)
{
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) 
    GFS_VALUE (cell, fd->velold[d]) = GFS_VALUE (cell, fd->velfaces[d]);
}

static void get_velfaces (FttCell * cell, FaceData * fd)
{
  GfsStateVector * s = GFS_STATE (cell);
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) 
    GFS_VALUE (cell, fd->velfaces[d]) = s->f[d].un;
}

static void get_cell_values (FttCell * cell, FaceData * fd)
{
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    GFS_VALUE (cell, fd->u[c]) = (GFS_VALUE (cell, fd->velfaces[2*c]) +
				  GFS_VALUE (cell, fd->velfaces[2*c + 1]))/2.;
}

static void advance_face_values (FttCell * cell, FaceData * fd)
{
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    GFS_VALUE (cell, fd->velfaces[d]) = ((1.0 + fd->beta)*GFS_VALUE (cell, fd->velfaces[d]) -
					 fd->beta*GFS_VALUE (cell, fd->velold[d]));

}

static FttVector rpos[FTT_NEIGHBORS] = {
#if FTT_2D
  {1.,0.,0.}, {-1.,0.,0.}, {0.,1.,0.}, {0.,-1.,0.}
#else  /* FTT_3D */
  {1.,0.,0.}, {-1.,0.,0.}, {0.,1.,0.}, {0.,-1.,0.}, {0.,0.,1.}, {0.,0.,-1.}
#endif /* FTT_3D */
};

/* d: direction of the face required */
/* d2: cell direction with respect to cellref */
static gdouble interpolate_value_skew (FttCell * cellref,
				       FttDirection d, 
				       FttDirection * d2, 
				       FaceData * fd)
{
  guint i;
  FttCell * cell;
  if (d2)
    cell = ftt_cell_neighbor (cellref, *d2);
  else
    cell = cellref;

  if (!cell) {
    g_warning ("warning: no bc implemented\n");
    return 0;
  }

  guint lref = ftt_cell_level (cellref);
  guint l    = ftt_cell_level (cell);
  if (l < lref) {
    /*    return GFS_VALUE (cell,fd->velfaces[d]); 0th order interpolation */
    FttVector posref, posinterp; 
    ftt_cell_pos (cellref, &posref);
    gdouble size = ftt_cell_size (cellref);
    posinterp.x = posref.x + size*(rpos[*d2].x + rpos[d].x/2.);
    posinterp.y = posref.y + size*(rpos[*d2].y + rpos[d].y/2.);
    posinterp.z = posref.z + size*(rpos[*d2].z + rpos[d].z/2.);
    return gfs_interpolate (cell, posinterp, fd->u[d/2]);
  }
  if (!FTT_CELL_IS_LEAF (cell)) { 
    FttCellChildren child;
    guint n = ftt_cell_children_direction (cell, FTT_OPPOSITE_DIRECTION(*d2), &child);
    gdouble vel = 0.;
    for (i = 0; i < n; i++)
      if (child.c[i])
        vel += GFS_VALUE (child.c[i],fd->velfaces[d]);
      return vel/n;
  }
  else 
    return GFS_VALUE (cell,fd->velfaces[d]);
}

/* b Adaptative boolean */
static gdouble transverse_advection (FttCell * cell, 
				     FttComponent oc,
				     FttDirection d,
				     gdouble un,
				     FaceData * fd,
				     gboolean b)
{
  gdouble uauxbot, uauxtop;
  gdouble vn, vntop, vnbot, vndiag;
  FttDirection daux;

  if (!b) {
    vn      = interpolate_value_skew (cell,2*oc,NULL,fd);//interpolate_alpha_skew (cell,2*oc,NULL,fd);
    vntop   = interpolate_value_skew (cell,2*oc,&d  ,fd);//interpolate_alpha_skew (cell,2*oc,&d  ,fd);
    vndiag  = interpolate_value_skew (cell,2*oc+1,&d ,fd);//interpolate_alpha_skew (cell,2*oc+1,&d,fd); 
    vnbot   = interpolate_value_skew (cell,2*oc+1,NULL,fd);//interpolate_alpha_skew (cell,2*oc+1,NULL,fd);
    daux    = 2*oc;
    uauxtop = interpolate_value_skew (cell, d, &daux, fd);
    daux    = 2*oc+1;
    uauxbot = interpolate_value_skew (cell, d, &daux, fd);
  } else {
    daux    = FTT_OPPOSITE_DIRECTION(d);
    vn      = interpolate_value_skew (cell,2*oc,&daux, fd);//interpolate_alpha_skew (cell,2*oc,&daux,fd);
    vntop   = interpolate_value_skew (cell,2*oc,&daux, fd);//interpolate_alpha_skew (cell,2*oc,&daux,fd);
    vndiag  = interpolate_value_skew (cell,2*oc+1,NULL,fd);//interpolate_alpha_skew (cell,2*oc+1,NULL,fd); 
    vnbot   = interpolate_value_skew (cell,2*oc,&daux ,fd);//interpolate_alpha_skew (cell,2*oc,&daux,fd);
    daux    = 2*oc;
    uauxtop = interpolate_value_skew (cell, FTT_OPPOSITE_DIRECTION(d), &daux, fd);
    daux    = 2*oc+1;
    uauxbot = interpolate_value_skew (cell, FTT_OPPOSITE_DIRECTION(d), &daux, fd);
  }

  return (uauxtop*(vn + vntop) - uauxbot*(vnbot + vndiag))/4.;
}

static void advection_term (FttCellFace * face, FaceData * fd)
{
  gdouble un, unext, unprev;

  GfsStateVector * s = GFS_STATE (face->cell);
  FttComponent c = face->d/2;
  FttDirection d,daux = face->d;
  gboolean cond;

  un = GFS_VALUE (face->cell,fd->velfaces[face->d]);
  if ((face->d % 2 ) != 0 ) {
    cond = TRUE;
    d = FTT_OPPOSITE_DIRECTION (daux);
    unext     = interpolate_value_skew (face->cell, d,    NULL , fd);
    unprev    = interpolate_value_skew (face->cell, daux, &daux, fd); 
  }
  else { 
    cond = FALSE;
    d = daux;
    unext     = interpolate_value_skew (face->cell, d, &d, fd);
    unprev    = interpolate_value_skew (face->cell, FTT_OPPOSITE_DIRECTION(d), NULL, fd);
  }

  s->f[face->d].v = ((un + unext)*unext - (un + unprev)*unprev)/4.;
#if FTT_2D
  s->f[face->d].v += transverse_advection (face->cell, 
					   FTT_ORTHOGONAL_COMPONENT (c), d, un, fd, cond);
#else  /* FTT_3D */
  static FttComponent orthogonal[FTT_DIMENSION][2] = {
    {FTT_Y, FTT_Z}, {FTT_X, FTT_Z}, {FTT_X, FTT_Y}
  };
  s->f[face->d].v += transverse_advection (face->cell, orthogonal[c][0], d, un, fd, cond);
  s->f[face->d].v += transverse_advection (face->cell, orthogonal[c][1], d, un, fd, cond); 
#endif
}

static gdouble transverse_diffusion (FttCell * cell, 
				     FttComponent oc,
				     FttDirection d,
				     gdouble un,
				     FaceData * fd)
{
  FttDirection daux    = 2*oc;
  gdouble uauxtop = interpolate_value_skew (cell, d, &daux, fd);
  daux    = 2*oc+1;
  gdouble uauxbot = interpolate_value_skew (cell, d, &daux, fd);
  return (uauxtop - un) - (un - uauxbot);
}

static void diffusion_term (FttCellFace * face, DataDif * data)
{
  gdouble size = ftt_cell_size (face->cell); /* fixme: I need to account for the metric */
  gdouble un, unext, unprev;

  gdouble flux      = 0.;  
  gdouble viscosity  = data->alpha ? gfs_function_face_value (data->alpha, face) : 1.;
  viscosity *= gfs_diffusion_cell (data->d->D, face->cell);

  GfsStateVector * s = GFS_STATE (face->cell);

  FttDirection od = FTT_OPPOSITE_DIRECTION(face->d);

  un      = interpolate_value_skew (face->cell, face->d, NULL, data->fd);

  if ( (face->d % 2 ) != 0 ) {
    unext   = interpolate_value_skew (face->cell, od     , NULL, data->fd);
    unprev  = interpolate_value_skew (face->cell, face->d, &(face->d) , data->fd); 
  }
  else {
    unext   = interpolate_value_skew (face->cell, face->d, &(face->d), data->fd);
    unprev  = interpolate_value_skew (face->cell, od,      NULL,    data->fd);
  } 
  flux = (unext - un) - (un - unprev);


  FttComponent c = face->d/2;
#if FTT_2D
  FttComponent oc = FTT_ORTHOGONAL_COMPONENT (c);
  flux += transverse_diffusion(face->cell, oc, face->d, un, data->fd);
#else
  static FttComponent orthogonal[FTT_DIMENSION][2] = {
    {FTT_Y, FTT_Z}, {FTT_X, FTT_Z}, {FTT_X, FTT_Y}
  };
  flux += transverse_diffusion(face->cell, orthogonal[c][0], face->d, un, data->fd);
  flux += transverse_diffusion(face->cell, orthogonal[c][1], face->d, un, data->fd);
#endif 

  s->f[face->d].v -= viscosity*flux/size;
}

static void update_vel (FttCell * cell, FaceData * fd)
{
  GfsStateVector * s = GFS_STATE (cell);
  gdouble size = ftt_cell_size (cell);

  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    GFS_VALUE (cell, fd->velfaces[d]) = (GFS_VALUE (cell, fd->velfaces[d]) + 
					 fd->beta*GFS_VALUE (cell, fd->velold[d]))/(1.+fd->beta); 
    s->f[d].un = (2*fd->beta*GFS_VALUE (cell, fd->velfaces[d]) + 
		  (0.5-fd->beta)*GFS_VALUE (cell, fd->velold[d])- s->f[d].v*(*fd->dt)/size)/(0.5+fd->beta);
    GFS_VALUE (cell, fd->velold[d]) = GFS_VALUE (cell, fd->velfaces[d]);
  }
}

/* Same as in source.c used here to obtain viscosity (make it more general?) */
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

static void reset_face_velocity (const FttCellFace * face)
{                               
  if (GFS_FACE_FRACTION_RIGHT (face) == 0.)
    return;

  if (ftt_face_type (face) == FTT_FINE_COARSE) 
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION(face->d)].un = 0;
}

static void correct_face_velocity (const FttCellFace * face, GfsVariable * u0)
{                               
 gdouble u;

 if (GFS_FACE_FRACTION_RIGHT (face) == 0.)
   return;

  u = GFS_STATE (face->cell)->f[face->d].un;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION(face->d)].un = u;
    break;
  case FTT_FINE_COARSE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION(face->d)].un += 
      u*gfs_domain_face_fraction (u0->domain, face)/
      (gfs_domain_face_fraction_right (u0->domain, face)*FTT_CELLS_DIRECTION (face->d));
    break;
  default:
    g_assert_not_reached ();
  }
}

static void obtain_face_fluxes (const FttCell * cell)
{                               
  FttCellChildren child;
  GfsStateVector * s = GFS_STATE (cell);

  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    FttCell * neighbor = ftt_cell_neighbor (cell, d);
    if (neighbor) {
      if (!FTT_CELL_IS_LEAF (neighbor)) {
        guint i, n = ftt_cell_children_direction (neighbor, FTT_OPPOSITE_DIRECTION(d), &child);
        s->f[d].v = 0;
        for (i = 0; i < n; i++)
          if (child.c[i])  
            s->f[d].v += GFS_STATE (child.c[i])->f[FTT_OPPOSITE_DIRECTION(d)].v;
        s->f[d].v /= n;
      }
      else if ((d % 2) > 0 && ftt_cell_level(cell) == ftt_cell_level(neighbor))
        s->f[d].v = GFS_STATE (neighbor)->f[FTT_OPPOSITE_DIRECTION(d)].v;
    }
  }
}

static GfsVariable ** face_velocity (GfsDomain * domain, GfsSkewSymmetric * ss)
{
  FttComponent c;
  static gchar name[][7] = {"Uface0","Uface1","Uface2","Uface3","Uface4","Uface5"};

  g_return_val_if_fail (domain != NULL, NULL);
  
  for (c = 0; c < FTT_NEIGHBORS; c++) {
    GfsVariable * v = gfs_variable_from_name (domain->variables, name[c]);
    g_return_val_if_fail (v != NULL, NULL);
    ss->velfaces[c] = v;
  }
  return ss->velfaces;
}

static GfsVariable ** face_velocity_old (GfsDomain * domain, GfsSkewSymmetric * ss)
{
  FttComponent c;
  static gchar name[][10] = {"Ufaceold0","Ufaceold1","Ufaceold2","Ufaceold3","Ufaceold4","Ufaceold5"};

  g_return_val_if_fail (domain != NULL, NULL);
  
  for (c = 0; c < FTT_NEIGHBORS; c++) {
    GfsVariable * v = gfs_variable_from_name (domain->variables, name[c]);
    g_return_val_if_fail (v != NULL, NULL);
    ss->velold[c] = v;
  }
  return ss->velold;
}


static void gfs_skew_symmetric_momentum (GfsSimulation * sim, FaceData * fd, GfsVariable **gmac)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  FttComponent c;
  FttDirection d;
  /* it is used for implementation of viscosity (improve?) */
  GfsSourceDiffusion * dif = source_diffusion_viscosity (fd->u[0]);

  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) advance_face_values, fd);
  
  /* boundary conditions */
  for (d = 0; d <  FTT_NEIGHBORS; d++)
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, fd->velfaces[d]);

  gfs_domain_cell_traverse (domain, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) get_cell_values, fd);

  gfs_domain_face_traverse (domain, FTT_XYZ,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttFaceTraverseFunc) advection_term, fd); 

  if (dif) { 
    DataDif dd = { dif , sim->physical_params.alpha, fd };
    gfs_domain_face_traverse (domain, FTT_XYZ,
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttFaceTraverseFunc) diffusion_term, &dd); 
  }

  /* regularize flux at faces */
  for (c = 0; c <  FTT_DIMENSION; c++)
    gfs_domain_face_bc (domain, c, fd->u[c]);

  gfs_domain_cell_traverse (domain, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) obtain_face_fluxes, NULL);

  gfs_domain_cell_traverse (domain, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) update_vel, fd);

  gfs_velocity_face_sources (domain, fd->u, (*fd->dt), sim->physical_params.alpha, gmac);

  gfs_domain_face_traverse (domain, FTT_XYZ,
                            FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, 
                            (FttFaceTraverseFunc) reset_face_velocity, NULL);

  gfs_domain_face_traverse (domain, FTT_XYZ,
                            FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, 
                            (FttFaceTraverseFunc) correct_face_velocity, fd->u[0]);
}

static void gfs_skew_symmetric_run (GfsSimulation * sim)
{
  GfsVariable * p,  * res = NULL, * gmac[FTT_DIMENSION]; 
  GfsDomain * domain;
  GSList * i;

  domain = GFS_DOMAIN (sim);

  p = gfs_variable_from_name (domain->variables, "P");

  g_assert (p);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++) 
    gmac[c] = gfs_temporary_variable (domain);

  gfs_variable_set_vector (gmac, FTT_DIMENSION);

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  i = domain->variables;
  while (i) {
    if (GFS_IS_VARIABLE_RESIDUAL (i->data))
      res = i->data;
    i = i->next;
  }

  gfs_simulation_set_timestep (sim);

  GfsVariable ** u = gfs_domain_velocity (domain);
  GfsVariable ** velfaces = face_velocity (domain,  GFS_SKEW_SYMMETRIC(sim));
  GfsVariable ** velold   = face_velocity_old (domain,  GFS_SKEW_SYMMETRIC(sim));

  FaceData fd = { velfaces, velold, u, p, &sim->advection_params.dt, GFS_SKEW_SYMMETRIC(sim)->beta};

  if (sim->time.i == 0) {
    /* provisional solution to initialize the face velocities at t=0 */
    gfs_domain_cell_traverse (domain, 
        FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
        (FttCellTraverseFunc) get_face_values, &fd);
  
    gfs_mac_projection (domain,
			&sim->projection_params, 
			&sim->advection_params,
			p, sim->physical_params.alpha, gmac, NULL);
 
    gfs_domain_cell_traverse (domain, 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) get_velfaces, &fd);

    gfs_advance_tracers (domain, sim->advection_params.dt*fd.beta);
  }

  while (sim->time.t < sim->time.end && sim->time.i < sim->time.iend) {
    
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    gfs_skew_symmetric_momentum (sim, &fd, gmac);

    sim->advection_params.dt = sim->advection_params.dt*2.;
    gfs_mac_projection (domain,
			&sim->projection_params, 
			&sim->advection_params,
			p, sim->physical_params.alpha, gmac, NULL);
    sim->advection_params.dt = sim->advection_params.dt/2.;

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_half_do, sim); 
    gfs_domain_cell_traverse (domain, 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) get_velfaces, &fd);

    gfs_domain_cell_traverse (domain, 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) get_cell_values, &fd);

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gfs_simulation_adapt (sim);

    sim->time.t = sim->tnext;
    sim->time.i++;

    gfs_simulation_set_timestep (sim);
    gfs_advance_tracers (domain, sim->advection_params.dt);

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);

  for (c = 0; c < FTT_DIMENSION; c++) 
    gts_object_destroy (GTS_OBJECT (gmac[c]));

}

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "skewsymmetric";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_skew_symmetric_class ();
  return NULL;
} 
