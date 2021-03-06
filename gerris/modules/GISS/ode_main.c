/* GfsSolidMovingODE: Header */

typedef struct {
  /*< private >*/
  GfsSolidMoving parent;
  /*< public >*/  
} GfsSolidMovingODE;

static GfsEventClass * gfs_solid_moving_ode_class (void);

#define GFS_SOLID_MOVING_ODE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSolidMovingODE,\
					         gfs_solid_moving_ode_class ())
#define GFS_IS_SOLID_MOVING_ODE(obj)         (gts_object_is_from_class (obj,\
						 gfs_solid_moving_ode_class ()))

/* GfsSolidMovingODE: Object */

static void solid_moving_ode_destroy (GtsObject * object)
{
  dBodyDestroy (BDInfo[GFS_SOLID_MOVING_ODE (object)->parent.bdnum].body);
  (* GTS_OBJECT_CLASS (gfs_solid_moving_ode_class ())->parent_class->destroy) (object);
}

static void solid_moving_ode_read (GtsObject ** o, GtsFile * fp)
{
  GfsSolidMovingODE * solid = GFS_SOLID_MOVING_ODE (*o);
  L = gfs_object_simulation (GFS_SURFACE (*o))->physical_params.L;
  Ln = pow(L,4);     //SPODE: for scaling   
  solid->parent.bdnum = BDNUM; //SPODE: store which body this solid belongs to
  BDInfo[BDNUM].body= dBodyCreate (world); //SPODE: store the information in the index

  (* GTS_OBJECT_CLASS (gfs_solid_moving_ode_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR) return;

  GfsVariable ** v = gfs_domain_velocity (GFS_DOMAIN (gfs_object_simulation (solid)));
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    {
      GfsSurfaceGenericBc * bc = GFS_SURFACE_GENERIC_BC (gts_object_new (GTS_OBJECT_CLASS (gfs_surface_bc_ode_class ())));
      bc->v = v[c];
      bc->v->surface_bc = bc;
      GFS_SURFACE_BC_ODE (bc)->c = c;
    } //SPODE: apply the surface boundary conditions to the solid boundary
#include "ode_solid_info.c" //SPODE: read the body information
  BDNUM++; //SPODE: prepared to read the next body
}

static gboolean solid_moving_ode_event (GfsEvent * event, GfsSimulation * sim)
{ 
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_solid_moving_ode_class ())->parent_class)->event) (event, sim)) {
    if ( (BDNUM != 2) && (GFS_SOLID_MOVING_ODE (event)->parent.bdnum!=BDNUM-1) ) return TRUE; 
    //SPODE: not to work until here is the last body, the solver should handle every the body in the same time.
    static gdouble bdy[7];
    static FttVector f, m;
    static gint BD_Flag=1;
    gint i,j;
    FttVector force[BDNUM+6];
    //SPODE: these variables stores the data in current time step, fx: the force on x direction; mx the momentum on x axis
    if(BD_Flag--==1) {
	gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, (FttCellTraverseFunc) getboundary_ODE, bdy);  
        //printf("dr=%g\n",2* L / pow(2,bdy[6])); 
    }
    /*SPODE: calculating the force and momentum*/
    gfs_domain_solid_force_ODE (GFS_DOMAIN (sim), NULL); 
    /*SPODE: calculating the force and momentum*/
//begin loop
    for (i = 1; i < BDNUM; i++) {
      if (BDInfo[i].cradius > 0.){
        f.x = BDInfo[i].pf.x + BDInfo[i].vf.x;
        f.y = BDInfo[i].pf.y + BDInfo[i].vf.y;
        f.z = BDInfo[i].pf.z + BDInfo[i].vf.z;
        m.x = BDInfo[i].pm.x + BDInfo[i].vm.x;
        m.y = BDInfo[i].pm.y + BDInfo[i].vm.y;
        m.z = BDInfo[i].pm.z + BDInfo[i].vm.z;
      } //SPODE: smoothing the force and momentum on each body
      if (AF_Flag == 1) { //Using Artificial Force
        for(j = 1; j < BDNUM; j++) {
          getforce_ODE (i, j, bdy, force);
          f.x += force[j].x;
          f.y += force[j].y;
          f.z += force[j].z;
        }
      }
      dBodyAddForce (BDInfo[i].body, f.x * Ln * BDInfo[i].f_tx + BDInfo[i].a_ex.x * BDInfo[i].mass.mass, 
				     f.y * Ln * BDInfo[i].f_ty + BDInfo[i].a_ex.y * BDInfo[i].mass.mass, 
				     f.z * Ln * BDInfo[i].f_tz + BDInfo[i].a_ex.z * BDInfo[i].mass.mass);
      dBodyAddTorque (BDInfo[i].body, m.x * Ln * BDInfo[i].f_rx + BDInfo[i].m_ex.x * BDInfo[i].mass._I(1,1), 
				     m.y * Ln * BDInfo[i].f_ry + BDInfo[i].m_ex.y * BDInfo[i].mass._I(2,2), 
				     m.z * Ln * BDInfo[i].f_rz + BDInfo[i].m_ex.z * BDInfo[i].mass._I(3,3)); //SP_ODE: update the force and momentum  
      if(sim->time.t>0.){
      if(BDInfo[i].OPInfo.buff_count==0){
      	BDInfo[i].OPInfo.pf.x = BDInfo[i].OPInfo.pf.y = BDInfo[i].OPInfo.pf.z = 0.;
      	BDInfo[i].OPInfo.vf.x = BDInfo[i].OPInfo.vf.y = BDInfo[i].OPInfo.vf.z = 0.;
      	BDInfo[i].OPInfo.f.x = BDInfo[i].OPInfo.f.y = BDInfo[i].OPInfo.f.z = 0.;
      	BDInfo[i].OPInfo.pm.x = BDInfo[i].OPInfo.pm.y = BDInfo[i].OPInfo.pm.z = 0.;
      	BDInfo[i].OPInfo.vm.x = BDInfo[i].OPInfo.vm.y = BDInfo[i].OPInfo.vm.z = 0.;
      	BDInfo[i].OPInfo.m.x = BDInfo[i].OPInfo.m.y = BDInfo[i].OPInfo.m.z = 0.;
      	BDInfo[i].OPInfo.c.x = BDInfo[i].OPInfo.c.y = BDInfo[i].OPInfo.c.z = 0.;
      	BDInfo[i].OPInfo.v.x = BDInfo[i].OPInfo.v.y = BDInfo[i].OPInfo.v.z = 0.;
      	BDInfo[i].OPInfo.e.x = BDInfo[i].OPInfo.e.y = BDInfo[i].OPInfo.e.z = 0.;
      	BDInfo[i].OPInfo.o.x = BDInfo[i].OPInfo.o.y = BDInfo[i].OPInfo.o.z = 0.;
        BDInfo[i].OPInfo.Mark.x = BDInfo[i].OPInfo.Mark.y = BDInfo[i].OPInfo.Mark.z = 0.;
      }
      if(BDInfo[i].OPInfo.buff_count++<BUFF){
        BDInfo[i].OPInfo.pf.x += BDInfo[i].pf.x/BUFF;
      	BDInfo[i].OPInfo.pf.y += BDInfo[i].pf.y/BUFF;
      	BDInfo[i].OPInfo.pf.z += BDInfo[i].pf.z/BUFF;
      	BDInfo[i].OPInfo.vf.x += BDInfo[i].vf.x/BUFF;
     	BDInfo[i].OPInfo.vf.y += BDInfo[i].vf.y/BUFF;
  	BDInfo[i].OPInfo.vf.z += BDInfo[i].vf.z/BUFF;
      	BDInfo[i].OPInfo.f.x += BDInfo[i].pf.x + BDInfo[i].vf.x;
      	BDInfo[i].OPInfo.f.y += BDInfo[i].pf.y + BDInfo[i].vf.y;
      	BDInfo[i].OPInfo.f.z += BDInfo[i].pf.z + BDInfo[i].vf.z;
      	BDInfo[i].OPInfo.pm.x += BDInfo[i].pm.x/BUFF;
      	BDInfo[i].OPInfo.pm.y += BDInfo[i].pm.y/BUFF;
      	BDInfo[i].OPInfo.pm.z += BDInfo[i].pm.z/BUFF;
      	BDInfo[i].OPInfo.vm.x += BDInfo[i].vm.x/BUFF;
      	BDInfo[i].OPInfo.vm.y += BDInfo[i].vm.y/BUFF;
      	BDInfo[i].OPInfo.vm.z += BDInfo[i].vm.z/BUFF;
      	BDInfo[i].OPInfo.m.x += BDInfo[i].pm.x + BDInfo[i].vm.x;
      	BDInfo[i].OPInfo.m.y += BDInfo[i].pm.y + BDInfo[i].vm.y;
      	BDInfo[i].OPInfo.m.z += BDInfo[i].pm.z + BDInfo[i].vm.z;
      	BDInfo[i].OPInfo.Mark.x += BDInfo[i].Mark.x/BUFF;
      	BDInfo[i].OPInfo.Mark.y += BDInfo[i].Mark.y/BUFF;
      	BDInfo[i].OPInfo.Mark.z += BDInfo[i].Mark.z/BUFF;
      	const dReal *tmp; 
      	tmp = dBodyGetPosition (BDInfo[i].body);  
      	BDInfo[i].OPInfo.c.x += tmp[0]/BUFF;
      	BDInfo[i].OPInfo.c.y += tmp[1]/BUFF;
      	BDInfo[i].OPInfo.c.z += tmp[2]/BUFF;
      	tmp = dBodyGetLinearVel (BDInfo[i].body);
      	BDInfo[i].OPInfo.v.x += tmp[0]/BUFF;
     	BDInfo[i].OPInfo.v.y += tmp[1]/BUFF;
      	BDInfo[i].OPInfo.v.z += tmp[2]/BUFF;
      	tmp = dBodyGetQuaternion (BDInfo[i].body);
	gdouble q0, q1, q2, q3; 
      	q0 = tmp[0]; q1 = tmp[1]; q2 = tmp[2]; q3 = tmp[3];
      	BDInfo[i].OPInfo.e.x += atan2(2.*(q0*q1+q2*q3), 1.-2.*(q1*q1+q2*q2))/PI*180./BUFF;
      	BDInfo[i].OPInfo.e.y += asin(2.*(q0*q2-q3*q1))/PI*180./BUFF;
      	BDInfo[i].OPInfo.e.z += atan2(2.*(q0*q3+q1*q2), 1.-2.*(q2*q2+q3*q3))/PI*180./BUFF;
      	tmp = dBodyGetAngularVel (BDInfo[i].body);
      	BDInfo[i].OPInfo.o.x += tmp[0]/BUFF;
      	BDInfo[i].OPInfo.o.y += tmp[1]/BUFF;
      	BDInfo[i].OPInfo.o.z += tmp[2]/BUFF;
      }
      if(BDInfo[i].OPInfo.buff_count==BUFF)
      	BDInfo[i].OPInfo.buff_count = 0;
      }
    }
//end loop
    if (Step_SP==0.) {//usual
	if(K>0) K--;
	if(K==0) dWorldStep (world, sim->advection_params.dt);
    }	
    if (Step_SP!=0.) {//Buff;
	if(K>0) K--;
	if(K==0) sim->advection_params.dt=Step_SP;
	dWorldStep (world, sim->advection_params.dt);
    }
    if (AF_Flag == 2) collision_ODE(bdy);//Using Direct Answer
    return TRUE;
   }
 return FALSE;
}

static void solid_moving_ode_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = solid_moving_ode_destroy;
  GTS_OBJECT_CLASS (klass)->read = solid_moving_ode_read;
  klass->event = solid_moving_ode_event;
}

static void solid_moving_ode_init (GfsSolidMovingODE * solid)
{ }

static GfsEventClass * gfs_solid_moving_ode_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo solid_moving_ode_info = {
      "GfsSolidMovingODE",
      sizeof (GfsSolidMovingODE),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) solid_moving_ode_class_init,
      (GtsObjectInitFunc) solid_moving_ode_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_solid_moving_class ()), 
				  &solid_moving_ode_info);
  }

  return klass;
}

