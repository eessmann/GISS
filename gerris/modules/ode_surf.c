/* GfsSurfaceBcODE: Header */

typedef struct _GfsSurfaceBcODE         GfsSurfaceBcODE;

struct _GfsSurfaceBcODE {
  /*< private >*/
  GfsSurfaceGenericBc parent;

  /*< public >*/
  FttComponent c;  
};

#define GFS_SURFACE_BC_ODE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSurfaceBcODE,\
					         gfs_surface_bc_ode_class ())
#define GFS_IS_SURFACE_BC_ODE(obj)         (gts_object_is_from_class (obj,\
						 gfs_surface_bc_ode_class ()))

static GfsSurfaceGenericBcClass * gfs_surface_bc_ode_class (void);

/* GfsSurfaceBcODE: Object */

static void surface_bc_ode_read (GtsObject ** o, GtsFile * fp)
{
  /* defined only through SolidMovingODE */
}

static void surface_bc_ode_write (GtsObject * o, FILE * fp)
{
  /* defined only through SolidMovingODE */
}

static gdouble p2p_distance(FttVector * p1, FttVector * p2, gint axis) //SP_ODE: calculae the distance 
{
  FttVector r;  
  r.x = (p1->x - p2->x) * (axis==0? 0.:1.);
  r.y = (p1->y - p2->y) * (axis==1? 0.:1.);
  r.z = (p1->z - p2->z) * (axis==2? 0.:1.);
  return sqrt(pow(r.x, 2.)+pow(r.y, 2.)+pow(r.z, 2.));
}

static gdouble p2c_distance (FttVector * p1, gint num) //SP_ODE: calculate the distance between a point and the center of a body
{
  FttVector * p2;
  p2 = (FttVector *)dBodyGetPosition (BDInfo[num].body);
  switch(BDInfo[num].type) { 
    case 0: //SP_ODE: normal body
      return p2p_distance(p1, p2,-1)-BDInfo[num].cradius;
    case 1: //SP_ODE: OutsideWall
      return BDInfo[num].cradius-p2p_distance(p1,p2,BDInfo[num].axis);
    case 2: //SP_ODE: InsideWall
      return p2p_distance(p1,p2,BDInfo[num].axis) - BDInfo[num].cradius;
    default:
      return -1.;
  }
}

static void surface_bc_ode (FttCell * cell, GfsSurfaceGenericBc * b)
{
  gint n;
  dVector3 v;
  GfsSurfaceBcODE * bc = GFS_SURFACE_BC_ODE (b);
  FttVector p1;
  p1.x = GFS_STATE (cell)->solid->ca.x * L;
  p1.y = GFS_STATE (cell)->solid->ca.y * L;
  p1.z = GFS_STATE (cell)->solid->ca.z * L;

if (cell->bdnum == 0) {
    gdouble distanceO, distanceN;
    distanceO = p2c_distance(&p1, 1);
    gint bdO = 1;
    for (n=2; n < BDNUM; n++) {
      distanceN = p2c_distance(&p1,n);
      if (distanceN < distanceO) {
	bdO = n;
	distanceO = distanceN;
      }
    }
    cell->bdnum = bdO;
  }//SP_ODE: the cell will be assigned to the closest body

  dBodyGetPointVel (BDInfo[cell->bdnum].body, p1.x, p1.y, p1.z, v);
/*  if (cell->bdnum == 0) cell->bdnum = bc->bdnum;
  dBodyGetPointVel (BDInfo[bc->bdnum].body, p1.x, p1.y, p1.z, v);*/
  GFS_STATE (cell)->solid->fv = v[bc->c]/L;
  cell->flags |= GFS_FLAG_DIRICHLET;
}

static void surface_bc_ode_class_init (GfsSurfaceGenericBcClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = surface_bc_ode_read;
  GTS_OBJECT_CLASS (klass)->write = surface_bc_ode_write;
  klass->bc = surface_bc_ode;
}

static GfsSurfaceGenericBcClass * gfs_surface_bc_ode_class (void)
{
  static GfsSurfaceGenericBcClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_surface_bc_info = {
      "GfsSurfaceBcODE",
      sizeof (GfsSurfaceBcODE),
      sizeof (GfsSurfaceGenericBcClass),
      (GtsObjectClassInitFunc) surface_bc_ode_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_surface_generic_bc_class ()),
				  &gfs_surface_bc_info);
  }

  return klass;
}

