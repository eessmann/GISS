#include "variable.h"

/* GfsRiver: Header */

#define GFS_RIVER_NVAR 3

typedef struct _GfsRiver GfsRiver;

struct _GfsRiver {
  /*< private >*/
  GfsSimulation parent;
  FttDirection d;
  gdouble cfl;

  /*< public >*/
  GfsVariable * v[GFS_RIVER_NVAR], * v1[GFS_RIVER_NVAR], * dv[FTT_DIMENSION][GFS_RIVER_NVAR];
  GfsVariable * flux[GFS_RIVER_NVAR];
  gdouble g, dt;
  GfsCenterGradient gradient;  
};

#define GFS_RIVER(obj)            GTS_OBJECT_CAST (obj,\
					           GfsRiver,\
					           gfs_river_class ())
#define GFS_IS_RIVER(obj)         (gts_object_is_from_class (obj,\
						   gfs_river_class ()))

GfsSimulationClass * gfs_river_class        (void);

/* GfsBcSubcritical: Header */

#define GFS_IS_BC_SUBCRITICAL(obj)         (gts_object_is_from_class (obj,\
						 gfs_bc_subcritical_class ()))

GfsBcClass * gfs_bc_subcritical_class  (void);
