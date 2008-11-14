#include "variable.h"

/* GfsRiver: Header */

#define GFS_RIVER_NVAR 3

typedef struct {
  /*< private >*/
  FttDirection d;
  gdouble cfl;

  /*< public >*/
  GfsVariable * v[GFS_RIVER_NVAR];
  GfsVariable * fv[FTT_NEIGHBORS][GFS_RIVER_NVAR];
  GfsVariable * flux[GFS_RIVER_NVAR];
  gdouble g, dt;
  GfsCenterGradient gradient;
} GfsRiverParams;

typedef struct _GfsRiver GfsRiver;

struct _GfsRiver {
  GfsSimulation parent;
  
  GfsRiverParams p;
};

#define GFS_RIVER(obj)            GTS_OBJECT_CAST (obj,\
					           GfsRiver,\
					           gfs_river_class ())
#define GFS_IS_RIVER(obj)         (gts_object_is_from_class (obj,\
						   gfs_river_class ()))

GfsSimulationClass * gfs_river_class        (void);
