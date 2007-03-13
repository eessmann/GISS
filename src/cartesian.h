/* GfsCartesianGrid: Header */

#include <gts.h>

typedef struct _GfsCartesianGrid      GfsCartesianGrid;

struct _GfsCartesianGrid {
  /*< private >*/
  GtsObject parent;
  guint N;      // Number of dimension
  guint * n;    // Size of each dimension
  gdouble ** x; // Position of each point in the grid
  gdouble * v;  // Data

  /*< public >*/
  /* add extra data here (if public) */
};

#define GFS_CARTESIAN_GRID(obj)            GTS_OBJECT_CAST (obj,\
                                       GfsCartesianGrid,\
                                      gfs_cartesian_grid_class ())
#define GFS_IS_CARTESIAN_GRID(obj)         (gts_object_is_from_class (obj,\
                                       gfs_cartesian_grid_class ()))

GtsObjectClass * gfs_cartesian_grid_class  (void);
GfsCartesianGrid * gfs_cartesian_grid_new    (GtsObjectClass * klass);
gboolean gfs_cartesian_grid_interpolate (GfsCartesianGrid * g, gdouble * p, gdouble * val);
