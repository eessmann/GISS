typedef struct _RSurface RSurface;

typedef struct { /* needs to be identical to typdirinfo in RStarTree.h */
  double m01, m02, m03;
  double m11, m13;
  double m22, m23, m33;
  double H0, H1, H2, H3, H4;
  float Hmin, Hmax;
  int n;
} RSurfaceSum;

typedef struct {
  float l, h;
} RSurfaceInterval;

typedef RSurfaceInterval  RSurfaceRect[2];

typedef int (* RSurfaceCheck) (RSurfaceRect rect, void * data);

RSurface * r_surface_new      (const char * fname, int size, FILE * fp);
RSurface * r_surface_open     (const char * fname, const char * mode, int size);
void       r_surface_update   (RSurface * rt);
void       r_surface_close    (RSurface * rt);
int        r_surface_insert   (RSurface * rt, double p[3], int id);

typedef int (* RSurfaceQuery) (double p[3], void * user_data);
void       r_surface_query    (RSurface * rt, 
			       double a[2], double b[2], double c[2], double d[2],
			       RSurfaceQuery q, void * user_data);
void       r_surface_query_region (RSurface * rt, 
				   double min[2], double max[2],
				   RSurfaceQuery q, void * user_data);
void       r_surface_sum_init (RSurfaceSum * sum);
void       r_surface_query_region_sum (RSurface * rt,
				       RSurfaceCheck includes,
				       RSurfaceCheck intersects,
				       void * data,
				       RSurfaceSum * sum);
const char * r_surface_name (RSurface * rt);
