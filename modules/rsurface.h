/* padding on 32 bits systems (to match automatic 64 bits padding) */

#if defined (__LP64__) || defined (__64BIT__) || defined (_LP64) || (__WORDSIZE == 64)
  #define PADDING_32_BITS
#else
  #define PADDING_32_BITS int padding
#endif

typedef struct _RSurface RSurface;

typedef struct { /* needs to be identical to typdirinfo in RStarTree.h */
  double m01, m02, m03;
  double m11, m13;
  double m22, m23, m33;
  double m44, m55, m66, m77;
  double m67, m76;
  double H0, H1, H2, H3, H4;
  double H5, H6;
  float Hmin, Hmax;
  int n;
  PADDING_32_BITS;
} RSurfaceSum;

typedef struct {
  float l, h;
} RSurfaceInterval;

typedef RSurfaceInterval  RSurfaceRect[2];

typedef int (* RSurfaceCheck) (RSurfaceRect rect, void * data, int depth);

RSurface * r_surface_new      (const char * fname, int size, FILE * fp);
RSurface * r_surface_open     (const char * fname, const char * mode, int size);
void       r_surface_update   (RSurface * rt);
int        r_surface_close    (RSurface * rt);
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
				       RSurfaceRect rect,
				       RSurfaceSum * sum);
const char * r_surface_name  (RSurface * rt);
int          r_surface_depth (RSurface * rt);
void         r_surface_info  (RSurface * rt);

typedef struct {
  int nlevel;       /* number of levels */
  int * nentries;   /* number of entries for each level */
  double ** aspect; /* aspect ratio for each entry */
} RSurfaceStats;

RSurfaceStats * r_surface_stats_new  (RSurface * rt, 
				      int level);
void            r_surface_stats_free (RSurfaceStats * s);
