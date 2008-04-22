typedef struct _RSurface RSurface;

RSurface * r_surface_new      (const char * fname, int size, FILE * fp);
RSurface * r_surface_open     (const char * fname, const char * mode, int size);
void       r_surface_close    (RSurface * rt);
int        r_surface_insert   (RSurface * rt, double p[3], int id);

typedef int (* RSurfaceQuery) (double p[3], void * user_data);
void       r_surface_query    (RSurface * rt, 
			       double a[2], double b[2], double c[2], double d[2],
			       RSurfaceQuery q, void * user_data);
void       r_surface_query_region (RSurface * rt, 
				   double min[2], double max[2],
				   RSurfaceQuery q, void * user_data);
const char * r_surface_name (RSurface * rt);
