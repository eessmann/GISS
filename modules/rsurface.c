#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "RStarTree/RStarTree.h"
#include "kdt/kdt.h"
#include "rsurface.h"

/* RSurface */

struct _RSurface {
  RSTREE t;
  char * name;
  Kdt * kdt;
};

RSurface * r_surface_new (const char * fname, int size, FILE * fp)
{
  return NULL;
}

RSurface * r_surface_open (const char * fname, const char * mode, int size)
{
  RSurface * rt = malloc (sizeof (RSurface));
  rt->t = NULL;
  rt->kdt = NULL;
  if (!strcmp (mode, "w")) {
    RemoveRST (fname);
    if (!CreateRST (fname, size, FALSE)) {
      free (rt);
      return NULL;
    }
  }
  rt->kdt = kdt_new ();
  if (kdt_open (rt->kdt, fname)) {
    kdt_destroy (rt->kdt);
    rt->kdt = NULL;
    if (!OpenRST (&rt->t, fname, !strcmp (mode, "w") ? "rw" : "r")) {
      free (rt);
      return NULL;
    }
  }
  rt->name = malloc (sizeof (char)*(strlen (fname) + 1));
  strcpy (rt->name, fname);
  return rt;
}

int r_surface_close (RSurface * rt)
{
  int status = 1;
  if (rt->kdt)
    kdt_destroy (rt->kdt);
  else
    status = CloseRST (&rt->t);
  free (rt->name);
  free (rt);
  return status;
}

int r_surface_insert  (RSurface * rt, double p[3], int id)
{
  assert (!rt->kdt); /* not implemented */
  typrect rect;
  typinfo info;
  boolean inserted;
  rect[0].l = p[0]; rect[0].h = p[0];
  rect[1].l = p[1]; rect[1].h = p[1];
  info.height = p[2];
  return (InsertRecord (rt->t, rect, &info, &inserted) && inserted);
}

void r_surface_query (RSurface * rt, 
		      double a[2], double b[2], double c[2], double d[2],
		      RSurfaceQuery q, void * user_data)
{
}

static boolean Intersects (RSTREE R,
			   typrect RSTrect,
			   typrect queryrect,
			   typrect unused)
{
  int maxdim= NumbOfDim -1;
  boolean inter;
  int d;
  
  d= -1;
  do {
    d++;
    inter= RSTrect[d].l <= queryrect[d].h &&
           RSTrect[d].h >= queryrect[d].l;
  } while (inter && d != maxdim);
  return inter;
}

static void ManageQuery (RSTREE R,
			 typrect rectangle,
			 refinfo infoptr,
			 void ** data,
			 boolean *modify,
			 boolean *finish)
{
  RSurfaceQuery q = data[0];
  double p[3];
  p[0] = rectangle[0].l; p[1] = rectangle[1].l;
  p[2] = infoptr->height;
  (*q) (p, data[1]);
  *modify = FALSE;
  *finish = FALSE;
}

void r_surface_query_region (RSurface * rt, 
			     double min[2], double max[2],
			     RSurfaceQuery q, void * user_data)
{
  assert (!rt->kdt); /* not implemented */
  typrect rect, unused;
  rect[0].l = min[0]; rect[0].h = max[0];
  rect[1].l = min[1]; rect[1].h = max[1];
  void * data[2];
  data[0] = q;
  data[1] = user_data;
  RegionQuery (rt->t, rect, unused, Intersects, Intersects, (QueryManageProc) ManageQuery, data);
}

void r_surface_sum_init (RSurfaceSum * sum)
{
  memset (sum, 0, sizeof (RSurfaceSum));
  sum->Hmin = 1e30;
  sum->Hmax = - 1e30;
}

/**
 * Fills @sum using @rect to normalise the results i.e. the sums are
 * expressed in a coordinate system centered on 
 * (rect[0].l + rect[0].h, rect[1].l + rect[1].h)/2 
 * and scaled by 
 * MAX(rect[0].h - rect[0].l, rect[1].h - rect[1].l).
 */
void r_surface_query_region_sum (RSurface * rt,
				 RSurfaceCheck includes,
				 RSurfaceCheck intersects,
				 void * data,
				 RSurfaceRect rect,
				 RSurfaceSum * sum)
{
  if (rt->kdt)
    kdt_query_sum (rt->kdt, (KdtCheck) includes, (KdtCheck) intersects, data, 
		   (KdtInterval *) rect, (KdtSum *) sum);
  else
    RegionQueryInfo (rt->t, (Check) includes, (Check) intersects, data, 
		     (typinterval *) rect, (typdirinfo *) sum);
}

const char * r_surface_name (RSurface * rt)
{
  return rt->name;
}

void r_surface_update (RSurface * rt)
{
  assert (!rt->kdt); /* not implemented */
  Update (rt->t);
}

int r_surface_depth (RSurface * rt)
{
  assert (!rt->kdt); /* not implemented */
  int height = -1;
  GetHeight (rt->t, &height);
  return height;
}

void r_surface_info (RSurface * rt)
{
  assert (!rt->kdt); /* not implemented */
  char name[100];
  int numbofdim, sizedirentry, sizedataentry, sizeinfo;
  int maxdirfanout, maxdatafanout, pagesize, numbofdirpages;
  int numbofdatapages, pagesperlevel[100], numbofrecords, height;
  boolean unique;
  
  InquireRSTDesc(rt->t,
		 name,
		 &numbofdim,
		 &sizedirentry,
		 &sizedataentry,
		 &sizeinfo,
		 &maxdirfanout,
		 &maxdatafanout,
		 &pagesize,
		 &numbofdirpages,
		 &numbofdatapages,
		 pagesperlevel,
		 &numbofrecords,
		 &height,
		 &unique);
  fprintf (stderr, 
	   "Number of dimensions:\t\t\t\t%d\n"
	   "Size (bytes) of a directory entry:\t\t%d\n"
	   "Size (bytes) of a data entry:\t\t\t%d\n"
	   "Size (bytes) of an information part:\t\t%d\n"
	   "Maximum # of entries for a directory node:\t%d\n"
	   "Maximum # of entries for a data node:\t\t%d\n"
	   "Total # of directory pages:\t\t\t%d\n"
	   "Total # of data pages:\t\t\t\t%d\n"
	   "Total # of data records:\t\t\t%d\n"
	   "Are records unique?:\t\t\t\t%s\n"
	   "Height of the tree:\t\t\t\t%d\n",
	   numbofdim,
	   sizedirentry,
	   sizedataentry,
	   sizeinfo,
	   maxdirfanout,
	   maxdatafanout,
	   numbofdirpages,
	   numbofdatapages,
	   numbofrecords,
	   unique ? "yes" : "no",
	   height);
  int i;
  fprintf (stderr, "Number of pages per level:\n");
  for (i = 0; i < height; i++)
    fprintf (stderr, "\tlevel %d:\t%d\n", i + 1, pagesperlevel[i]);
}

static int includes (RSurfaceRect RSTrect, RSurfaceStats * s, int depth)
{
  return 0;
}

static int intersects (RSurfaceRect RSTrect, RSurfaceStats * s, int depth)
{
  double w = RSTrect[0].h - RSTrect[0].l, h = RSTrect[1].h - RSTrect[1].l;
  double ratio = 0.;
  if (w > 0. && h > 0.)
    ratio = w < h ? w/h : h/w;
  s->aspect[depth][s->nentries[depth]++] = ratio;
  return (depth < s->nlevel - 1);
}

static int compare_double (const void * p1, const void * p2)
{
  return (*(double * ) p1 < *(double * ) p2);
}

RSurfaceStats * r_surface_stats_new (RSurface * rt, int level)
{
  assert (!rt->kdt); /* not implemented */
  char name[100];
  int numbofdim, sizedirentry, sizedataentry, sizeinfo;
  int maxdirfanout, maxdatafanout, pagesize, numbofdirpages;
  int numbofdatapages, pagesperlevel[100], numbofrecords, height;
  boolean unique;
  
  InquireRSTDesc(rt->t,
		 name,
		 &numbofdim,
		 &sizedirentry,
		 &sizedataentry,
		 &sizeinfo,
		 &maxdirfanout,
		 &maxdatafanout,
		 &pagesize,
		 &numbofdirpages,
		 &numbofdatapages,
		 pagesperlevel,
		 &numbofrecords,
		 &height,
		 &unique);

  RSurfaceStats * s = malloc (sizeof (RSurfaceStats));
  s->nlevel = level <= 0 ? height - 1 : (level < height - 1 ? level : height - 1);
  s->nentries = calloc (sizeof (int), s->nlevel);
  s->aspect = malloc (sizeof (double *)*s->nlevel);
  int i;
  for (i = 0; i < s->nlevel; i++)
    s->aspect[i] = malloc (sizeof (double)*pagesperlevel[i]);

  RSurfaceRect rect = {{-0.5,-0.5},{0.5,0.5}};
  RSurfaceSum sum;
  r_surface_sum_init (&sum);
  r_surface_query_region_sum (rt, (RSurfaceCheck) includes, (RSurfaceCheck) intersects, s,
			      rect, &sum);
  for (i = 1; i < s->nlevel; i++)
    qsort (s->aspect[i], s->nentries[i], sizeof (double), compare_double);

  return s;
}

void r_surface_stats_free (RSurfaceStats * s)
{
  int i;
  for (i = 0; i < s->nlevel; i++)
    free (s->aspect[i]);
  free (s->aspect);
  free (s->nentries);
  free (s);
}
