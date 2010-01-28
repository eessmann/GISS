#include <stdio.h>
#include <string.h>
#include "RStarTree/RStarTree.h"
#include "rsurface.h"

/* RSurface */

struct _RSurface {
  RSTREE t;
  char * name;
};

RSurface * r_surface_new (const char * fname, int size, FILE * fp)
{
  return NULL;
}

RSurface * r_surface_open (const char * fname, const char * mode, int size)
{
  RSurface * rt = malloc (sizeof (RSurface));
  if (!strcmp (mode, "w")) {
    RemoveRST (fname);
    if (!CreateRST (fname, size, FALSE)) {
      free (rt);
      return NULL;
    }
  }
  rt->t = NULL;
  if (!OpenRST (&rt->t, fname, !strcmp (mode, "w") ? "rw" : "r")) {
    free (rt);
    return NULL;
  }
  rt->name = malloc (sizeof (char)*(strlen (fname) + 1));
  strcpy (rt->name, fname);
  return rt;
}

void r_surface_close (RSurface * rt)
{
  CloseRST (&rt->t);
  free (rt->name);
  free (rt);
}

int r_surface_insert  (RSurface * rt, double p[3], int id)
{
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
  typrect rect, unused;
  rect[0].l = min[0]; rect[0].h = max[0];
  rect[1].l = min[1]; rect[1].h = max[1];
  void * data[2];
  data[0] = q;
  data[1] = user_data;
  RegionQuery (rt->t, rect, unused, Intersects, Intersects, ManageQuery, data);
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
  RegionQueryInfo (rt->t, (Check) includes, (Check) intersects, data, (typinterval *) rect,
		   (typdirinfo *) sum);
}

const char * r_surface_name (RSurface * rt)
{
  return rt->name;
}

void r_surface_update (RSurface * rt)
{
  Update (rt->t);
}

int r_surface_depth (RSurface * rt)
{
  int height = -1;
  GetHeight (rt->t, &height);
  return height;
}

void r_surface_info (RSurface * rt)
{
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
