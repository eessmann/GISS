#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */

#include "ftt.h"
#include "rsurface.h"

typedef struct {
  FILE ** fp;
  int * size;
  char * name, * ext;
  int maxdepth;
} Params;

static int includes (RSurfaceRect RSTrect, Params * p, int depth)
{
  if (RSTrect[0].l == RSTrect[0].h && RSTrect[1].l == RSTrect[1].h)
    p->size[depth]++;
  return 0;
}

#define MAXEXT 10

static int intersects (RSurfaceRect RSTrect, Params * p, int depth)
{
  if (p->fp[depth] == NULL) {
    snprintf (p->ext, MAXEXT, "-l%d", depth);
    p->fp[depth] = fopen (p->name, "w");
  }
  fprintf (p->fp[depth], "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
	   RSTrect[0].l, RSTrect[1].l,
	   RSTrect[0].h, RSTrect[1].l,
	   RSTrect[0].h, RSTrect[1].h,
	   RSTrect[0].l, RSTrect[1].h,
	   RSTrect[0].l, RSTrect[1].l);
  p->size[depth]++;
  return (depth < p->maxdepth);
}

int main (int argc, char** argv)
{
  int c = 0;
  int verbose = 0, gnuplot = 0;

  /* parse options using getopt */
  while (c != EOF) {
#ifdef HAVE_GETOPT_LONG
    static struct option long_options[] = {
      {"gnuplot", no_argument, NULL, 'g'},
      {"verbose", no_argument, NULL, 'v'},
      {"help", no_argument, NULL, 'h'},
      { NULL }
    };
    int option_index = 0;
    switch ((c = getopt_long (argc, argv, "vhg",
			      long_options, &option_index))) {
#else /* not HAVE_GETOPT_LONG */
    switch ((c = getopt (argc, argv, "vhg"))) {
#endif /* not HAVE_GETOPT_LONG */
    case 'g': /* gnuplot */
      gnuplot = 1;
      break;
    case 'v': /* verbose */
      verbose = 1;
      break;
    case 'h': /* help */
      fprintf (stderr,
	       "Usage: rsurfacedraw [OPTION] BASENAME MAXDEPTH\n"
	       "\n"
	       "Draws gnuplot representation of the bounding boxes of the R*-tree.\n"
	       "\n"
	       "  -g    --gnuplot     write gnuplot commands on standard output\n"
	       "  -v    --verbose     display statistics\n"
	       "  -h    --help        display this help and exit\n"
	       "\n"
	       "Report bugs to %s\n",
	       FTT_MAINTAINER);
      return 0; /* success */
      break;
    case '?': /* wrong options */
      fprintf (stderr, "Try `rsurfacedraw --help' for more information.\n");
      return 1; /* failure */
    }
  }
  if (optind >= argc) { /* missing BASENAME */
    fprintf (stderr, 
	     "rsurfacedraw: missing BASENAME\n"
	     "Try `rsurfacedraw --help' for more information.\n");
    return 1; /* failure */
  }
  if (optind + 1 >= argc) { /* missing MAXDEPTH */
    fprintf (stderr, 
	     "rsurfacedraw: missing MAXDEPTH\n"
	     "Try `rsurfacedraw --help' for more information.\n");
    return 1; /* failure */
  }

  RSurface * rs = r_surface_open (argv[optind], "r", 0);
  RSurfaceRect rect = {{-0.5,-0.5},{0.5,0.5}};
  RSurfaceSum s;
  Params p;
  p.maxdepth = atoi (argv[optind + 1]);
  p.name = malloc (sizeof (char)*(strlen (argv[optind]) + MAXEXT + 1));
  strcpy (p.name, argv[optind]);
  p.ext = &p.name[strlen (p.name)];
  p.fp = calloc (p.maxdepth + 1, sizeof (FILE *));
  p.size = calloc (p.maxdepth + 1, sizeof (int));
  r_surface_sum_init (&s);
  r_surface_query_region_sum (rs, (RSurfaceCheck) includes, (RSurfaceCheck) intersects, &p, 
			      rect, &s);
  r_surface_close (rs);

  if (verbose) {
    int i;
    for (i = 1; i <= p.maxdepth; i++)
      if (p.size[i] > 0) {
	fprintf (stderr, "level %d: %d", i, p.size[i]);
	if (i < p.maxdepth && p.size[i + 1] > 0)
	  fprintf (stderr, " average # of entries: %d\n", p.size[i + 1]/p.size[i]);
	else
	  fputc ('\n', stderr);
      }
  }
  
  if (gnuplot) {
    int i;
    printf ("set size ratio -1\n");
    printf ("plot '%s-l1' w l t 'Level 1'", argv[optind]);
    for (i = 2; i <= p.maxdepth; i++)
      if (p.size[i] > 0)
	printf (", '%s-l%d' w l t 'Level %d'", argv[optind], i, i);
    putchar ('\n');
  }

  return 0.;
}
