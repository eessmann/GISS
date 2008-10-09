#include <stdio.h>
#include <stdlib.h>
#include "rsurface.h"

static int query (double p[3], void * data)
{
  printf ("%.8f %.8f %g\n", p[0], p[1], p[2]);
  return 0;
}

int main (int argc, char** argv)
{
  if (argc != 2) {
    fprintf (stderr, "Usage: %s basename\n", argv[0]);
    return -1;
  }

  RSurface * rs = r_surface_open (argv[1], "r", 0);
  double min[2], max[2];
  int count = 0;
  while (scanf ("%lf %lf %lf %lf", &min[0], &min[1], &max[0], &max[1]) == 4) {
    r_surface_query_region (rs, min, max, query, NULL);
    if (count > 0 && count % 1000 == 0)
      fprintf (stderr, "\r%d", count);
    count++;
  }
  if (count >= 1000)
    fputc ('\n', stderr);
  r_surface_close (rs);

  return 0.;
}
