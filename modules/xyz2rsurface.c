#include <stdio.h>
#include <stdlib.h>
#include "rsurface.h"

int main (int argc, char** argv)
{
  if (argc != 3) {
    fprintf (stderr, "Usage: %s basename pagesize\n", argv[0]);
    return -1;
  }

  RSurface * rs = r_surface_open (argv[1], "w", atoi (argv[2]));
  if (rs == NULL) {
    fprintf (stderr, "xyz2rsurface: cannot open files `%s*'\n", argv[1]);
    return 1;
  }
  double x[3];
  int id, count = 0;
  while (scanf ("%d %lf %lf %lf", &id, &x[0], &x[1], &x[2]) == 4) {
    if (!r_surface_insert (rs, x, id)) {
      fprintf (stderr, "\nxyz2rsurface: error inserting point (%g,%g,%g)\n",
	       x[0], x[1], x[2]);
      return 1;
    }
    count++;
  }
  r_surface_close (rs);

  return 0.;
}
