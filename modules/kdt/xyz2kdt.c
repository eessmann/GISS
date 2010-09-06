/* Gerris - The GNU Flow Solver
 * Copyright (C) 2010 National Institute of Water and Atmospheric Research
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.  
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>

#include "kdt.h"

static void progress (float complete, void * data)
{
  struct timeval * start = data, now;
  gettimeofday (&now, NULL);
  double remaining = ((double) (now.tv_usec - start->tv_usec) + 
		      1e6*(double) (now.tv_sec - start->tv_sec))*(1. - complete);
  int hours = remaining/1e6/3600.;
  int mins = remaining/1e6/60. - 60.*hours;
  int secs = remaining/1e6 - 3600.*hours - 60.*mins;
  fprintf (stderr, "\rxyz2kdt: %3.0f%% complete %02d:%02d:%02d remaining", 
	   complete*100.,
	   hours, mins, secs);
}

int main (int argc, char * argv[])
{
  if (argc != 3) {
    fprintf (stderr, "usage: xyz2kdt basename blksize\n");
    return 1;
  }

  char name[] = "XXXXXX";
  int fd = mkstemp (name);
  assert (fd >= 0);
  assert (unlink (name) == 0);

  fprintf (stderr, "xyz2kdt: reading points...");
  int blksize = atoi (argv[2]);
  KdtHeap h;
  kdt_heap_create (&h, fd, blksize/sizeof (KdtPoint));
  KdtPoint p;
  while (scanf ("%lf %lf %lf", &p.x, &p.y, &p.z) == 3)
    kdt_heap_put (&h, &p);
  kdt_heap_flush (&h);
  kdt_heap_free (&h);
  lseek (fd, 0, SEEK_SET);

  struct timeval start;
  gettimeofday (&start, NULL);
  Kdt * kdt = kdt_new ();
  kdt_create_from_file (kdt, argv[1], blksize, fd,
			progress, &start);
  fprintf (stderr, "\n");
  kdt_destroy (kdt);
  close (fd);

  return 0;
}
