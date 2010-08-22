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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "kdt.h"

int main (int argc, char * argv[])
{
  if (argc != 4) {
    fprintf (stderr, "usage: xyz2kdt basename blksize input\n");
    return 1;
  }

  int fd = open (argv[3], O_RDWR);
  if (fd < 0) {
    perror ("xyz2rsurface: could not open input file");
    return 1;
  }
  remove (argv[3]);

  Kdt * kdt = kdt_new ();
  kdt_create_from_file (kdt, argv[1], atoi (argv[2]), fd);
  kdt_destroy (kdt);
  close (fd);

  return 0;
}
