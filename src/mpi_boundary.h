/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001 National Institute of Water and Atmospheric Research
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

#ifndef __MPI_BOUNDARY_H__
#define __MPI_BOUNDARY_H__

#include <mpi.h>
#include "boundary.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef struct _GfsBoundaryMpi         GfsBoundaryMpi;
typedef struct _GfsBoundaryMpiClass    GfsBoundaryMpiClass;

struct _GfsBoundaryMpi {
  /*< private >*/
  GfsBoundary parent;

  MPI_Comm comm;
  gint process, id;

  MPI_Request request[2];
  guint nrequest;

  GArray * sndbuf, * rcvbuf;
  unsigned int sndcount, rcvcount;
};

struct _GfsBoundaryMpiClass {
  GfsBoundaryClass parent_class;
};

#define GFS_BOUNDARY_MPI(obj)            GTS_OBJECT_CAST (obj,\
					           GfsBoundaryMpi,\
					           gfs_boundary_mpi_class ())
#define GFS_BOUNDARY_MPI_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						       GfsBoundaryMpiClass,\
						       gfs_boundary_mpi_class())
#define GFS_IS_BOUNDARY_MPI(obj)         (gts_object_is_from_class (obj,\
						   gfs_boundary_mpi_class ()))
     
GfsBoundaryMpiClass * gfs_boundary_mpi_class    (void);
GfsBoundaryMpi *      gfs_boundary_mpi_new      (GfsBoundaryMpiClass * klass,
					       GfsBox * box,
					       FttDirection d,
					       gint process,
					       gint id);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MPI_BOUNDARY_H__ */
