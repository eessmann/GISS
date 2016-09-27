 #!/bin/bash
 
module load mpi/openmpi-x86_64

export GISS=/Deployment/GISS 
export PATH=$PATH:$GISS/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GISS/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$GISS/lib/pkgconfig

export ODE=/Deployment/ODE
export PATH=$PATH:$ODE/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ODE/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$ODE/lib/pkgconfig

export HYPRE=/Deployment/hypre
export PATH=$PATH:$HYPRE/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HYPRE/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HYPRE/lib/pkgconfig
