 #!/bin/bash
module avail
module load mpi/openmpi-x86_64
module list

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


cd /GISS/Third_Party/ode
bash bootstrap
bash configure --enable-shared=yes --enable-double-precision --enable-builtin-threading-impl --enable-libccd --prefix=/Deployment/ODE
make -j6
make install

cd /GISS/Third_Party/hypre/src
bash configure --prefix=$HYPRE --enable-shared
make -j6 
make install

cd /GISS/Third_Party/gts
bash configure --prefix=/Deployment/GISS
make -j6
make install

cd /GISS/Third_Party/gerris_Inuse
/sbin/ldconfig -v
CPPFLAGS=-I$HYPRE/include LDFLAGS="-L$HYPRE/lib" ./configure --prefix=$GISS
make -j6
make install
