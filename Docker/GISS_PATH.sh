 #!/bin/bash
export GHC=/Deployment/GHC 
export PATH=$PATH:$GHC/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GHC/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$GHC/lib/pkgconfig

export ODE=/Deployment/ODE
export PATH=$PATH:$ODE/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ODE/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$ODE/lib/pkgconfig

export HYPRE=/Deployment/hypre
export PATH=$PATH:$HYPRE/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HYPRE/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HYPRE/lib/pkgconfig
