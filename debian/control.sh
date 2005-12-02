GTS_VERSION=`pkg-config gts --modversion`

cat <<EOF > debian/control
Source: gerris-snapshot
Section: math
Priority: extra
Maintainer: Stephane Popinet <popinet@users.sf.net>
Build-Depends: debhelper (>> 4), libgts-snapshot-dev (>= $GTS_VERSION), libglib2.0-dev, autotools-dev

Package: gerris-snapshot
Section: math
Priority: extra
Architecture: any
Depends: \${shlibs:Depends}, libgts-snapshot (>= $GTS_VERSION)
Conflicts: gerris
Replaces: gerris
Description: Gerris Flow Solver (development snapshot)
 Gerris is a system for the solution of the partial differential
 equations describing fluid flow.
 .
 A brief summary of its main (current) features:
 .
    * Quadtree-based (Octree in 3D) spatial discretisation with
      automatic and dynamic refinement.
    * Multigrid Poisson solver.
    * Second-order Godunov type advection scheme.
    * Solves the time-dependent incompressible Euler, Stokes ans Navier-Stokes
      equations.
    * Support for complex solid boundaries (automatic locally-refined
      mesh generation).
 .
 See http://gfs.sf.net for more information and documentation.
EOF
