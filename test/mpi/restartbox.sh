#! /bin/sh
# Tests that a simulation can be correctly restarted from a saved 
# simulation file (with a refined box inside)

modules=$HOME/adds/lib/gerris
output=`mktemp /tmp/output.XXXXXX`
sim=`mktemp /tmp/sim.XXXXXX`
input=`mktemp /tmp/input.XXXXXX`
cat > $input <<EOF
16 32 FlSimulation FlBox FlGEdge { rootlevel = 2 x = -0.375 y = -0.375 } {
  FlTime { i = 0 t = 0 end = 0.2 }
  FlAdvectionParams {
  cfl      = 0.8
  gradient = fl_center_gradient
  flux     = fl_face_velocity_advection_flux
}
  FlTracerAdvectionParams {
  cfl      = 0.8
  gradient = fl_center_van_leer_gradient
  flux     = fl_face_advection_flux
}
  FlApproxProjectionParams {
  tolerance = 0.001
  nrelax    = 4
  minlevel  = 0
  nitermax  = 100
}
  FlProjectionParams {
  tolerance = 0.001
  nrelax    = 4
  minlevel  = 0
  nitermax  = 100
}
  GModule $modules/libperiodic_flow.so
  GModule $modules/libtesting.so
  FlRefine 5
  FlRefineBox 6
  FlInitNonStationary { }
  FlOutputSimulation { start = 0.1 } $sim { variables = P,C,U,V }
  FlOutputNonStationaryError { start = 0.1 step = 0.01 } stdout
}
FlBox { id = 12 pid = 11 size = 256 } NULL NULL NULL NULL
FlBox { id = 2 pid = 0 size = 256 } NULL NULL NULL NULL
FlBox { id = 7 pid = 5 size = 256 } NULL NULL NULL NULL
FlBox { id = 11 pid = 10 size = 64 } NULL NULL NULL NULL
FlBox { id = 6 pid = 4 size = 64 } NULL NULL NULL NULL
FlBox { id = 10 pid = 9 size = 64 } NULL NULL NULL NULL
FlBox { id = 15 pid = 14 size = 256 } NULL NULL NULL NULL
FlBox { id = 5 pid = 3 size = 64 } NULL NULL NULL NULL
FlBox { id = 16 pid = 15 size = 64 } NULL NULL NULL NULL
FlBox { id = 9 pid = 8 size = 64 } NULL NULL NULL NULL
FlBox { id = 14 pid = 13 size = 64 } NULL NULL NULL NULL
FlBox { id = 4 pid = 2 size = 64 } NULL NULL NULL NULL
FlBox { id = 8 pid = 7 size = 64 } NULL NULL NULL NULL
FlBox { id = 13 pid = 12 size = 64 } NULL NULL NULL NULL
FlBox { id = 3 pid = 1 size = 64 } NULL NULL NULL NULL
FlBox { id = 1 pid = 6 size = 64 } NULL NULL NULL NULL
1 3 3
4 1 0
6 1 3
7 1 1
2 12 3
2 15 0
3 2 0
7 2 3
3 13 3
5 3 0
4 5 3
10 4 3
9 4 0
5 16 3
5 15 1
6 13 2
10 6 0
14 6 1
7 9 0
14 7 3
12 8 0
15 8 3
16 8 1
11 8 2
9 15 3
11 9 3
10 16 2
11 10 0
14 11 0
13 12 0
14 12 2
16 13 0
EOF
../../src/gerris $input > $output

output_restart=`mktemp /tmp/output.XXXXXX`
../../src/gerris $sim > $output_restart

rm -f $input $sim

if ./compare.sh $output_restart $output 0.15; then
    rm -f $output $output_restart
    exit 1
fi
rm -f $output $output_restart

