#! /bin/sh
# Tests that a simulation can be correctly restarted from a saved 
# simulation file

modules=$HOME/adds/lib/gerris
output=`mktemp /tmp/output.XXXXXX`
sim=`mktemp /tmp/sim.XXXXXX`
input=`mktemp /tmp/input.XXXXXX`
cat > $input <<EOF
1 2 FlDomain FlBox FlGEdge {} {
  FlApproxProjectionParams {
    tolerance = 1e-3
    nrelax    = 4
    minlevel  = 0
    nitermax  = 100
  }
  FlProjectionParams {
    tolerance = 1e-3
    nrelax    = 4
    minlevel  = 0
    nitermax  = 100
  }
  FlTime { end = 0.2 }
  FlRefine 6
  GModule $modules/libperiodic_flow.so
  FlInitNonStationary {}
  FlOutputSimulation { start = 0.1 } $sim {}
  FlOutputNonStationaryError { start = 0.1 step = 0.01 } stdout
}
FlBox {}
1 1 0
1 1 2
EOF
../../src/gerris $input > $output

output_restart=`mktemp /tmp/output.XXXXXX`
../../src/gerris $sim > $output_restart

rm -f $input $sim

if ./compare.sh $output_restart $output 0.1; then
    rm -f $output $output_restart
    exit 1
fi
rm -f $output $output_restart

