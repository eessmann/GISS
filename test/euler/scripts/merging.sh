#! /bin/sh

usage()
{
    cat <<EOF
Usage: merging.sh [OPTIONS]
Options:
	[--convective]
EOF
    exit $1
}

while test $# -gt 0; do
  case "$1" in
  -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
  *) optarg= ;;
  esac

  case $1 in
    --convective)
      convective="GfsAdvectionParams {
    flux = gfs_face_velocity_convective_flux
  }"
      ;;
    --*)
      usage 1 1>&2
      ;;
  esac
  shift
done

if test -e merging; then 
    rm -r -f merging
fi
mkdir merging
mkdir merging/run1
mkdir merging/run2
mkdir merging/run3

tmax=4.0
cd merging/run1
cat <<EOF > merging.sim
1 0 GfsSimulation GfsBox GfsGEdge {} {
  GfsTime { end = $tmax }
  $convective
  GModule testing
  GfsRefine 6
  GfsRefineSphere 7  { r = 0.25 }
  GfsRefineSphere 8  { r = 0.125 }
  GfsRefineSphere 9  { r = 0.0625 }
  GfsRefineSphere 10 { r = 0.03125 }
  GfsAddGaussianVortex {} { x = 0 y =  0.01 scale = 0.01 }
  GfsAddGaussianVortex {} { x = 0 y = -0.01 scale = 0.01 }
  GfsInitVorticity {}
  GfsOutputTime { istep = 1 } stdout
  GfsOutputScalarNorm { istep = 1 } stdout { v = Divergence }
  GfsOutputProjectionStats { istep = 1 } stdout
  GfsOutputSimulation { step = 0.2 } sim-%3.1f {}
  GfsOutputTiming { start = end } stdout
}
GfsBox {}
EOF
gerris2D merging.sim > log

cd ../run2
cat <<EOF > merging.sim
1 0 GfsSimulation GfsBox GfsGEdge {} {
  GfsTime { end = $tmax }
  $convective
  GModule testing
  GfsRefine 6
  GfsRefineSphere 7  { r = 0.25 }
  GfsRefineSphere 8  { r = 0.125 }
  GfsRefineSphere 9  { r = 0.09375 }
  GfsRefineSphere 10 { r = 0.0625 }
  GfsRefineSphere 11 { r = 0.03125 }
  GfsAddGaussianVortex {} { x = 0 y =  0.01 scale = 0.01 }
  GfsAddGaussianVortex {} { x = 0 y = -0.01 scale = 0.01 }
  GfsInitVorticity {}
  GfsOutputTime { istep = 1 } stdout
  GfsOutputScalarNorm { istep = 1 } stdout { v = Divergence }
  GfsOutputProjectionStats { istep = 1 } stdout
  GfsOutputSimulation { step = 0.2 } sim-%3.1f {}
  GfsOutputTiming { start = end } stdout
}
GfsBox {}
EOF
gerris2D merging.sim > log

cd ../run3
cat <<EOF > merging.sim
1 0 GfsSimulation GfsBox GfsGEdge {} {
  GfsTime { end = $tmax }
  $convective
  GModule testing
  GfsRefine 6
  GfsRefineSphere 7  { r = 0.25 }
  GfsRefineSphere 8  { r = 0.125 }
  GfsRefineSphere 9  { r = 0.09375 }
  GfsRefineSphere 10 { r = 0.0625 }
  GfsRefineSphere 11 { r = 0.046875 }
  GfsRefineSphere 12 { r = 0.03125 }
  GfsAddGaussianVortex {} { x = 0 y =  0.01 scale = 0.01 }
  GfsAddGaussianVortex {} { x = 0 y = -0.01 scale = 0.01 }
  GfsInitVorticity {}
  GfsOutputTime { istep = 1 } stdout
  GfsOutputScalarNorm { istep = 1 } stdout { v = Divergence }
  GfsOutputProjectionStats { istep = 1 } stdout
  GfsOutputSimulation { step = 0.2 } sim-%3.1f {}
  GfsOutputTiming { start = end } stdout
}
GfsBox {}
EOF
gerris2D merging.sim > log
cd ..

uerror=`mktemp /tmp/uerror.XXXXXX`
for time in `seq 0.2 0.2 $tmax`; do
    formatted=`echo $time | awk '{printf ("%3.1f", $1);}'`
    run1=`gfscompare2D -r -n -v run1/sim-$formatted run2/sim-$formatted U 2>&1 | awk '{ if ($2 == "err") printf ("%g %g %g ", $4, $6, $8);}'`
    run1="$run1 `gfscompare2D -r -n -v run2/sim-$formatted run3/sim-$formatted U 2>&1 | awk '{ if ($2 == "err") printf (\"%g %g %g \", $4, $6, $8); }'`"
    echo $time $run1 >> $uerror
done

echo "@description \"merging.sh\""

echo "@WITH G0"
echo "@G0 ON"
echo "@G0 TYPE logy"
echo "@XAXIS LABEL \"Time\""
echo "@xaxis label char size 0.6"
echo "@xaxis ticklabel char size 0.6"
echo "@YAXIS LABEL \"Total error\""
echo "@yaxis label char size 0.6"
echo "@yaxis ticklabel char size 0.6"

echo "@LEGEND STRING 0 \"L1 norm (run1)\""
echo "@TARGET S0"
echo "@TYPE xy" 
awk '{print $1 " " $2}' < $uerror
echo "&"

echo "@LEGEND STRING 1 \"L2 norm (run1)\""
echo "@TARGET S1"
echo "@TYPE xy" 
awk '{print $1 " " $3}' < $uerror
echo "&"

echo "@LEGEND STRING 2 \"Lmax norm (run1)\""
echo "@TARGET S2"
echo "@TYPE xy" 
awk '{print $1 " " $4}' < $uerror
echo "&"

echo "@LEGEND STRING 3 \"L1 norm (run2)\""
echo "@TARGET S3"
echo "@TYPE xy" 
awk '{print $1 " " $8}' < $uerror
echo "&"

echo "@LEGEND STRING 4 \"L2 norm (run2)\""
echo "@TARGET S4"
echo "@TYPE xy" 
awk '{print $1 " " $9}' < $uerror
echo "&"

echo "@LEGEND STRING 5 \"Lmax norm (run2)\""
echo "@TARGET S5"
echo "@TYPE xy" 
awk '{print $1 " " $10}' < $uerror
echo "&"

echo "@WITH G1"
echo "@G1 ON"
echo "@G1 TYPE logy"
echo "@XAXIS LABEL \"Time\""
echo "@xaxis label char size 0.6"
echo "@xaxis ticklabel char size 0.6"
echo "@YAXIS LABEL \"Refined error\""
echo "@yaxis label char size 0.6"
echo "@yaxis ticklabel char size 0.6"

echo "@LEGEND STRING 0 \"L1 norm (run1)\""
echo "@TARGET S0"
echo "@TYPE xy" 
awk '{print $1 " " $5}' < $uerror
echo "&"

echo "@LEGEND STRING 1 \"L2 norm (run1)\""
echo "@TARGET S1"
echo "@TYPE xy" 
awk '{print $1 " " $6}' < $uerror
echo "&"

echo "@LEGEND STRING 2 \"Lmax norm (run1)\""
echo "@TARGET S2"
echo "@TYPE xy" 
awk '{print $1 " " $7}' < $uerror
echo "&"

echo "@LEGEND STRING 3 \"L1 norm (run2)\""
echo "@TARGET S3"
echo "@TYPE xy" 
awk '{print $1 " " $11}' < $uerror
echo "&"

echo "@LEGEND STRING 4 \"L2 norm (run2)\""
echo "@TARGET S4"
echo "@TYPE xy" 
awk '{print $1 " " $12}' < $uerror
echo "&"

echo "@LEGEND STRING 5 \"Lmax norm (run2)\""
echo "@TARGET S5"
echo "@TYPE xy" 
awk '{print $1 " " $13}' < $uerror
echo "&"

echo "@WITH G2"
echo "@G2 ON"
echo "@G2 TYPE xy"
echo "@XAXIS LABEL \"Time\""
echo "@xaxis label char size 0.6"
echo "@xaxis ticklabel char size 0.6"
echo "@YAXIS LABEL \"Total order\""
echo "@yaxis label char size 0.6"
echo "@yaxis ticklabel char size 0.6"

echo "@LEGEND STRING 0 \"L1 norm\""
echo "@TARGET S0"
echo "@TYPE xy" 
awk '{print $1 " " log($2/$8)/log(2)}' < $uerror
echo "&"

echo "@LEGEND STRING 1 \"L2 norm\""
echo "@TARGET S1"
echo "@TYPE xy" 
awk '{print $1 " " log($3/$9)/log(2)}' < $uerror
echo "&"

echo "@LEGEND STRING 2 \"Lmax norm\""
echo "@TARGET S2"
echo "@TYPE xy" 
awk '{print $1 " " log($4/$10)/log(2)}' < $uerror
echo "&"

echo "@WITH G3"
echo "@G3 ON"
echo "@G3 TYPE xy"
echo "@XAXIS LABEL \"Time\""
echo "@xaxis label char size 0.6"
echo "@xaxis ticklabel char size 0.6"
echo "@YAXIS LABEL \"Refined order\""
echo "@yaxis label char size 0.6"
echo "@yaxis ticklabel char size 0.6"

echo "@LEGEND STRING 0 \"L1 norm\""
echo "@TARGET S0"
echo "@TYPE xy" 
awk '{print $1 " " log($5/$11)/log(2)}' < $uerror
echo "&"

echo "@LEGEND STRING 1 \"L2 norm\""
echo "@TARGET S1"
echo "@TYPE xy" 
awk '{print $1 " " log($6/$12)/log(2)}' < $uerror
echo "&"

echo "@LEGEND STRING 2 \"Lmax norm\""
echo "@TARGET S2"
echo "@TYPE xy" 
awk '{print $1 " " log($7/$13)/log(2)}' < $uerror
echo "&"

rm -f $uerror
