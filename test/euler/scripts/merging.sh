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
      convective="AdvectionParams {
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
cat <<EOF > merging.gfs
1 0 GfsSimulation GfsBox GfsGEdge {} {
  Time { end = $tmax }
  $convective
  Refine {
    double r = sqrt (x*x + y*y); 
    return r < 0.03125 ? 10 : r < 0.0625 ? 9 : r < 0.125 ? 8 : r < 0.25 ? 7 : 6;
  }
  InitVorticity {} {
    double gaussian (double xo, double yo, double scale) {
      return 2.*M_PI*exp (- 2.*((x - xo)*(x - xo) + (y - yo)*(y - yo))/(scale*scale));
    }
    return gaussian (0, 0.01, 0.01) + gaussian (0, -0.01, 0.01);
  }
  OutputTime { istep = 1 } stdout
  OutputScalarNorm { istep = 1 } stdout { v = Divergence }
  OutputProjectionStats { istep = 1 } stdout
  OutputSimulation { step = 0.2 } sim-%3.1f {}
  OutputTiming { start = end } stdout
}
GfsBox {}
EOF
gerris2D merging.gfs > log

cd ../run2
cat <<EOF > merging.gfs
1 0 GfsSimulation GfsBox GfsGEdge {} {
  Time { end = $tmax }
  $convective
  Refine {
    double r = sqrt (x*x + y*y); 
    return r < 0.03125 ? 11 : r < 0.0625 ? 10 : r < 0.09375 ? 9 : r < 0.125 ? 8 : r < 0.25 ? 7 : 6;
  }
  InitVorticity {} {
    double gaussian (double xo, double yo, double scale) {
      return 2.*M_PI*exp (- 2.*((x - xo)*(x - xo) + (y - yo)*(y - yo))/(scale*scale));
    }
    return gaussian (0, 0.01, 0.01) + gaussian (0, -0.01, 0.01);
  }
  OutputTime { istep = 1 } stdout
  OutputScalarNorm { istep = 1 } stdout { v = Divergence }
  OutputProjectionStats { istep = 1 } stdout
  OutputSimulation { step = 0.2 } sim-%3.1f {}
  OutputTiming { start = end } stdout
}
GfsBox {}
EOF
gerris2D merging.gfs > log

cd ../run3
cat <<EOF > merging.gfs
1 0 GfsSimulation GfsBox GfsGEdge {} {
  Time { end = $tmax }
  $convective
  Refine {
    double r = sqrt (x*x + y*y); 
    return r < 0.03125 ? 12 : r < 0.046875 ? 11 : r < 0.0625 ? 10 : r < 0.09375 ? 9 : r < 0.125 ? 8 : r < 0.25 ? 7 : 6;
  }
  InitVorticity {} {
    double gaussian (double xo, double yo, double scale) {
      return 2.*M_PI*exp (- 2.*((x - xo)*(x - xo) + (y - yo)*(y - yo))/(scale*scale));
    }
    return gaussian (0, 0.01, 0.01) + gaussian (0, -0.01, 0.01);
  }
  OutputTime { istep = 1 } stdout
  OutputScalarNorm { istep = 1 } stdout { v = Divergence }
  OutputProjectionStats { istep = 1 } stdout
  OutputSimulation { step = 0.2 } sim-%3.1f {}
  OutputTiming { start = end } stdout
}
GfsBox {}
EOF
gerris2D merging.gfs > log
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
