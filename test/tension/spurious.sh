#!/bin/sh

if test -f circle.gts; then
    :
else
    ../poisson/shapes ellipse | transform --scale 0.8 > circle.gts
fi

param=`mktemp /tmp/spurious.XXXXXX`

for adapt in "" "AdaptFunction { istep = 1 } { cmax = 0.1 maxlevel = 5 } { return _Tx*_Tx + _Ty*_Ty; }"; do
  if test "$adapt"; then
    echo "Adaptive refinement"
  else
    echo "Constant refinement"
  fi
  for La in 12000 1200 120; do
      mu=`echo $La | awk '{print sqrt (0.4/$1)}'`
      tmax=`echo $mu | awk '{print 1000.*$1*0.4}'`
      du=`echo $mu | awk '{print 1e-7/$1}'`
      cat <<EOF > $param
1 2 GfsSimulation GfsBox GfsGEdge {} {
  Time {
    end = $tmax
    dtmax = 1e-3
  }
  ApproxProjectionParams { tolerance = 1e-6 }
  ProjectionParams { tolerance = 1e-6 }
  Refine 5
  VariableTracer {} T { scheme = vof }
  SourceTension {} U V T 1
  AdvectionParams { scheme = none }
  SourceDiffusion {} U $mu
  SourceDiffusion {} V $mu
  InitFraction {} T circle.gts
  $adapt
  EventStop { istep = 1 } U $du
  OutputScalarStats {istep = 1} spurious-$La { v = P }
  OutputScalarNorm {istep = 1} spurious-$La { v = Velocity }
  OutputScalarSum { istep = 1 } spurious-$La { v = T }
}
GfsBox {}
1 1 right 
1 1 top
EOF
      gerris2D $param
      if awk -v La=$La -v mu=$mu '{if ($1 == "Velocity") {max1 = max; max = $9;}} END {
      printf ("Laplace number: %6g max(U)*mu/sigma: %g\n", La, max1*mu);
      if (max1*mu > 2.5e-3)
        exit 1;
    }' < spurious-$La; then
	  :
      else
	  rm -f $param
	  exit 1
      fi
  done
done

rm -f $param
