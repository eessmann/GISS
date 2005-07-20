#!/bin/sh

if test -f wave.gts; then
    :
else
    file=`mktemp /tmp/wave.XXXXXX`
    awk 'BEGIN{
      print "0.51 -0.51"
      for (x = 0.51; x >= -0.51; x -= 0.01) 
        print x " " 0.01*cos (2.*3.14159265359*x);
      print "-0.51 -0.51"
    }' > $file
    shapes $file > wave.gts
    rm -f $file
fi

param=`mktemp /tmp/wave.XXXXXX`

for level in 3 4 5 6 7; do
    rm -f wave-$level
    cat <<EOF > $param
1 1 GfsSimulation GfsBox GfsGEdge {} {
  Time {
    end = 2.2426211256
    dtmax = 1.01430173026e-3
  }
  ApproxProjectionParams { tolerance = 1e-6 }
  ProjectionParams { tolerance = 1e-6 }
  Refine $level
  VariableTracer {} T { scheme = vof }
  SourceTension {} U V T 1
  AdvectionParams { scheme = none }
  SourceDiffusion {} U 0.0182571749236
  SourceDiffusion {} V 0.0182571749236
  InitFraction {} T wave.gts
  OutputProgress { istep = 1 } stdout
  OutputSimulation { istep = 3 } $param-%ld.gfs { binary = 1 }
  EventScript { istep = 3 } {
    echo \$GfsTime | awk -v level=$level '{printf ("%g ", \$1*11.1366559937) >> ("wave-" level) }'
    gfs2oogl2D -g -c T < $param-\$GfsIter.gfs | sort -n +1 -2 | awk -v level=$level 'BEGIN {
      n = exp (level*log (2.));
      x = 1./(n + 1.);
      c1 = c2 = 0.;
    }
    {
      if (\$1 > 0. && \$1 < x) {
        if (c1 + c2 < 1. && c2 + \$4 > 1.) {
          print y2+(0.5-c2)/n >> ("wave-" level);
          exit 0;
        }
        c1 = c2;
        c2 = \$4;
        y2 = \$2;
      }
    }'
    rm -f $param-\$GfsIter.gfs
  }
}
GfsBox {}
1 1 right 
EOF
    if gerris2D $param; then
	if awk -v level=$level 'BEGIN {s = 0.; n = 0; } {
          t = $1; y = $2;
          getline < "prosperetti"
          dt = abs ($1 - t);
          if (dt > 1e-5) {
            print "Level " level ": times do not match " $1 " " t " " dt
            exit 1
          }
          s += (y - $2)*(y - $2);
          n += 1;
        }
        END {
          s = 100.*sqrt (s/n)/0.01;
          printf ("Level %d: %g%%: ", level, s);

          max[3] = 43.; max[4] = 26.; max[5] = 8.; max[6] = 1.9; max[7] = 0.7;
          if (s > max[level]) {
            print "FAIL [" max[level] "]"
            exit 1
          }
          else
            print "PASS";
        }' < wave-$level; then
	    :
	else
	    exit 1;
        fi
    else
	exit 1;
    fi
done

rm -f $param
