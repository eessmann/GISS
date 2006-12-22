levels="3 4 5 6"

if ! $donotrun; then
    if test -f wave.gts; then :
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

    for level in $levels; do
	rm -f wave-$level
	if sed "s/LEVEL/$level/g" < $1 | gerris2D -; then
	    :
	else
	    exit 1;
	fi
    done
fi

rm -f convergence
for level in $levels; do
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
          s = sqrt (s/n)/0.01;
          printf ("%d %g\n", level, s);
        }' < wave-$level >> convergence; then
	:
    else
	exit 1;
    fi
done

awk 'BEGIN{first=1}{ 
  if (first) printf("Gerris & %.5f",$2);
  else printf(" & %.5f",$2);
  first=0;
}
END{printf("\n");}' < convergence > convergence.tex

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20
    set output 'amplitude.eps'
    set xlabel 'tau'
    set ylabel 'Relative amplitude'
    plot 'prosperetti' w l t "Prosperetti", 'wave-6' every 10 w p t "Gerris"
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('convergence',1,2) - Curve('convergence.ref',1,2)).max() > 1e-8:
    exit(1)
EOF
else
   exit 1
fi
