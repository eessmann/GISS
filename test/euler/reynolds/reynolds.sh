rm -f reynolds

for level in 5 6 7; do
  if sed "s/LEVEL/$level/g" < $1 | gerris2D - | awk -v m=$2 -v level=$level '{
    time = $3
    ke = $5
    if (time == 0)
      ke0 = ke;
  }END{
    a = -log(ke/ke0)/time
    nu = a/(4.*(2.*m*3.14159265359)^2)
    print level " " 1./nu
  }' >> reynolds; then :
  else
      exit 1
  fi
done

if cat <<EOF | gnuplot ; then :
    set term postscript eps color solid 20
    set output 'divmax.eps'
    set xlabel 'Time'
    set ylabel 'Divergence Max'
    plot [0:2]'div5' u 3:9 t "5" w l lw 2, 'div6' u 3:9 t "6" w l lw 2, 'div7' u 3:9 t "7" w l lw 2, 'div5.ref' u 3:9 t "5 (ref)" w l lw 2, 'div6.ref' u 3:9 t "6 (ref)" w l lw 2, 'div7.ref' u 3:9 t "7 (ref)" w l lw 2
    set output 'divL2.eps'
    set ylabel 'Divergence L2'
    plot [0:2]'div5' u 3:7 t "5" w l lw 2, 'div6' u 3:7 t "6" w l lw 2, 'div7' u 3:7 t "7" w l lw 2, 'div5.ref' u 3:7 t "5 (ref)" w l lw 2, 'div6.ref' u 3:7 t "6 (ref)" w l lw 2, 'div7.ref' u 3:7 t "7 (ref)" w l lw 2
    set output 'kinetic.eps'
    set ylabel 'Kinetic energy'
    plot [0:2]'kinetic5' u 3:5 t "5" w l lw 2, 'kinetic6' u 3:5 t "6" w l lw 2, 'kinetic7' u 3:5 t "7" w l lw 2
    set output 'reynolds.eps'
    set xlabel 'Level'
    set ylabel 'Effective Reynolds number'
    set logscale y
    plot 'reynolds' u 1:2 t "" w lp lw 2, 'reynolds.ref' u 1:2 t "ref" w lp lw 2
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('div5',3,9) - Curve('div5.ref',3,9)).max() > 0.05*Curve('div5.ref',3,9).mean() or\
   (Curve('div6',3,9) - Curve('div6.ref',3,9)).max() > 0.05*Curve('div6.ref',3,9).mean() or\
   (Curve('div7',3,9) - Curve('div7.ref',3,9)).max() > 0.05*Curve('div7.ref',3,9).mean() or\
   (Curve('div5',3,7) - Curve('div5.ref',3,7)).max() > 0.05*Curve('div5.ref',3,7).mean() or\
   (Curve('div6',3,7) - Curve('div6.ref',3,7)).max() > 0.05*Curve('div6.ref',3,7).mean() or\
   (Curve('div7',3,7) - Curve('div7.ref',3,7)).max() > 0.05*Curve('div7.ref',3,7).mean() or\
   (Curve('reynolds',1,2) - Curve('reynolds.ref',1,2)).max() > 0.0 :
   exit(1)
EOF
else
   exit 1
fi
