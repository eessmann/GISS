levels="5 6 7 8"

if ! $donotrun; then
    rm -f fit
    for level in $levels; do
	if gerris2D -D LEVEL=$level oscillation.gfs >> fit; then :
	else
	    exit 1
	fi
    done
fi

rm -f fit-*
if awk '{
  level = $1; a = $2; b = $3; c = $4;
  for (t = 0; t <= 1.; t += 0.005)
    print t, 2.*a*exp(-b*t) >> "fit-" level;
}' < fit ; then :
else
    exit 1
fi

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20

    D = 0.2
    n = 2.
    sigma = 1.
    rhol = 1.
    rhog = 1./1000.
    r0 = 0.1
    omega0 = sqrt((n**3-n)*sigma/((rhol+rhog)*r0**3))

    set output 'k.eps'
    set xlabel 'Time'
    set ylabel 'Kinetic energy'
    set logscale y
    plot [0:1][8e-5:]'k-8' u 3:5 t "256x256" w l, 'k-7' u 3:5 t "128x128" w l, 'k-6' u 3:5 t "64x64" w l, 'k-5' u 3:5 t "32x32" w l, 'fit-8' t "fit" w l lt 7, 'fit-7' t "" w l lt 7, 'fit-6' t "" w l lt 7, 'fit-5' t "" w l lt 7

    set output 'laplace.eps'
    set xlabel 'Diameter (grid points)'
    set ylabel 'Equivalent Laplace number'
    set logscale y
    set logscale x 2
    set grid
    empirical_constant = 30.
    plot 'fit' u (D*2.**(\$1)):(1./(\$3**2.*D**3.))*empirical_constant**2. t "" w lp pt 5 ps 2

    unset logscale
    set output 'frequency.eps'
    set xlabel 'Diameter (grid points)'
    set ylabel 'Frequency error (%)'
    unset grid
    set xzeroaxis
    plot 'fit' u (D*2.**(\$1)):(\$4/2./omega0-1.)*100. t "" w lp pt 5 ps 2

    
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('fit',1,3) - Curve('fit.ref',1,3)).max() > 1e-2 or\
   (Curve('fit',1,4) - Curve('fit.ref',1,4)).max() > 1e-2:
  exit(1)
EOF
else
   exit 1
fi
