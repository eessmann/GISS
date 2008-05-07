levels="3 4 5 6"

if ! $donotrun; then
    rm -f error
    for level in $levels; do
	if gerris2D -DLEVEL=$level poiseuille.gfs >> error; then :
	else
	    exit 1
	fi
    done
fi

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20
    set output 'convergence.eps'
    set xlabel 'Number of grid points'
    set ylabel 'Maximum error'
    set logscale
    set grid
    set xtics 2
    plot [6:80]'error' u (2**\$1):4 w lp t 'Gerris' pt 5, 1./x**2./5. t 'second order'
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('error',1,4) - Curve('error.ref',1,4)).max() > 1e-6:
    exit(1)
EOF
else
   exit 1
fi
