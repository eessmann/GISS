if test x$donotrun != xtrue; then
    rm -f p U V
    for L in `seq 9`; do
	if gerris2D -DLEVEL=$L $1 ; then :
	else
	    exit 1
	fi
    done
    if echo Save solution.eps '{ format = EPS line_width = 0.25 }' | gfsview-batch2D solution.gfs `basename $1 .gfs`.gfv ; then :
    else
	exit 1
    fi
else
    exit 1
fi
if cat <<EOF | gnuplot ; then :
    set title 'Convergence for groundwater flow with uniform permeability'
    set xlabel 'Refinement level'
    set ylabel 'Error norm'
    set logscale y
    set terminal postscript eps color lw 3 solid 20
    set output 'convergence.eps'
    set style data lp
    plot 'p' t '1 (p)', \
        '' u 1:3 t '2 (p)', \
        '' u 1:4 t 'max (p)', \
        'U' t '1 (U)' w l lt 4, \
        '' u 1:3 t '2 (U)' w l lt 5, \
        '' u 1:4 t 'max (U)' w l lt 6, \
        'V' t '1 (V)' w p lt 4, \
        '' u 1:3 t '2 (V)' w p lt 5, \
        '' u 1:4 t 'max (V)' w p lt 6
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('p',1,2) - Curve('p.ref',1,2)).max() > 1e-8 or\
   (Curve('p',1,3) - Curve('p.ref',1,3)).max() > 1e-6:
    print (Curve('p',1,2) - Curve('p.ref',1,2)).max()
    print (Curve('p',1,3) - Curve('p.ref',1,3)).max()
    exit(1)
if (Curve('U',1,2) - Curve('U.ref',1,2)).max() > 1e-8 or\
   (Curve('U',1,3) - Curve('U.ref',1,3)).max() > 1e-6:
    print (Curve('U',1,2) - Curve('U.ref',1,2)).max()
    exit(1)
EOF
else
    exit 1
fi
