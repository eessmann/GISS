if ! $donotrun; then
    shapes ellipse | transform --scale 1.6 | transform --tx -0.5 --ty 0.5 > circle.gts
    for La in 12000 1200 120; do
	mu=`echo $La | awk '{print sqrt (0.8/$1)}'`
	tmax=`echo $mu | awk '{print 1000.*$1*0.8}'`
	if sed "s/LEVEL/6/g" < $1 |\
           sed "s/MU/$mu/g" |\
           sed "s/TMAX/$tmax/g" | gerris2D - > La\=$La; then :
	else
	    exit 1
	fi
    done

    rm -f convergence
    for level in 4 5 6; do
        La=12000
        mu=`echo $La | awk '{print sqrt (0.8/$1)}'`
	tmax=`echo $mu | awk '{print 1000.*$1*0.8}'`
	if sed "s/LEVEL/$level/g" < $1 |\
           sed "s/MU/$mu/g" |\
           sed "s/TMAX/$tmax/g" | gerris2D - | awk -v level=$level '{
             max = $2
           }END{print 0.8*2**level " " max}' >> convergence; then :
	else
	    exit 1
	fi
    done
fi

if cat <<EOF | gfsview-batch2D sim-4 vectors.gfv ; then :
Save vectors.eps { format = EPS line_width = 1 }
EOF
else
    exit 1
fi

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20
    set output 'laplace.eps'
    set xlabel 'Tau'
    set ylabel 'Ca'
    set logscale y
    plot 'La=120' w l, 'La=1200' w l, 'La=12000' w l
    set output 'convergence.eps'
    set xlabel 'D'
    set ylabel 'Ca'
    unset key
    unset logscale
    plot [][0:]'convergence' w lp ps 3
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
