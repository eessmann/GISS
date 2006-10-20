if ! $donotrun; then
    shapes -n 1000 ellipse | transform --scale 1.6 | transform --tx -0.5 --ty 0.5 > circle.gts
    for La in 120 1200 12000; do
	mu=`echo $La | awk '{print sqrt (0.8/$1)}'`
	tmax=`echo $mu | awk '{print 0.8*0.8/$1}'`
	if sed "s/LEVEL/5/g" < $1 |\
           sed "s/MU/$mu/g" |\
	   sed "s/LAPLACE/$La/g" |\
	   sed "s/DT/0/" |\
           sed "s/end = TMAX/iend = 1/g" | gerris2D - |\
           sed "s/iend = 1/end = $tmax/" | gerris2D - > /dev/null; then :
	else
	    exit 1
	fi
    done

    rm -f convergence
    La=12000
    for level in 4 6 7; do
        mu=`echo $La | awk '{print sqrt (0.8/$1)}'`
	tmax=`echo $mu | awk '{print 0.8*0.8/$1}'`
	if sed "s/LEVEL/$level/g" < $1 |\
           sed "s/MU/$mu/g" |\
	   sed "s/LAPLACE/$La/g" |\
	   sed "s/DT/1e-9/" |\
           sed "s/end = TMAX/iend = 1/g" | gerris2D - |\
           sed "s/iend = 1/end = $tmax/" | gerris2D - > sim-$level; then : 
	else
	    exit 1
	fi
    done
    for level in 4 5 6 7; do
	if awk -v level=$level < E-$La-$level '{
             max2 = $3
             maxi = $4
           }END{print 0.8*2**level " " max2 " " maxi}' >> convergence; then : 
	else
	    exit 1
	fi
    done
fi

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20
    set output 'laplace.eps'
    set xlabel 'Tau'
    set ylabel 'U(D/sigma)^1/2'
    set logscale y
    plot 'La-120-5' w l t "La=120", 'La-1200-5' w l t "La=1200", 'La-12000-5' w l t "La=12000"
    set output 'convergence.eps'
    set xlabel 'D'
    set ylabel 'Shape error'
    set logscale x
    plot [10:120]'convergence' u 1:2 w lp t "RMS" ps 3, 'convergence' u 1:3 w lp t "Max" ps 3, 0.2/(x*x) t "Second order"
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('convergence',1,3) - Curve('convergence.ref',1,3)).max() > 1e-6:
    exit(1)
EOF
else
   exit 1
fi
