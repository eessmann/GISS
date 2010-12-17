levels="7 8 9 10 11 12 13"

if test x$donotrun != xtrue; then
	for i in $levels; do
	    if gerris2D -DLEVEL=$i source.gfs 2> log-$i; then :
	    else
		exit 1
	    fi
	done
fi

rm -rf error
    for i in $levels; do
	tail -n 1 error-$i >> error
    done

    if gnuplot <<EOF ; then :
set term postscript eps color enhanced 
set output 'error.eps'
set logscale
set xlabel 'Spatial resolution'
set ylabel 'Relative error norms'
plot 'error' u (2**\$1):2 t 'L1' w lp ps 2, 'error' u (2**\$1):2 t 'L2' w lp ps 2
EOF
    else
	exit 1
    fi

if echo "Save stdout { width = 800 height = 200 }" | \
    gfsview-batch2D end.gfs source.gfv | \ 
    convert -colors 256 ppm:- velfield.eps; then :
else
     exit 1
fi

if echo "Save stdout { width = 800 height = 200 }" | \
    gfsview-batch2D end.gfs error.gfv | \ 
    convert -colors 256 ppm:- localerror.eps; then :
else
     exit 1
fi

