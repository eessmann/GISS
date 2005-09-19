if ! $donotrun; then
    for level in 3 4 5 6 7 8 9; do
	if sed "s/LEVEL/$level/g" < $1 | \
	    gerris2D - > res-$level; then :
	else
	    exit 1
	fi
    done
fi

rm -f error
for level in 3 4 5 6 7 8; do
    next=`expr $level + 1`
    echo -n "$level " >> error
    if gfscompare2D -C -c -v sim-$level sim-$next P 2>&1 | \
	awk '{if ($1 == "total") print $4 " " $6 " " $8;}' >> error; then :
    else
	exit 1
    fi
done

if echo "Save solution.eps { format = EPS line_width = 0.25}" | gfsview-batch2D sim-9 solution.gfv; then :
else
    exit 1
fi

if awk '
BEGIN { n = 0 }
{
  l[n] = $1; n1[n] = $2; n2[n] = $3; ni[n++] = $4;
}
END {
  for (i = 1; i < n; i++)
    print l[i] " " log(n1[i-1]/n1[i])/log(2.) " " log(n2[i-1]/n2[i])/log(2.) " " log(ni[i-1]/ni[i])/log(2.);
}' < error > order; then :
else
    exit 1
fi

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20
    set output 'residual.eps'
    set xlabel 'CPU time'
    set ylabel 'Maximum residual'
    set logscale y
    plot 'res-7.ref' u 1:3 t 'ref' w lp, 'res-7' u 1:3 t '' w lp
    set output 'rate.eps'
    set xlabel 'V-cycle'
    set ylabel 'Cumulative residual reduction factor'
    unset logscale
    plot 'res-7.ref' u 2:4 t 'ref' w lp, 'res-7' u 2:4 t '' w lp
    set output 'error.eps'
    set xlabel 'Level'
    set ylabel 'Error norms'
    set key
    set logscale y
    plot 'error.ref' u 1:2 t '1 (ref)' w lp, \
         'error.ref' u 1:3 t '2 (ref)' w lp, \
         'error.ref' u 1:4 t 'max (ref)' w lp, \
         'error' u 1:2 t '1' w lp, \
         'error' u 1:3 t '2' w lp, \
         'error' u 1:4 t 'max' w lp
    set output 'order.eps'
    set xlabel 'Level'
    set ylabel 'Order'
    set key
    unset logscale
    set xtics 0,1
    set ytics 0,1
    set grid
    plot [][0:3] 'order.ref' u 1:2 t '1 (ref)' w lp, \
                 'order.ref' u 1:3 t '2 (ref)' w lp, \
                 'order.ref' u 1:4 t 'max (ref)' w lp, \
                 'order' u 1:2 t '1' w lp, \
                 'order' u 1:3 t '2' w lp, \
                 'order' u 1:4 t 'max' w lp
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
c = Curve()
for p in Curve('res-7.ref',1,3).l:
    c.l.append((p[0]+0.1, p[1]))
if (Curve('res-7',1,3) - c).max() > 1e-8 or\
   (Curve('error',1,4) - Curve('error.ref',1,4)).max() > 1e-6:
    exit(1)
EOF
else
   exit 1
fi
