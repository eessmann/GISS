if test x$donotrun != xtrue; then
    awk 'BEGIN{for (x = 0./100.; x <= 0.5; x += 1./100.) print 0,x,0;}' > yprofile
    for level in 5 6 7; do
	if sed "s/LEVEL/$level/g" < $1 | \
           gerris2D -; then :
	else
	    exit 1
	fi
    done
fi

if cat <<EOF | gnuplot; then :
set term postscript eps color enhanced lw 2 18 solid
set output 'efield.eps'
#
# electric distribution
# 

set xlabel 'r'
set ylabel 'E'
set key top right
perm=1.
R0=0.1
rhoinic=0.5
rho(x)= x<R0 ?  0. : 0.5*R0*R0*rhoinic/perm/x
plot \
          'prof-5' u 1:2 t 'R_{o}/h_{min} = 3.2', \
          'prof-6' u 1:2 t 'R_{o}/h_{min} = 6.4', \
          'prof-7' u 1:2 t 'R_{o}/h_{min} = 12.8', \
          rho(x) t 'steady electric field'
EOF
else
    exit 1
fi

if cat <<EOF | gnuplot; then :
set term postscript eps color enhanced lw 2 18 solid
set output 'charge.eps'
#
# time evolution of global amount of charge
# 

set xlabel 't'
set ylabel 'Q (%  Error)'
set key bottom right

plot \
          'rhoe-5' u 1:3 t 'R_{o}/h_{min} = 3.2', \
          'rhoe-6' u 1:3 t 'R_{o}/h_{min} = 6.4', \
          'rhoe-7' u 1:3 t 'R_{o}/h_{min} = 12.8' 
EOF
else
    exit 1
fi


