if test x$donotrun != xtrue; then
    if gerris2D bump.gfs; then :
    else
	exit 1
    fi
fi

if gnuplot <<EOF ; then :
set term postscript eps color enhanced lw 2 20 solid

#
# time evolution of charge density at origin
# 

set output 'profile.eps'
set xlabel 't'
set ylabel '{/Symbol r}_{e}'
a=0.05
perm=2.
K=1.
rhoinic=1./a/sqrt(2*pi)
rho(x)=(rhoinic*exp(-x*K/perm))
plot \
          'timevol' u 1:3 t '', \
          rho(x) t '{/Symbol r}_{e}(0,t)'

#
# time evolution of bump times 0, 2, 4 and 6
# 

set output 'figure.eps'
set xlabel 'r'
set ylabel '{/Symbol r}_{e}'
time0=0
time2=2
time4=4
time6=6
rhoinic=1./a/sqrt(2*pi)
rho0(x)=(rhoinic*exp(-0.5*x*x/a/a)*exp(-time0*K/perm))
rho2(x)=(rhoinic*exp(-0.5*x*x/a/a)*exp(-time2*K/perm))
rho4(x)=(rhoinic*exp(-0.5*x*x/a/a)*exp(-time4*K/perm))
rho6(x)=(rhoinic*exp(-0.5*x*x/a/a)*exp(-time6*K/perm))
plot [0:0.2]\
          't_0' u 1:2 t '' lt 1, \
          't_2' u 1:2 t '' lt 2, \
          't_4' u 1:2 t '' lt 3, \
          't_6' u 1:2 t '' lt 4, \
          rho0(x) t '{/Symbol r}_{e}(x,0)' lt 1,\
          rho2(x) t '{/Symbol r}_{e}(x,2)' lt 2,\
          rho4(x) t '{/Symbol r}_{e}(x,4)' lt 3,\
          rho6(x) t '{/Symbol r}_{e}(x,6)' lt 4

#
# Evolution of error norms
# 

set output 'error.eps'
set xlabel 't'
set ylabel 'Error'
plot 'norms.ref' u 3:9 w l t 'Max (ref)', 'norms' u 3:9 w p t 'Max', \
     'norms.ref' u 3:7 w l t 'L2 (ref)',  'norms' u 3:7 w p t 'L2'
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
if (Curve('norms.ref',3,7) - Curve('norms',3,7)).max() > 1e-5 or \
   (Curve('norms.ref',3,9) - Curve('norms',3,9)).max() > 1e-5:
    print (Curve('norms.ref',3,7) - Curve('norms',3,7)).max()
    print (Curve('norms.ref',3,9) - Curve('norms',3,9)).max()
    exit(1)
EOF
else
   exit 1
fi
