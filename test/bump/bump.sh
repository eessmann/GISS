if test x$donotrun != xtrue; then
    awk 'BEGIN{for (x = -0.5 ; x <= 0.5; x += 0.01) print x,0,0;}' >  xprofile
   gerris2D  bump.gfs
fi

if cat <<EOF | gnuplot; then :
set term postscript eps color enhanced lw 2 18 solid
set output 'profile.eps'

#
# time evolution charge density at origin
# 

set xlabel 't'
set ylabel '{/Symbol r}_{e}'
set key bottom right
a=0.05
perm=2.
K=1.
rhoinic=1./a/sqrt(2*pi)
rho(x)=(rhoinic*exp(-x*K/perm))
plot \
          'timevol' u 1:3 t '', \
          rho(x) t '{/Symbol r}_{e}(0,t)'
EOF
else
    exit 1
fi

if cat <<EOF | gnuplot; then :
set term postscript eps color enhanced lw 2 18 solid
set output 'figure.eps'

#
# time evolution of bump times 0, 2, 4 and 6
# 

set xlabel 'x'
set ylabel '{/Symbol r}_{e}'
set key top right
a=0.05
perm=2.
K=1.
time0=0
time2=2
time4=4
time6=6
rhoinic=1./a/sqrt(2*pi)
rho0(x)=(rhoinic*exp(-0.5*x*x/a/a)*exp(-time0*K/perm))
rho2(x)=(rhoinic*exp(-0.5*x*x/a/a)*exp(-time2*K/perm))
rho4(x)=(rhoinic*exp(-0.5*x*x/a/a)*exp(-time4*K/perm))
rho6(x)=(rhoinic*exp(-0.5*x*x/a/a)*exp(-time6*K/perm))
plot \
          't_0' u 1:2 t '', \
          't_2' u 1:2 t '', \
          't_4' u 1:2 t '', \
          't_6' u 1:2 t '', \
          rho0(x) t '{/Symbol r}_{e}(x,0)',\
          rho2(x) t '{/Symbol r}_{e}(x,2)',\
          rho4(x) t '{/Symbol r}_{e}(x,4)',\
          rho6(x) t '{/Symbol r}_{e}(x,6)'
EOF
else
    exit 1
fi
