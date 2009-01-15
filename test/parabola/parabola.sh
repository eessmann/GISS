levels="5 6 7 8 9"

if ! $donotrun; then
    for level in $levels; do
	if gerris2D -DLEVEL=$level parabola.gfs; then :
	else
	    exit 1
	fi
    done
fi

for level in $levels; do
    awk -v level=$level 'BEGIN{
      s1 = 0.;
      s2 = 0.;
      smax = 0.;
      n = 0;
      h0 = 10.;
    }{
      n++;
      s1 += $5;
      s2 += $7*$7;
      if ($9 > smax) smax = $9;
    }END { print level, s1/n/h0, sqrt(s2/n)/h0, smax/h0; }' < error-$level
done > error

for level in $levels; do
    if  paste U-$level vol-$level | awk -v level=$level 'BEGIN {
      sum = 0.; 
      n = 0; 

      h0 = 10.
      a = 3000.
      tau = 1e-3
      B = 5.
      G = 9.81
      p = sqrt (8.*G*h0)/a
      s = sqrt (p*p - tau*tau)/2.
    } {
          u0 = $5/$10;
          t = $3;
          ref = B*exp (-tau*t/2.)*sin (s*t);
          sum += (u0 - ref)*(u0 - ref);
          n += 1;
        }
        END {
          print level, sqrt (sum/n);
        }'; then
	:
    else
	exit 1;
    fi
done > error-u

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

if awk '
BEGIN { n = 0 }
{
  l[n] = $1; n2[n++] = $2;
}
END {
  for (i = 1; i < n; i++)
    print l[i] " " log(n2[i-1]/n2[i])/log(2.);
}' < error-u > order-u; then :
else
    exit 1
fi

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20

    set output 'error.eps'
    set xlabel 'Level'
    set ylabel 'Relative error norms'
    set key
    set logscale y
    set xtics 0,1
    set grid
    plot 'error.ref' u 1:2 t '1 (ref)' w lp, \
         'error.ref' u 1:3 t '2 (ref)' w lp, \
         'error.ref' u 1:4 t 'max (ref)' w lp, \
         'error' u 1:2 t '1' w lp, \
         'error' u 1:3 t '2' w lp, \
         'error' u 1:4 t 'max' w lp
   
    set output 'error-u.eps'
    set ylabel 'u0 relative error L2 norm'
    plot 'error-u' u 1:(\$2/4.) t '' w lp
   
    set output 'order.eps'
    set xlabel 'Level'
    set ylabel 'Order'
    unset logscale
    set ytics 0,1
    plot [][-1:2] 'order.ref' u 1:2 t '1 (ref)' w lp, \
                  'order.ref' u 1:3 t '2 (ref)' w lp, \
                  'order.ref' u 1:4 t 'max (ref)' w lp, \
                  'order' u 1:2 t '1' w lp, \
                  'order' u 1:3 t '2' w lp, \
                  'order' u 1:4 t 'max' w lp

    set output 'order-u.eps'
    plot [][-1:2]'order-u' u 1:2 t '' w lp

    set output 'u0.eps'

    h0 = 10.
    a = 3000.
    tau = 1e-3
    B = 5.
    G = 9.81
    p = sqrt (8.*G*h0)/a
    s = sqrt (p*p - tau*tau)/2.
    u0(t) = B*exp (-tau*t/2.)*sin (s*t)

    set xtics auto
    set ytics auto
    unset grid
    set ylabel 'u0'
    set xlabel 'Time'
    plot u0(x) t 'Analytical', '< paste U-7 vol-7' u 3:(\$5/\$10) every 2 w p t 'Gerris'

    set output 'elevation.eps'
    set xlabel 'x (m)'
    set ylabel 'z (m)'
    t = 1500
    psi(x) = a*a*B*B*exp (-tau*t)/(8.*G*G*h0)*(- s*tau*sin (2.*s*t) + \
      (tau*tau/4. - s*s)*cos (2.*s*t)) - B*B*exp(-tau*t)/(4.*G) - \
      exp (-tau*t/2.)/G*(B*s*cos (s*t) + tau*B/2.*sin (s*t))*x + h0
    bed(x) = h0*(x/a)**2
    set key top center
    plot [-5000:5000] \
      'sim-6-1500.txt' u 1:7:8 w filledcu lc 3 t 'Numerical', \
      psi(x) > bed(x) ? psi(x) : bed(x) lc 2 t 'Analytical', \
      bed(x) lw 3 lc 1 lt 1 t 'Bed profile'
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('error',1,4) - Curve('error.ref',1,4)).max() > 1e-4:
    print (Curve('error',1,4) - Curve('error.ref',1,4)).max()
    exit(1)
EOF
else
   exit 1
fi
