if test x$donotrun != xtrue; then
    awk 'BEGIN{for (i=0; i < 100; i++) {printf("%4.3f %4.3f %4.3f\n",0.,-0.75+i/100,0.) > "xprofile";  }}'
    awk 'BEGIN{for (i=0; i < 100; i++) {printf("%4.3f %4.3f %4.3f\n",-0.5+i/100,-0.25,0.) > "yprofile";  }}'

    if gerris2D $1; then :
    else
	exit 1
    fi
fi

cat <<EOF| gnuplot
       set term pos enhanced eps color lw 3 solid 20
       set out 'xprof.eps'
       set xlabel 'Y'
       set ylabel 'U'
       plot [-0.5:0.5]'xprof.ghia' u 1:2 t 'Ghia et al.' w p ps 2 pt 9,\
             'xprof' w l t 'Gerris stretch'
       
       set out 'yprof.eps'
       set xlabel 'X'
       set ylabel 'V'
       
       plot  [-0.5:0.5]'yprof.ghia' u 1:2 t 'Ghia et al.' w p ps 2 pt 9,\
             'yprof' w l t 'Gerris stretch'

EOF

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('xprof',1,2) - Curve('xprof.ghia',1,2)).normi() > 2e-2 or \
   (Curve('yprof',1,2) - Curve('yprof.ghia',1,2)).normi() > 1.2e-2:
    print (Curve('xprof',1,2) - Curve('xprof.ghia',1,2)).normi()
    print (Curve('yprof',1,2) - Curve('yprof.ghia',1,2)).normi()
    exit(1)
EOF
else
   exit 1
fi