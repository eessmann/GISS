if test x$donotrun != xtrue; then
   gerris2D swirl.gfs
fi

if cat <<EOF | gnuplot ; then :
    set term postscript eps lw 3 color solid 20 enhanced
    set output 'profiles.eps'
    set xlabel '{/Symbol z}'
    set key 6,0.8
    plot [0:6]'analytical' u 1:2 w l t '-F({/Symbol z})', 'nu' u 1:2 t 'Gerris', \
              'analytical' u 1:3 w l t 'G({/Symbol z})', 'nu' u 1:3 t 'Gerris'
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
print (Curve('nu',1,2) - Curve('analytical',1,2)).normi()
print (Curve('nu',1,3) - Curve('analytical',1,3)).normi()
if (Curve('nu',1,2) - Curve('analytical',1,2)).normi() > 0.011024 or \
   (Curve('nu',1,3) - Curve('analytical',1,3)).normi() > 0.009161 :
    exit(1)
EOF
else
   exit 1
fi
