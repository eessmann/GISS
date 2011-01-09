if test x$donotrun != xtrue; then
    awk '{printf("%4.3f %4.3f 0.\n",0.,$1)}' < prof.ref > profile
    if ( gerris2D -DLEVEL=7 -DSOLVER=$solver $1 ) ; then :
    else
	exit 1
    fi     	
    
else
    exit 1
fi


if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('prof.dat',3,5) - Curve('prof.ref',1,2)).max() > 2.5e-6:
   print (Curve('prof.dat',3,5) - Curve('prof.ref',1,2)).max()
   exit(1)
EOF
else
   exit 1
fi
