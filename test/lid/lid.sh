if ! $donotrun; then
    if gerris2D $1; then :
    else
	exit 1
    fi
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('xprof',3,6) - Curve('xprof.ghia',1,2)).normi() > 2e-2 or \
   (Curve('yprof',2,7) - Curve('yprof.ghia',1,2)).normi() > 1.7e-2:
    exit(1)
EOF
else
   exit 1
fi
