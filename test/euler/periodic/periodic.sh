if ! $donotrun; then
    rm -f r0 r1 r2

    for r in 0 1 2; do
	for level in 5 6 7; do
	    if sed "s/LEVEL/$level/g" < $1 | \
		sed "s/BOX/$r/g" | \
		gerris2D - | \
		awk -v level=$level '{
                  print level " " $7 " " $9
                }' >> r$r; then :
	    else
		exit 1
	    fi
	done
    done
fi

if cat <<EOF | python > minion1.tex; then :
from check import *
from sys import *
from math import *

print r"""\begin{tabular}{|c|c|c|c|c|c|}\hline
        & \multicolumn{5}{c|}{\$L_2\$} \\\ \hline
        & \$L=5\$   & \$O_2\$ & \$L=6\$    & \$O_2\$ & \$L=7\$  \\\ \hline"""

def order(r,rr):
   for i in range(0,len(r.l)-1):
     y0,y1 = r.l[i][1],r.l[i+1][1]
     yr0,yr1 = rr.l[i][1],rr.l[i+1][1]
     print '& %.2e / {\color{blue}%.2e} & %4.2f / {\color{blue}%4.2f}' % \
       (y0, yr0, log(y0/y1)/log(2.), log(yr0/yr1)/log(2.)),
   print '& %.2e / {\color{blue}%.2e}' % (r.l[i+1][1],rr.l[i+1][1]), r'\\\ \hline'

print 'Uniform',
order(Curve('r0',1,2),Curve('r0.ref',1,2))
print '\$r=1\$',
order(Curve('r1',1,2),Curve('r1.ref',1,2))
print '\$r=2\$',
order(Curve('r2',1,2),Curve('r2.ref',1,2))

print r"""\hline
        & \multicolumn{5}{c|}{\$L_\infty\$} \\\ \hline
        & \$L=5\$   & \$O_\infty\$  & \$L=6\$   & \$O_\infty\$  & \$L=7\$ \\\ \hline"""

print 'Uniform',
order(Curve('r0',1,3),Curve('r0.ref',1,3))
print '\$r=1\$',
order(Curve('r1',1,3),Curve('r1.ref',1,3))
print '\$r=2\$',
order(Curve('r2',1,3),Curve('r2.ref',1,3))

print r'\end{tabular}'
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
for r in ['r0','r1','r2']:
  if (Curve(r,1,2) - Curve(r+'.ref',1,2)).max() > 1e-6 or\
     (Curve(r,1,3) - Curve(r+'.ref',1,3)).max() > 1e-6:
    exit(1)
EOF
else
   exit 1
fi
