#!/bin/sh

bottom=.255
top=.49

tac <<EOF | shapes - > dumbell.gts
-0.51 -0.51
-0.51 $bottom
-0.1 $bottom
-0.1 $top
0.1 $top
0.1 $bottom
0.51 $bottom
0.51 -0.51
EOF

for solver in gerris hypre; do

    awk -v solver=$solver '{
         if ($1 == "PARAMS") {
            if ( solver == "hypre") {
               print "  GModule hypre { tolerance = 1e-30 ncyclemax = 100 verbose = 1}"}
            if ( solver == "gerris") {
               print "  ApproxProjectionParams { nitermax = 1000 minlevel = 1 tolerance = 1e-30 }"}
         }
         else {print $0}}' dumbell.gfs > tmp.gfs

    if gerris2D tmp.gfs | awk '{
           if ($1 == "residual.infty:" && $3 > 6.621e-02)
              exit (1);
             }'; then
               :
    else
	exit 1
    fi
    rm -f tmp.gfs

done