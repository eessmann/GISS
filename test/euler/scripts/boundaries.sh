#! /bin/sh

if test -e boundaries; then
    :
else
    mkdir boundaries
fi

PATH=$PATH:../../poisson:../../../tools:..
cd boundaries
/bin/sh -c "shapes almgren > boundaries.gts"

cat <<EOF | gerris2D -
1 0 GfsSimulation GfsBox GfsGEdge {} {
  Time { end = 0 }
  ApproxProjectionParams {
    tolerance = 1e-6
  }
  Refine 7
  GtsSurfaceFile boundaries.gts
  Init {} { U = 1 }
  OutputSimulation { start = end } boundaries_128.gfs {
    variables = U,V,P
  }
}
GfsBox { left = BoundaryInflowConstant 1 right = BoundaryOutflow }
EOF

cat <<EOF | gerris2D -
1 0 GfsSimulation GfsBox GfsGEdge {} {
  Time { end = 0 }
  ApproxProjectionParams {
    tolerance = 1e-6
  }
  Refine 8
  GtsSurfaceFile boundaries.gts
  Init {} { U = 1 }
  OutputSimulation { start = end } boundaries_256.gfs {
    variables = U,V,P
  }
}
GfsBox { left = BoundaryInflowConstant 1 right = BoundaryOutflow }
EOF

cat <<EOF | gerris2D -
1 0 GfsSimulation GfsBox GfsGEdge {} {
  Time { end = 0 }
  ApproxProjectionParams {
    tolerance = 1e-6
  }
  Refine 9
  GtsSurfaceFile boundaries.gts
  Init {} { U = 1 }
  OutputSimulation { start = end } boundaries_512.gfs {
    variables = U,V,P
  }
}
GfsBox { left = BoundaryInflowConstant 1 right = BoundaryOutflow }
EOF

error256=`gfscompare2D -n -v boundaries_256.gfs boundaries_512.gfs U 2>&1 | awk '{if ($1 == "total") print $0;}'`
error128=`gfscompare2D -n -v boundaries_128.gfs boundaries_256.gfs U 2>&1 | awk '{if ($1 == "total") print $0;}'`
error256full=`gfscompare2D -f 7 -n -v boundaries_256.gfs boundaries_512.gfs U 2>&1 | awk '{if ($1 == "total") print $0;}'`
error128full=`gfscompare2D -f 7 -n -v boundaries_128.gfs boundaries_256.gfs U 2>&1 | awk '{if ($1 == "total") print $0;}'`

cat <<EOF
% command: boundaries.sh
% legend: Projection test of Almgren et al., 1997.

\begin{table}[htbp]
\begin{center}
\begin{tabular}{||l|c|c|c||c|c|c||} \hline
           & \multicolumn{3}{c||}{All cells} & \multicolumn{3}{c||}{Full 128 cells} \\\\ \hline
           & 128-256  & Rate & 256-512  & 128-256  & Rate & 256-512  \\\\ \hline
EOF

L1_128=`echo $error128 | awk '{ printf ("%.2e", $4)}'`
L1_256=`echo $error256 | awk '{ printf ("%.2e", $4)}'`
L1_128full=`echo $error128full | awk '{ printf ("%.2e", $4)}'`
L1_256full=`echo $error256full | awk '{ printf ("%.2e", $4)}'`
rate=`awk -v e1=$L1_128 -v e2=$L1_256 'BEGIN{printf ("%.2f", log(e1/e2)/log(2))}'`
ratefull=`awk -v e1=$L1_128full -v e2=$L1_256full 'BEGIN{printf ("%.2f", log(e1/e2)/log(2))}'`
echo "\$L_1\$      & $L1_128 & $rate & $L1_256 & $L1_128full & $ratefull & $L1_256full \\\\"

L1_128=`echo $error128 | awk '{ printf ("%.2e", $6)}'`
L1_256=`echo $error256 | awk '{ printf ("%.2e", $6)}'`
L1_128full=`echo $error128full | awk '{ printf ("%.2e", $6)}'`
L1_256full=`echo $error256full | awk '{ printf ("%.2e", $6)}'`
rate=`awk -v e1=$L1_128 -v e2=$L1_256 'BEGIN{printf ("%.2f", log(e1/e2)/log(2))}'`
ratefull=`awk -v e1=$L1_128full -v e2=$L1_256full 'BEGIN{printf ("%.2f", log(e1/e2)/log(2))}'`
echo "\$L_2\$      & $L1_128 & $rate & $L1_256 & $L1_128full & $ratefull & $L1_256full \\\\"

L1_128=`echo $error128 | awk '{ printf ("%.2e", $8)}'`
L1_256=`echo $error256 | awk '{ printf ("%.2e", $8)}'`
L1_128full=`echo $error128full | awk '{ printf ("%.2e", $8)}'`
L1_256full=`echo $error256full | awk '{ printf ("%.2e", $8)}'`
rate=`awk -v e1=$L1_128 -v e2=$L1_256 'BEGIN{printf ("%.2f", log(e1/e2)/log(2))}'`
ratefull=`awk -v e1=$L1_128full -v e2=$L1_256full 'BEGIN{printf ("%.2f", log(e1/e2)/log(2))}'`
echo "\$L_\infty\$ & $L1_128 & $rate & $L1_256 & $L1_128full & $ratefull & $L1_256full \\\\ \hline"

cat <<EOF
\end{tabular}
\end{center}
\caption{Errors and convergence rates for the \$x\$-component of the velocity.}
\end{table}
EOF

error256=`gfscompare2D -n -v boundaries_256.gfs boundaries_512.gfs V 2>&1 | awk '{if ($1 == "total") print $0;}'`
error128=`gfscompare2D -n -v boundaries_128.gfs boundaries_256.gfs V 2>&1 | awk '{if ($1 == "total") print $0;}'`
error256full=`gfscompare2D -f 7 -n -v boundaries_256.gfs boundaries_512.gfs V 2>&1 | awk '{if ($1 == "total") print $0;}'`
error128full=`gfscompare2D -f 7 -n -v boundaries_128.gfs boundaries_256.gfs V 2>&1 | awk '{if ($1 == "total") print $0;}'`

cat <<EOF

\begin{table}[htbp]
\begin{center}
\begin{tabular}{||l|c|c|c||c|c|c||} \hline
           & \multicolumn{3}{c||}{All cells} & \multicolumn{3}{c||}{Full 128 cells} \\\\ \hline
           & 128-256  & Rate & 256-512  & 128-256  & Rate & 256-512  \\\\ \hline
EOF

L1_128=`echo $error128 | awk '{ printf ("%.2e", $4)}'`
L1_256=`echo $error256 | awk '{ printf ("%.2e", $4)}'`
L1_128full=`echo $error128full | awk '{ printf ("%.2e", $4)}'`
L1_256full=`echo $error256full | awk '{ printf ("%.2e", $4)}'`
rate=`awk -v e1=$L1_128 -v e2=$L1_256 'BEGIN{printf ("%.2f", log(e1/e2)/log(2))}'`
ratefull=`awk -v e1=$L1_128full -v e2=$L1_256full 'BEGIN{printf ("%.2f", log(e1/e2)/log(2))}'`
echo "\$L_1\$      & $L1_128 & $rate & $L1_256 & $L1_128full & $ratefull & $L1_256full \\\\"

L1_128=`echo $error128 | awk '{ printf ("%.2e", $6)}'`
L1_256=`echo $error256 | awk '{ printf ("%.2e", $6)}'`
L1_128full=`echo $error128full | awk '{ printf ("%.2e", $6)}'`
L1_256full=`echo $error256full | awk '{ printf ("%.2e", $6)}'`
rate=`awk -v e1=$L1_128 -v e2=$L1_256 'BEGIN{printf ("%.2f", log(e1/e2)/log(2))}'`
ratefull=`awk -v e1=$L1_128full -v e2=$L1_256full 'BEGIN{printf ("%.2f", log(e1/e2)/log(2))}'`
echo "\$L_2\$      & $L1_128 & $rate & $L1_256 & $L1_128full & $ratefull & $L1_256full \\\\"

L1_128=`echo $error128 | awk '{ printf ("%.2e", $8)}'`
L1_256=`echo $error256 | awk '{ printf ("%.2e", $8)}'`
L1_128full=`echo $error128full | awk '{ printf ("%.2e", $8)}'`
L1_256full=`echo $error256full | awk '{ printf ("%.2e", $8)}'`
rate=`awk -v e1=$L1_128 -v e2=$L1_256 'BEGIN{printf ("%.2f", log(e1/e2)/log(2))}'`
ratefull=`awk -v e1=$L1_128full -v e2=$L1_256full 'BEGIN{printf ("%.2f", log(e1/e2)/log(2))}'`
echo "\$L_\infty\$ & $L1_128 & $rate & $L1_256 & $L1_128full & $ratefull & $L1_256full \\\\ \hline"

cat <<EOF
\end{tabular}
\end{center}
\caption{Errors and convergence rates for the \$y\$-component of the velocity.}
\end{table}
EOF
