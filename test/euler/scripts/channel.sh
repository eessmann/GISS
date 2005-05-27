#! /bin/sh

if test -e channel; then
    :
else
    mkdir channel
fi

PATH=$PATH:../../poisson:../../../tools:..
cd channel
/bin/sh -c "shapes channel | transform --revert --scale 4 --tx 1.5 > channel.gts"

cat <<EOF | gerris2D - > log128
4 3 GfsSimulation GfsBox GfsGEdge {} {
  Time { end = 1 }
  Refine 5
  GtsSurfaceFile channel.gts
  Init {} { U = 1 }
  OutputSolidStats {} stdout
  OutputTime { istep = 10 } stdout
  OutputProjectionStats { istep = 10 } stdout
  OutputSimulation { start = end } channel_128.gfs {
    variables = U,V
  }
}
GfsBox { left = BoundaryInflowConstant 1 }
GfsBox {}
GfsBox {}
GfsBox { right = BoundaryOutflow }
1 2 right
2 3 right
3 4 right
EOF

cat <<EOF | gerris2D - > log256
4 3 GfsSimulation GfsBox GfsGEdge {} {
  Time { end = 1 }
  Refine 6
  GtsSurfaceFile channel.gts
  Init {} { U = 1 }
  OutputSolidStats {} stdout
  OutputTime { istep = 10 } stdout
  OutputProjectionStats { istep = 10 } stdout
  OutputSimulation { start = end } channel_256.gfs {
    variables = U,V
  }
}
GfsBox { left = BoundaryInflowConstant 1 }
GfsBox {}
GfsBox {}
GfsBox { right = BoundaryOutflow }
1 2 right
2 3 right
3 4 right
EOF

cat <<EOF | gerris2D - > log512
4 3 GfsSimulation GfsBox GfsGEdge {} {
  Time { end = 1 }
  Refine 7
  GtsSurfaceFile channel.gts
  Init {} { U = 1 }
  OutputSolidStats {} stdout
  OutputTime { istep = 10 } stdout
  OutputProjectionStats { istep = 10 } stdout
  OutputSimulation { start = end } channel_512.gfs {
    variables = U,V
  }
}
GfsBox { left = BoundaryInflowConstant 1 }
GfsBox {}
GfsBox {}
GfsBox { right = BoundaryOutflow }
1 2 right
2 3 right
3 4 right
EOF

error256=`gfscompare2D -n -v channel_256.gfs channel_512.gfs U 2>&1 | awk '{if ($1 == "total") print $0;}'`
error128=`gfscompare2D -n -v channel_128.gfs channel_256.gfs U 2>&1 | awk '{if ($1 == "total") print $0;}'`
error256full=`gfscompare2D -f 5 -n -v channel_256.gfs channel_512.gfs U 2>&1 | awk '{if ($1 == "total") print $0;}'`
error128full=`gfscompare2D -f 5 -n -v channel_128.gfs channel_256.gfs U 2>&1 | awk '{if ($1 == "total") print $0;}'`

cat <<EOF
% command: channel.sh
% legend: Channel test of Almgren et al., 1997.

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

error256=`gfscompare2D -n -v channel_256.gfs channel_512.gfs V 2>&1 | awk '{if ($1 == "total") print $0;}'`
error128=`gfscompare2D -n -v channel_128.gfs channel_256.gfs V 2>&1 | awk '{if ($1 == "total") print $0;}'`
error256full=`gfscompare2D -f 5 -n -v channel_256.gfs channel_512.gfs V 2>&1 | awk '{if ($1 == "total") print $0;}'`
error128full=`gfscompare2D -f 5 -n -v channel_128.gfs channel_256.gfs V 2>&1 | awk '{if ($1 == "total") print $0;}'`

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
