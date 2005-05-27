rm -f reynolds

for level in 5 6 7; do
  sed "s/LEVEL/$level/" < stationary1.gfs | gerris2D - | awk -v m=1 -v level=$level '{
    time = $3
    ke = $5
    if (time == 0)
      ke0 = ke;
  }END{
    a = -log(ke/ke0)/time
    nu = a/(4.*(2.*m*3.14159265359)^2)
    print level " " 1./nu
  }' >> reynolds
done

cat <<EOF | gnuplot
    set term postscript eps
    set output 'divmax.eps'
    set xlabel 'Time'
    set ylabel 'Divergence Max'
    plot [0:2]'div5' u 3:9 t "5" w l, 'div6' u 3:9 t "6" w l, 'div7' u 3:9 t "7" w l
    set output 'divL2.eps'
    set ylabel 'Divergence L2'
    plot [0:2]'div5' u 3:7 t "5" w l, 'div6' u 3:7 t "6" w l, 'div7' u 3:7 t "7" w l
    set output 'kinetic.eps'
    set ylabel 'Kinetic energy'
    plot [0:2]'kinetic5' u 3:5 t "5" w l, 'kinetic6' u 3:5 t "6" w l, 'kinetic7' u 3:5 t "7" w l
    set output 'reynolds.eps'
    set xlabel 'Level'
    set ylabel 'Effective Reynolds number'
    set nokey
    set logscale y
    plot 'reynolds' u 1:2 w lp
EOF

rm -f div? kinetic? reynolds
