# !/bin/sh

delaunay=$HOME/adds/src/gts-mainline/examples/delaunay

param=`mktemp /tmp/energy.XXXXXX`
geom=`mktemp /tmp/energy.XXXXXX`
energy=`mktemp /tmp/energy.XXXXXX`

points=`mktemp /tmp/energy.XXXXXX`
points1=`mktemp /tmp/energy.XXXXXX`
awk 'BEGIN {
  for (x = -0.51; x <= 0.51; x += 0.01)
    for (y = -0.51; y <= 0.51; y += 0.01)
       print x " " y " " 1.*exp(-100.*(x*x)) - 0.5001;
}' > $points
echo `wc -l $points | awk '{print $1}'`" 0 0" > $points1
cat $points >> $points1
$delaunay < $points1 | transform --revert > $geom
rm -f $points $points1

cat <<EOF > $param
1 0 GfsOcean GfsBox GfsGEdge {} {
  Time { end = 10 dtmax = 1e-2 }
  Refine 5
  GtsSurfaceFile $geom
  ApproxProjectionParams { tolerance = 1e-9 }
  AdvectionParams { scheme = none }
  Init {} { P = { return 1e-3*cos (2.*M_PI*x); } }
  OutputProgress { istep = 1 } stderr
  OutputEnergy { istep = 1 } $energy
}
GfsBox {
  front = Boundary
}
EOF

gerris2D3 $param

if awk '{
  if ($5 + $7 > 5e-7)
    exit 1;
}' < $energy; then
  status=0;
else
  status=1;
fi

rm -f $param $geom $energy

exit $status
