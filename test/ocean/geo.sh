# !/bin/sh

param=`mktemp /tmp/geo.XXXXXX`
error=`mktemp /tmp/geo.XXXXXX`

cat <<EOF > $param
1 0 GfsOcean GfsBox GfsGEdge { variables = PS,Div } {
 # dt = 1000 s
 Time { iend = 1580 dtmax = 0.10285 }
 Refine 6
 # Lx = Ly = 1000 km
 # H0 = 1000 m
 # g = 0.01 m/s^2
 PhysicalParams { g = 9.4534734306584e-4 }
 ApproxProjectionParams { tolerance = 1e-6 }
 AdvectionParams { scheme = none }
 Init {} {
   # e-folding radius = 100 km
   # umax = 1 m/s = sqrt(200)*exp(-1/2)
   U = { return   5.667583815e-4*200.*y*exp (-100.*(x*x + y*y)); }
   V = { return - 5.667583815e-4*200.*x*exp (-100.*(x*x + y*y)); }
   P = { return   5.667583815e-4*exp (-100.*(x*x + y*y)); }
 }
 # f0 = 1.0285e-4 s-1
 SourceCoriolis {} U 1
 OutputErrorNorm { istep = 10 } { awk '{print \$3/1.0285e-4/3600./24. " " \$9*1000e6*1.0285e-4*1.0285e-4/0.01;}' > $error } { v = P } {
   s = { return 5.667583815e-4*exp (-100.*(x*x + y*y)); }
   unbiased = 1
   v = E
 }
 OutputProgress { istep = 10 } stdout
}
GfsBox {
  front = Boundary
}
EOF

gerris2D3 $param

if awk 'BEGIN{emax=0.}{if ($2 > emax) emax=$2;}
END{
  print "maximum geostrophic error: " emax;
  if (emax > 0.015)
    exit 0;
  else
    exit 1;
}' < $error; then
  rm -f $param $error
  exit 1;
fi

cat <<EOF > $param
1 0 GfsOcean GfsBox GfsGEdge { variables = PS,Div } {
 # dt = 1000 s
 Time { iend = 1580 dtmax = 0.10285 }
 Refine 6
 # Lx = Ly = 1000 km
 # H0 = 1000 m
 # g = 0.01 m/s^2
 PhysicalParams { g = 9.4534734306584e-4 }
 ApproxProjectionParams { tolerance = 1e-6 }
 Init {} {
   # e-folding radius = 100 km
   # umax = 1 m/s = sqrt(200)*exp(-1/2)
   U = { return   5.667583815e-4*200.*y*exp (-100.*(x*x + y*y)); }
   V = { return - 5.667583815e-4*200.*x*exp (-100.*(x*x + y*y)); }
   P = { return   5.667583815e-4*exp (-100.*(x*x + y*y)); }
 }
 # f0 = 1.0285e-4 s-1
 # beta = 1.607e-11 m-1s-1
 SourceCoriolis {} U { return 1. + 0.156246961595*(y + 0.5); }
 OutputEnergy { istep = 10 } { awk '{print \$3/1.0285e-4/3600./24. " " (\$5+\$7)/(1.009e-06+5.002e-06)}' > $error }
 OutputProgress { istep = 10 } stdout
}
GfsBox {
  front = Boundary
}
EOF

gerris2D3 $param

if awk 'BEGIN{emin=1.}{if ($2 < emin) emin=$2;}
END{
  print "minimum geostrophic energy: " emin;
  if (emin < 0.96)
    exit 0;
  else
    exit 1;
}' < $error; then
  rm -f $param $error
  exit 1;
fi

rm -f $param $error
