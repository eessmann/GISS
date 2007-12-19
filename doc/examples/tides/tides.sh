# GMT path
PATH=$PATH:/usr/lib/gmt/bin

# The Gerris ocean model uses the following non-dimensional variables
#
# Dimensional variables:
# gravity:
# g* = 9.81 m/s^2
# reference depth:
# H* = 5000 m
# reference length:
# L* = 500 km
# Coriolis f-plane:
# f* = 1e-4 s^-1
# Velocity:
# u* = ... m/s
# Time:
# t* = ... s
# Surface elevation:
# eta* = ... m
# Position:
# x* = ... m
#
# Corresponding non-dimensional variables
# g = g*H*/L*^2f*^2 = 19.62
# u = u*/L*f*
# t = t*f*
# x = x*/L*
# eta = eta*/H*
# pressure:
# p = eta g
#
# -- create bathymetry --
# The domain is centered on lon,lat (deg): 174,-40.8
# default domain size is 500 km
# --angle=40 means that the resulting projection will also be rotated by 40 degrees
# --rel is the relative error allowed when simplifying bathymetry i.e. 2.5% of the total depth
# bathymetry: contains three columns lon (deg), lat (deg), depth (m)
#
# bat2gts also creates two scripts: xy2lolat and lolat2xy which
# convert coordinates in input between a lon,lat coordinates system
# to/from the x,y, coordinate system in which Gerris works.
bat2gts --lat=-40.8 --long=174 --angle=40 --rel=0.025 < bathymetry | transform --tz -0.4994 > bath.gts

# M2 tidal coefficients
#
# coefficients contains the amplitude and phase of the M2
# tide as a function of space i.e. four columns:
# lon (deg) lat (deg) amplitude (m) phase (degree)

# just counts the number of lines in the input file
lines=`wc -l coefficients | awk '{print $1}'`

# transform lon,lat into x,y coordinates
# Gerris defines the tidal M2 mode as:
# AM2*cos (omega*t) + BM2*sin (omega*t)
#
# We first compute A_amp from the amplitude and phase.
# We also need to rescale the tidal amplitude according to:
# eta = eta*g/H*
./lolat2xy < coefficients | awk -v lines=$lines '
BEGIN {
  print lines " 0 0"
} {
  H = 5000
  g = 19.62
  amp = $3
  phi = $4*3.14159265357/180.
  print $1 " " $2 " " g*amp/H*cos(phi)
}' | delaunay > AM2.gts

# Now compute BM2
./lolat2xy < coefficients | awk -v lines=$lines '
BEGIN {
  print lines " 0 0"
} {
  H = 5000
  g = 19.62
  amp = $3
  phi = $4*3.14159265357/180.
  print $1 " " $2 " " g*amp/H*sin(phi)
}' | delaunay > BM2.gts

# Run the simulation
gerris2D3 -m tides.gfs | gfsview2D3 tides.gfv

# Use batch mode of gfsview to generate figures
echo "Save amplitude.eps { format = EPS }" | gfsview-batch2D3 end.gfs amplitude.gfv
echo "Save phase.eps { format = EPS }" | gfsview-batch2D3 end.gfs phase.gfv
echo "Save ellipses.eps { format = EPS }" | gfsview-batch2D3 end.gfs ellipses.gfv
echo "Save residual.eps { format = EPS }" | gfsview-batch2D3 end.gfs residual.gfv
