# !/bin/sh

cat <<EOF
Level		angle of maximum C	maximum C
EOF

wave=`mktemp /tmp/wave.XXXXXX`
res=`mktemp /tmp/wave.XXXXXX`
for l in `seq 5 1 7`; do
  cat <<EOF > $wave
1 0 GfsOcean GfsBox GfsGEdge {} {
 Time { end = 37.80501984 dtmax = 0.1 }
 PhysicalParams { g = 5.87060327757e-3 }
 Init {} {
   P = {
     #include <gsl/gsl_sf_bessel.h>
     #link -lgsl -lgslcblas
     double theta = atan2(y,x);
     double r = sqrt (x*x + y*y);
     double D = 8.83906519983e-2;
     double k = 3.;
     double A = 1./2555510.;
     return A*cos (k*theta)*gsl_sf_bessel_Inu (k, r/D);
   }
   U = {
     #include <gsl/gsl_sf_bessel.h>
     #link -lgsl -lgslcblas
     #define Ik(k,r,D) (gsl_sf_bessel_Inu ((k) - 1., (r)/(D))/(D)\
                         - (k)/(r)*gsl_sf_bessel_Inu ((k), (r)/(D)))
     double theta = atan2(y,x);
     double r = sqrt (x*x + y*y);
     double D = 8.83906519983e-2;
     double k = 3.;
     double sigma = 0.4986;
     double A = 1./2555510.;
     double ur = -A*D*D/5.87060327757e-3*sin (k*theta)*(sigma*Ik (k, r, D) - 
       	                                               k/r*gsl_sf_bessel_Inu (k, r/D));
     double vt =  A*D*D/5.87060327757e-3*cos (k*theta)*(Ik (k, r, D) - 
       	                                               k*sigma/r*gsl_sf_bessel_Inu (k, r/D));
     return ur*cos (theta) - vt*sin (theta);
   }
   V = {
     #include <gsl/gsl_sf_bessel.h>
     #link -lgsl -lgslcblas
     #define Ik(k,r,D) (gsl_sf_bessel_Inu ((k) - 1., (r)/(D))/(D)\
                         - (k)/(r)*gsl_sf_bessel_Inu ((k), (r)/(D)))
     double theta = atan2(y,x);
     double r = sqrt (x*x + y*y);
     double D = 8.83906519983e-2;
     double k = 3.;
     double sigma = 0.4986;
     double A = 1./2555510.;
     double ur = -A*D*D/5.87060327757e-3*sin (k*theta)*(sigma*Ik (k, r, D) - 
       	                                               k/r*gsl_sf_bessel_Inu (k, r/D));
     double vt =  A*D*D/5.87060327757e-3*cos (k*theta)*(Ik (k, r, D) - 
       	                                               k*sigma/r*gsl_sf_bessel_Inu (k, r/D));
     return ur*sin (theta) + vt*cos (theta);
   }
 }
 Refine $l
 GtsSurfaceFile basin.gts
 AdvectionParams { scheme = none }
 ApproxProjectionParams { tolerance = 1e-9 }
 SourceCoriolis {} U 1.
EOF
  for a in `awk 'BEGIN{
                        for (x = -30; x <= 0.; x += 1.) print x;
                        for (x = 0.1; x <= 4.9; x += 0.1) print x;
                        for (x = 5.; x <= 30.; x += 1.) print x;
                      }'`; do
cat <<EOF >> $wave
 EventScript { start = end } { echo -n "$a " }
 OutputCorrelation { start = end } stdout { v = P } {
   s = {
     #include <gsl/gsl_sf_bessel.h>
     #link -lgsl -lgslcblas
     double theta = atan2(y,x) + $a*M_PI/180.;
     double r = sqrt (x*x + y*y);
     double D = 8.83906519983e-2;
     double k = 3.;
     double A = 1./2555510.;
     return A*cos (k*theta)*gsl_sf_bessel_Inu (k, r/D);
   }
   unbiased = 1
 }
EOF
  done >> $wave
cat <<EOF >> $wave
  OutputProgress { istep = 1 } stderr
}
GfsBox { front = Boundary }
EOF
  gerris2D3 $wave | awk '{ print $1 " " $5}' > $res
  if awk -v l=$l 'BEGIN { min1 = 0. } {
    if ($2 > min1) {
      theta = $1;
      min1 = $2;
    }
  } END {
    printf ("%d\t\t%g\t\t\t%g\n", l, theta, min1);
    if (l == 5 && (theta > 4. || min1 < 0.953))
      exit 0;
    else if (l == 6 && (theta > 1. || min1 < 0.995))
      exit 0;
    else if (l == 7 && (theta > 0.6 || min1 < 0.998))
      exit 0;
    exit 1;
  }' < $res; then
    rm -f $wave $res
    exit 1;
  fi
done
rm -f $wave $res
