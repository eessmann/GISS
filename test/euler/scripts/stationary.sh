#! /bin/sh
# Estimation of numerical dissipation
# see @TechReport{rider95,
#  author = 	 {W. J. Rider},
#  title = 	 {Approximate projection methods for incompressible flows: Implementation, variants and robustness},
#  institution =  {Los Alamos National Laboratory},
#  year = 	 1995,
#  number =	 {LA-UR-2000},
#  local_url =    {rider95.ps.gz},
#  url =          {http://www-xdiv.lanl.gov/XHM/personnel/wjr/Web_papers/proj/proj.ps.Z}
#  pages=81-85
# } for details

description ()
{
    awk 'BEGIN { n = 0 }
         {
           if ($1 == "@description") {
	     if (n > 0) { 
               for (i = 2; i <= NF; i++)
	         printf ("%s ", $i);
	       printf ("\n");
	     }
             n++
           }
         }' < $1 | sed 's/\"//g'
}

if test -z "$1"; then
    echo "usage: stationary FILE.xmgr"
    exit 1
fi

sim=`mktemp /tmp/sim.XXXXXX`
maxlevel=`description $1 | awk '{if ($1 == "GfsRefine") print $2;}'`
levelbox=`description $1 | awk -v level=$maxlevel '{if ($1 == "GfsRefineBox") print $2 - level;}'`
#maxlevel=6
m=`description $1 | awk '{if ($1 == "GfsInitStationary") print $3;}'`
levels=`seq 5 1 $maxlevel`
for level in $levels; do
    description $1 | awk -v level=$level -v levelbox=$levelbox '{
	    if ($1 == "GfsRefine")
		print "GfsRefine " level;
            else if ($1 == "GfsRefineBox")
		print "GfsRefineBox " level + levelbox;
            else 
		print $0;
        }' > $sim
    if gerris2D $sim > /tmp/log$level; then
	:
    else
	rm -f $sim /tmp/log*
	exit 1
    fi
done
rm -f $sim

awk '{if ($1 == "##") print $0}' < $1

echo "@description \"./stationary.sh\""

echo "@WITH G0"
echo "@G0 ON" 
echo "@G0 TYPE xy"
echo "@XAXIS LABEL \"Time\""
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@YAXIS LABEL \"Divergence (L2 norm)\""
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 

n=0
for level in $levels; do
    label=`awk -v level=$level 'BEGIN{size = exp(level*log(2)); 
                                  print size "x" size;}'`
    echo "@LEGEND STRING $n \"$label\""
    echo "@TARGET S$n"
    echo "@TYPE xy" 
    awk '{if ($3 == "t:") 
	    time = $4; 
	  else if ($1 == "Divergence")
	    print time " " $5;
         }' < /tmp/log$level
    echo "&" 
    n=`expr $n + 1`
done

echo "@WITH G1"
echo "@G1 ON" 
echo "@G1 TYPE xy"
echo "@XAXIS LABEL \"Time\""
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@YAXIS LABEL \"Divergence (Lmax norm)\""
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 

n=0
for level in $levels; do
    label=`awk -v level=$level 'BEGIN{size = exp(level*log(2)); 
                                  print size "x" size;}'`
    echo "@LEGEND STRING $n \"$label\""
    echo "@TARGET S$n"
    echo "@TYPE xy" 
    awk '{if ($3 == "t:") 
	    time = $4; 
	  else if ($1 == "Divergence")
	    print time " " $7;
         }' < /tmp/log$level
    echo "&" 
    n=`expr $n + 1`
done

echo "@WITH G2" 
echo "@G2 ON" 
echo "@G2 TYPE xy" 
echo "@XAXIS LABEL \"Time\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@YAXIS LABEL \"Kinetic energy\""
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 

n=0
for level in $levels; do
    label=`awk -v level=$level 'BEGIN{size = exp(level*log(2));
                                  print size "x" size;}'`
    echo "@LEGEND STRING $n \"$label\""
    echo "@TARGET S$n"
    echo "@TYPE xy" 
    awk '{if ($3 == "t:") 
	    time = $4; 
	  else if ($1 == "Velocity2")
	    print time " " $3;
         }' < /tmp/log$level
    echo "&" 
    n=`expr $n + 1`
done

echo "@WITH G3"
echo "@G2 ON" 
echo "@G2 TYPE xy" 
echo "@XAXIS LABEL \"Level\""
echo "@xaxis label char size 0.6"
echo "@xaxis ticklabel char size 0.6"
echo "@YAXIS LABEL \"Effective Reynolds number\""
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@TARGET S0"
echo "@TYPE xy" 

for level in $levels; do
    echo $level `awk -v m=$m '{
           if ($3 == "t:")
	     time = $4;
	   else if ($1 == "Velocity2") {
	     if (time == 0)
	       ke0 = $3;
	     else
	       ke = $3;
	   }
         }END{
	   a = -log(ke/ke0)/time;
	   nu = a/(4.*(2.*m*3.14159265359)^2);
	   print 1./nu;
         }' < /tmp/log$level`
done

echo "&" 

rm -f /tmp/log*
