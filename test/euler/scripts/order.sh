#! /bin/sh

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
    echo "usage: order.sh FILE.xmgr"
    exit 1
fi

error=`mktemp /tmp/error.XXXXXX`
order=`mktemp /tmp/order.XXXXXX`
sim=`mktemp /tmp/sim.XXXXXX`
maxlevel=`description $1 | awk '{if ($1 == "GfsRefine") print $2;}'`
levelbox=`description $1 | awk -v level=$maxlevel '{if ($1 == "GfsRefineBox") print $2 - level;}'`
#maxlevel=6
for level in `seq 3 1 $maxlevel`; do
    description $1 | awk -v level=$level -v levelbox=$levelbox '{
	    if ($1 == "GfsRefine")
		print "GfsRefine " level;
            else if ($1 == "GfsRefineBox")
		print "GfsRefineBox " level + levelbox;
            else 
		print $0;
        }' > $sim
    if res=`gerris2D $sim | awk '{if ($1 == "domain:") print $3 " " $4 " " $5;}'`; then
	:
    else
	rm -f $sim $error $order
	exit 1
    fi
    echo $level $res >> $error
done
rm -f $sim

awk 'BEGIN {n = 0}
{
    level[n] = $1;
    e1[n] = $2; e2[n] = $3; ei[n++] = $4;
}
END {
    for (i = 1; i < n; i++)
	print level[i] " " log(e1[i-1]/e1[i])/log(2) " " log(e2[i-1]/e2[i])/log(2) " " log(ei[i-1]/ei[i])/log(2);
}' < $error > $order

awk '{if ($1 == "@description") print $0}' < $1

echo "@WITH G0"
echo "@G0 ON" 
echo "@G0 TYPE logy" 
echo "@XAXIS LABEL \"Level\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@YAXIS LABEL \"Error\""
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@LEGEND STRING 0 \"first\""
echo "@LEGEND STRING 1 \"second\""
echo "@LEGEND STRING 2 \"infinite\""
echo "@TARGET S0" 
echo "@TYPE xy" 
awk '{print $1 " " $2}' < $error
echo "&" 
echo "@TARGET S1" 
echo "@TYPE xy" 
awk '{print $1 " " $3}' < $error
echo "&" 
echo "@TARGET S2" 
echo "@TYPE xy" 
awk '{print $1 " " $4}' < $error
echo "&" 

rm -f $error

echo "@WITH G1" 
echo "@G1 ON" 
echo "@G1 TYPE xy" 
echo "@XAXIS LABEL \"Level\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@YAXIS LABEL \"Order\""
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@LEGEND STRING 0 \"first\""
echo "@LEGEND STRING 1 \"second\""
echo "@LEGEND STRING 2 \"infinite\""
echo "@TARGET S0" 
echo "@TYPE xy" 
awk '{print $1 " " $2}' < $order
echo "&" 
echo "@TARGET S1" 
echo "@TYPE xy" 
awk '{print $1 " " $3}' < $order
echo "&" 
echo "@TARGET S2" 
echo "@TYPE xy" 
awk '{print $1 " " $4}' < $order
echo "&" 

rm -f $order
