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
    echo "usage: time.sh FILE.xmgr"
    exit 1
fi

echo "@description \"./time.sh\""
number=0
error=`mktemp /tmp/error.XXXXXX`
sim=`mktemp /tmp/sim.XXXXXX`
maxlevel=`description $1 | awk '{if ($1 == "GfsRefine") print $2;}'`
levelbox=`description $1 | awk -v level=$maxlevel '{if ($1 == "GfsRefineBox") print $2 - level;}'`
#maxlevel=6
for level in `seq 5 1 $maxlevel`; do
    description $1 | awk -v level=$level -v levelbox=$levelbox '{
	    if ($1 == "GfsRefine")
		print "GfsRefine " level;
            else if ($1 == "GfsRefineBox")
		print "GfsRefineBox " level + levelbox;
            else 
		print $0;
        }' > $sim
    if gerris2D $sim > $error.$level; then
	:
    else
	rm -f $sim $error
	exit 1
    fi
done
rm -f $sim

echo "@WITH G0"
echo "@G0 ON"
echo "@G0 TYPE xy"
echo "@XAXIS LABEL \"Time\""
echo "@xaxis label char size 0.6"
echo "@xaxis ticklabel char size 0.6"
echo "@YAXIS LABEL \"Error (L2 norm)\""
echo "@yaxis label char size 0.6"
echo "@yaxis ticklabel char size 0.6"
number=0
for level in `seq 5 1 $maxlevel`; do
    echo "@LEGEND STRING $number \"`awk -v level=$level 'BEGIN{size = exp(level*log(2)); print size "x" size;}'`\""
    echo "@TARGET S$number"
    echo "@TYPE xy" 
    awk '{ if ($1 == "domain:")
             print $2 " " $4}' < $error.$level
    echo "&"
    number=`expr $number + 1`
done

echo "@WITH G1"
echo "@G1 ON"
echo "@G1 TYPE xy"
echo "@XAXIS LABEL \"Time\""
echo "@xaxis label char size 0.6"
echo "@xaxis ticklabel char size 0.6"
echo "@YAXIS LABEL \"Vorticity (L2 norm)\""
echo "@yaxis label char size 0.6"
echo "@yaxis ticklabel char size 0.6"
number=0
for level in `seq 5 1 $maxlevel`; do
    echo "@LEGEND STRING $number \"`awk -v level=$level 'BEGIN{size = exp(level*log(2)); print size "x" size;}'`\""
    echo "@TARGET S$number"
    echo "@TYPE xy" 
    awk '{if ($3 == "t:") time = $4; 
          else if ($1 == "Vorticity") print time " " $5;}' < $error.$level
    echo "&"
    number=`expr $number + 1`
done

echo "@WITH G2"
echo "@G2 ON"
echo "@G2 TYPE xy"
echo "@XAXIS LABEL \"Time\""
echo "@xaxis label char size 0.6"
echo "@xaxis ticklabel char size 0.6"
echo "@YAXIS LABEL \"Divergence (L2 norm)\""
echo "@yaxis label char size 0.6"
echo "@yaxis ticklabel char size 0.6"
number=0
for level in `seq 5 1 $maxlevel`; do
    echo "@LEGEND STRING $number \"`awk -v level=$level 'BEGIN{size = exp(level*log(2)); print size "x" size;}'`\""
    echo "@TARGET S$number"
    echo "@TYPE xy" 
    awk '{if ($3 == "t:") time = $4;
          else if ($1 == "Divergence")
            print time " " $5;}' < $error.$level
    echo "&"
    number=`expr $number + 1`
done

echo "@WITH G3"
echo "@G3 ON"
echo "@G3 TYPE xy"
echo "@XAXIS LABEL \"Time\""
echo "@xaxis label char size 0.6"
echo "@xaxis ticklabel char size 0.6"
echo "@YAXIS LABEL \"Divergence (Lmax norm)\""
echo "@yaxis label char size 0.6"
echo "@yaxis ticklabel char size 0.6"
number=0
for level in `seq 5 1 $maxlevel`; do
    echo "@LEGEND STRING $number \"`awk -v level=$level 'BEGIN{size = exp(level*log(2)); print size "x" size;}'`\""
    echo "@TARGET S$number"
    echo "@TYPE xy" 
    awk '{if ($3 == "t:") time = $4; 
          else if ($1 == "Divergence")
            print time " " $7;}' < $error.$level
    echo "&"
    number=`expr $number + 1`
done

rm -f $error*
