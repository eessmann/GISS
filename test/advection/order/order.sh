#! /bin/sh

if test -z "$1"; then
    echo "usage: order.sh COMMAND [MAXLEVEL]"
    echo "  MAXLEVEL: maximum level to test (default is 7)"
    exit 1
fi

if test -z "$2"; then
    maxlevel=7
else
    maxlevel=$2
fi

PATH=$PATH:..
command=$1
rm -f /tmp/error
for level in `seq 3 1 $maxlevel`; do
    res=`/bin/sh -c "$command -N -e -l $level 2>&1" | awk '{if ($1 == "Error" && $3 != 0 && $3 % 4 == 0) print $3 " " $5 " " $7 " " $9;}'`
    echo $level $res >> /tmp/error
done

awk 'BEGIN {n = 0} {
    level[n] = $1;
    e11[n] = $3; e21[n] = $4; ei1[n] = $5;
    e12[n] = $7; e22[n] = $8; ei2[n++] = $9;
}
END {
    for (i = 1; i < n; i++)
	print level[i] " " sqrt (e11[i-1]/e11[i]) " " sqrt (e21[i-1]/e21[i]) " " sqrt (ei1[i-1]/ei1[i]) " " sqrt (e12[i-1]/e12[i]) " " sqrt (e22[i-1]/e22[i]) " " sqrt (ei2[i-1]/ei2[i]);
}' < /tmp/error > /tmp/order

echo "@description \"$command\""

echo "@WITH G0" 
echo "@G0 ON" 
echo "@G0 TYPE logy" 
echo "@XAXIS LABEL \"Level\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@YAXIS LABEL \"Error (Pi)\""
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@LEGEND STRING 0 \"first\""
echo "@LEGEND STRING 1 \"second\""
echo "@LEGEND STRING 2 \"infinite\""
echo "@TARGET S0" 
echo "@TYPE xy" 
awk '{print $1 " " $3}' < /tmp/error
echo "&" 
echo "@TARGET S1" 
echo "@TYPE xy" 
awk '{print $1 " " $4}' < /tmp/error
echo "&" 
echo "@TARGET S2" 
echo "@TYPE xy" 
awk '{print $1 " " $5}' < /tmp/error
echo "&" 

echo "@WITH G2"
echo "@G2 ON" 
echo "@G2 TYPE logy" 
echo "@XAXIS LABEL \"Level\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@YAXIS LABEL \"Error (2*Pi)\""
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@LEGEND STRING 0 \"first\""
echo "@LEGEND STRING 1 \"second\""
echo "@LEGEND STRING 2 \"infinite\""
echo "@TARGET S0" 
echo "@TYPE xy" 
awk '{print $1 " " $7}' < /tmp/error
echo "&" 
echo "@TARGET S1" 
echo "@TYPE xy" 
awk '{print $1 " " $8}' < /tmp/error
echo "&" 
echo "@TARGET S2" 
echo "@TYPE xy" 
awk '{print $1 " " $9}' < /tmp/error
echo "&" 

rm -f /tmp/error

echo "@WITH G1" 
echo "@G1 ON" 
echo "@G1 TYPE xy" 
echo "@XAXIS LABEL \"Level\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@YAXIS LABEL \"Order (Pi)\""
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@LEGEND STRING 0 \"first\""
echo "@LEGEND STRING 1 \"second\""
echo "@LEGEND STRING 2 \"infinite\""
echo "@TARGET S0" 
echo "@TYPE xy" 
awk '{print $1 " " $2}' < /tmp/order
echo "&" 
echo "@TARGET S1" 
echo "@TYPE xy" 
awk '{print $1 " " $3}' < /tmp/order
echo "&" 
echo "@TARGET S2" 
echo "@TYPE xy" 
awk '{print $1 " " $4}' < /tmp/order
echo "&" 

echo "@WITH G3"
echo "@G3 ON" 
echo "@G3 TYPE xy" 
echo "@XAXIS LABEL \"Level\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@YAXIS LABEL \"Order (2*Pi)\""
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@LEGEND STRING 0 \"first\""
echo "@LEGEND STRING 1 \"second\""
echo "@LEGEND STRING 2 \"infinite\""
echo "@TARGET S0" 
echo "@TYPE xy" 
awk '{print $1 " " $5}' < /tmp/order
echo "&" 
echo "@TARGET S1" 
echo "@TYPE xy" 
awk '{print $1 " " $6}' < /tmp/order
echo "&" 
echo "@TARGET S2" 
echo "@TYPE xy" 
awk '{print $1 " " $7}' < /tmp/order
echo "&" 

rm -f /tmp/order
