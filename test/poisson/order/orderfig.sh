#! /bin/sh
# $1: output of poisson
# $2: description (command line used to run order.sh)

if test -z "$1"; then
    echo "usage: orderfig.sh INPUT DESCRIPTION"
    echo "  INPUT: output of poisson"
    echo "  DESCRIPTION: command line used to run poisson"
    exit 1;
fi

echo "@description \"$2\""

echo "@WITH G0" 
echo "@G0 ON" 
echo "@G0 TYPE logy" 
echo "@XAXIS LABEL \"Level\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@YAXIS LABEL \"Total error norm\"" 
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@LEGEND STRING 0 \"first\"" 
echo "@LEGEND STRING 1 \"second\"" 
echo "@LEGEND STRING 2 \"infty\"" 
echo "@TARGET G0.S0"
echo "@TYPE xy"
awk '{print $1 " " $2}' < $1
echo "&" 
echo "@TARGET G0.S1" 
echo "@TYPE xy" 
awk '{print $1 " " $3}' < $1
echo "&" 
echo "@TARGET G0.S2" 
echo "@TYPE xy" 
awk '{print $1 " " $4}' < $1
echo "&" 

echo "@WITH G1" 
echo "@G1 ON" 
echo "@G1 TYPE xy" 
echo "@XAXIS LABEL \"Level\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@YAXIS LABEL \"Total order\"" 
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@LEGEND STRING 0 \"first\"" 
echo "@LEGEND STRING 1 \"second\"" 
echo "@LEGEND STRING 2 \"infty\"" 
echo "@TARGET G1.S0" 
echo "@TYPE xy" 
awk -f order.awk < $1 | awk '{print $1 " " $2}' 
echo "&" 
echo "@TARGET G1.S1" 
echo "@TYPE xy" 
awk -f order.awk < $1 | awk '{print $1 " " $3}' 
echo "&" 
echo "@TARGET G1.S2" 
echo "@TYPE xy" 
awk -f order.awk < $1 | awk '{print $1 " " $4}' 
echo "&" 

echo "@WITH G2" 
echo "@G2 ON" 
echo "@G2 TYPE logy" 
echo "@XAXIS LABEL \"Level\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@YAXIS LABEL \"Refined error norm\"" 
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@LEGEND STRING 0 \"first\"" 
echo "@LEGEND STRING 1 \"second\"" 
echo "@LEGEND STRING 2 \"infty\"" 
echo "@TARGET G2.S0"
echo "@TYPE xy"
awk '{print $1 " " $5}' < $1
echo "&" 
echo "@TARGET G2.S1" 
echo "@TYPE xy" 
awk '{print $1 " " $6}' < $1
echo "&" 
echo "@TARGET G2.S2" 
echo "@TYPE xy" 
awk '{print $1 " " $7}' < $1
echo "&" 

echo "@WITH G3" 
echo "@G3 ON" 
echo "@G3 TYPE xy" 
echo "@XAXIS LABEL \"Level\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@YAXIS LABEL \"Refined order\"" 
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@LEGEND STRING 0 \"first\"" 
echo "@LEGEND STRING 1 \"second\"" 
echo "@LEGEND STRING 2 \"infty\"" 
echo "@TARGET G3.S0" 
echo "@TYPE xy" 
awk -f order.awk < $1 | awk '{print $1 " " $5}' 
echo "&" 
echo "@TARGET G3.S1" 
echo "@TYPE xy" 
awk -f order.awk < $1 | awk '{print $1 " " $6}' 
echo "&" 
echo "@TARGET G3.S2" 
echo "@TYPE xy" 
awk -f order.awk < $1 | awk '{print $1 " " $7}' 
