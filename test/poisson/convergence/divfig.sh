#! /bin/sh
# $1: stderr output of poisson
# $2: description (command line used to run poisson)

if test -z "$1"; then
    echo "usage: divfig.sh INPUT DESCRIPTION"
    echo "  INPUT: output of poisson"
    echo "  DESCRIPTION: command line used to run poisson"
    exit 1;
fi

echo "@description \"$2\""

echo "@with g0" 
echo "@g0 on" 
echo "@g0 type logy" 
echo "@xaxis label \"Iteration\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@yaxis label \"Divergence\"" 
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@legend string 0 \"first\""
echo "@legend string 1 \"second\""
echo "@legend string 2 \"infinite\""
echo "@target G0.S0" 
echo "@type xy" 
awk -f divergence.awk < $1 | awk '{print $1 " " $3}' 
echo "&" 
echo "@target G0.S1" 
echo "@type xy" 
awk -f divergence.awk < $1 | awk '{print $1 " " $4}' 
echo "&" 
echo "@target G0.S2" 
echo "@type xy" 
awk -f divergence.awk < $1 | awk '{print $1 " " $5}' 
echo "&" 

echo "@with g1" 
echo "@g1 on" 
echo "@g1 type logy" 
echo "@xaxis label \"Time\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@yaxis label \"Divergence\"" 
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@legend string 0 \"first\"" 
echo "@legend string 1 \"second\"" 
echo "@legend string 2 \"infinite\"" 
echo "@target G1.S0" 
echo "@type xy" 
awk -f divergence.awk < $1 | awk '{print $2 " " $3}' 
echo "&" 
echo "@target G1.S1" 
echo "@type xy" 
awk -f divergence.awk < $1 | awk '{print $2 " " $4}' 
echo "&" 
echo "@target G1.S2" 
echo "@type xy" 
awk -f divergence.awk < $1 | awk '{print $2 " " $5}' 
echo "&" 

echo "@with g2" 
echo "@g2 on" 
echo "@g2 type xy" 
echo "@xaxis label \"Iteration\"" 
echo "@xaxis label char size 0.6" 
echo "@xaxis ticklabel char size 0.6" 
echo "@yaxis label \"Convergence rate\"" 
echo "@yaxis label char size 0.6" 
echo "@yaxis ticklabel char size 0.6" 
echo "@legend string 0 \"first\"" 
echo "@legend string 1 \"second\"" 
echo "@legend string 2 \"infinite\"" 
echo "@target G2.S0" 
echo "@type xy" 
awk -f divergence_rate.awk < $1 | awk '{print $1 " " $3}' 
echo "&" 
echo "@target G2.S1" 
echo "@type xy" 
awk -f divergence_rate.awk < $1 | awk '{print $1 " " $4}' 
echo "&" 
echo "@target G2.S2" 
echo "@type xy" 
awk -f divergence_rate.awk < $1 | awk '{print $1 " " $5}' 
echo "&" 

echo "@with g3" 
echo "@g3 on" 
echo "@g3 hidden true"
