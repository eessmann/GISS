#! /bin/sh

if test -z $1; then
    echo "usage: ./monitor.sh OUTPUT"
    exit 1
fi

echo "@WITH G0"
echo "@G0 ON"
echo "@G0 TYPE logy"
echo "@XAXIS LABEL \"Time\""
echo "@xaxis label char size 0.6"
echo "@xaxis ticklabel char size 0.6"
echo "@YAXIS LABEL \"Divergence\""
echo "@yaxis label char size 0.6"
echo "@yaxis ticklabel char size 0.6"
echo "@LEGEND STRING 0 \"L1 norm\""
echo "@LEGEND STRING 1 \"L2 norm\""
echo "@LEGEND STRING 2 \"Lmax norm\""

echo "@TARGET S0"
echo "@TYPE xy"
awk '{if ($2 == "time:") time = $3; else if ($1 == "div" && $2 == "after") print time " " $4;}' < $1
echo "&"

echo "@TARGET S1"
echo "@TYPE xy"
awk '{if ($2 == "time:") time = $3; else if ($1 == "div" && $2 == "after") print time " " $6;}' < $1
echo "&"

echo "@TARGET S2"
echo "@TYPE xy"
awk '{if ($2 == "time:") time = $3; else if ($1 == "div" && $2 == "after") print time " " $8;}' < $1
echo "&"

echo "@WITH G1"
echo "@G1 ON"
echo "@G1 TYPE xy"
echo "@XAXIS LABEL \"Time\""
echo "@xaxis label char size 0.6"
echo "@xaxis ticklabel char size 0.6"
echo "@YAXIS LABEL \"Kinetic Energy\""
echo "@yaxis label char size 0.6"
echo "@yaxis ticklabel char size 0.6"

echo "@TARGET S0"
echo "@TYPE xy"
awk '{if ($2 == "time:") time = $3; else if ($1 == "velocity") print time " " $5*$5;}' < $1
echo "&"

echo "@WITH G2"
echo "@G2 ON"
echo "@G2 TYPE xy"
echo "@XAXIS LABEL \"Time\""
echo "@xaxis label char size 0.6"
echo "@xaxis ticklabel char size 0.6"
echo "@YAXIS LABEL \"Vorticity\""
echo "@yaxis label char size 0.6"
echo "@yaxis ticklabel char size 0.6"
echo "@LEGEND STRING 0 \"L1 norm\""
echo "@LEGEND STRING 1 \"L2 norm\""
echo "@LEGEND STRING 2 \"Lmax norm\""

echo "@TARGET S0"
echo "@TYPE xy"
awk '{if ($2 == "time:") time = $3; else if ($1 == "vorticity") print time " " $3;}' < $1
echo "&"

echo "@TARGET S1"
echo "@TYPE xy"
awk '{if ($2 == "time:") time = $3; else if ($1 == "vorticity") print time " " $5;}' < $1
echo "&"

echo "@TARGET S2"
echo "@TYPE xy"
awk '{if ($2 == "time:") time = $3; else if ($1 == "vorticity") print time " " $7;}' < $1
echo "&"
