#! /bin/sh
# $1: command line
# $2: options

if test -z "$1"; then
    echo "usage: runfig.sh COMMAND [OPTIONS]"
    exit 1;
fi

/bin/sh -c "$1 > /tmp/tutu 2>&1"
if test $2; then
    ./divfig.sh /tmp/tutu "$1" > /tmp/divfig
    xmgr -graph 0 -autoscale xy -graph 1 -autoscale xy -graph 2 -autoscale xy /tmp/divfig -p divfig.par
else
    ./divfig.sh /tmp/tutu "$1"
fi
