#! /bin/sh
# $1: command line
# $2: options

if test -z "$1"; then
    echo "usage: runorderfig.sh COMMAND [OPTIONS]"
    exit 1;
fi

/bin/sh -c "$1 2>&1 | awk -f order1.awk > /tmp/orderfig"
if test $2; then
    ./orderfig.sh /tmp/orderfig "$1" > /tmp/order
    xmgr -graph 0 -autoscale xy -graph 1 -autoscale x -graph 2 -autoscale xy -graph 3 -autoscale x /tmp/order -p orderfig.par
else
    ./orderfig.sh /tmp/orderfig "$1"
fi
