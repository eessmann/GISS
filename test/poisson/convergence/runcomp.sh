#! /bin/sh
# $1: xmgr file 
# $2: options

if test -z "$1"; then
    echo "usage: runcomp.sh REFERENCE [OPTIONS]"
    echo "  REFERENCE: an xmgr reference file"
    exit 1;
fi

command=`awk -f ../command.awk < $1`

options=`transform -h 2>&1 | awk -f ../options.awk`
options="$options `../shapes -h 2>&1 | awk -f ../options.awk`"
options="$options `../poisson -h 2>&1 | awk -f ../options.awk`"
for replace in $options; do
    command=`echo $command | sed $replace`;
done

./runfig.sh "$command" $2
