#! /bin/sh

if test -z "$1"; then
    command="./periodic -l 5"
    tolerance=1e-4
else
    command=$1
    tolerance=$2
fi

log=`mktemp /tmp/log.XXXXXX`
log1=`mktemp /tmp/log1.XXXXXX`

$command 2>&1 | awk '{if ($1 == "error") print $3 " " $5 " " $7 " " $9;}' > $log
mpirun -np 1 $command -M 2>&1 | awk '{if ($1 == "error") print $3 " " $5 " " $7 " " $9;}' > $log1
if cat $log $log1 | sort -k 1,2 | awk '{printf ("%s ", $0); getline; print $0;}' | awk -v tolerance=$tolerance '{
    er = $4 - $8; 
    if (er < 0.)
	er = -er; 
    if (er > tolerance)
	exit 1;
}'; then
    failed=""
else
    failed=1
fi

rm -f $log $log1
if test $failed; then
    exit 1;
fi
