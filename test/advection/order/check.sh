#! /bin/sh

if test -z "$1"; then
    echo "usage: check.sh REFERENCES TOLERANCE"
    echo "  REFERENCES: reference xmgr files. example: \"reference/*.xmgr\""
    echo "  TOLERANCE:  1 no tolerance: 0: maximum tolerance"
    exit 1
fi

if test -e test; then
    :
else
    mkdir test;
fi

failed=0; all=0;
startdate=`date +%s`
for file in $1; do
    testfile=test/`basename $file`
    sname=`basename $file | awk '{print substr ($1, 1, index ($1, ".") - 1)}'`
    command=`awk '\
	BEGIN{ FS=" |\""; }\
	{\
	    if ($1 == "@description") {\
		for (i = 3; i <= NF; i++)\
		    printf ("%s ", $i);\
		exit 0;\
	    }\
	}\' < $file`
    maxlevel=`awk 'BEGIN{max=0} /^[[:space:]]+[[:digit:]|^.|-]+[[:space:]]+[[:digit:]|^.|-]+/ {if ($1 > max) max = $1;} END {print max}' < $file`
    ./order.sh "$command" $maxlevel > $testfile
    avg=`awk -f avgorder.awk < $testfile`
    avgref="$avg `awk -f avgorder.awk < $file`"
    if echo $avgref | awk -v tol=$2 -f checkorder.awk; then
	all=`expr $all + 1`;
	echo "PASS: $sname [$avg]"
    else
	all=`expr $all + 1`;
	failed=`expr $failed + 1`;
	echo "FAIL: $sname [$avg]"
    fi
    expr `date +%s` - $startdate > timestamp
done

if test "$failed" -eq 0; then
    banner="All $all tests passed";
else
    banner="$failed of $all tests failed";
fi
dashes=`echo "$banner" | sed s/./=/g`;
echo "$dashes";
echo "$banner";
echo "$dashes";
if test "$failed" -eq 0; then
    exit 0;
else
    exit 1;
fi
