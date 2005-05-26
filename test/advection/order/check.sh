#! /bin/sh

if test -z "$1"; then
    echo "usage: check.sh DIRECTORIES ORDER"
    exit 1
fi

failed=0; all=0;
startdate=`date +%s`
for file in $1; do
    sname=`basename $file`
    if sh order.sh 3 8 $file/advection.gfs $2 > $file/result.xmgr; then
	all=`expr $all + 1`;
	echo "PASS: $sname"
    else
	all=`expr $all + 1`;
	failed=`expr $failed + 1`;
	echo "FAIL: $sname"
    fi
done
expr `date +%s` - $startdate > timestamp

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
