#! /bin/sh

if test -z "$1"; then
    echo "usage: check.sh DIRECTORIES"
    exit 1
fi

failed=0; all=0;
startdate=`date +%s`
PATH=$PATH:../../../../poisson:../../..
for test in $1; do
    cd $test
    command=`awk '{if ($1 == "command:") print substr ($0, 10);}' < description.txt`
    rm -f *.gts
    if /bin/sh -c "$command > log 2>&1"; then
	all=`expr $all + 1`;
	echo "PASS: `basename $test`"
    else
	all=`expr $all + 1`;
	failed=`expr $failed + 1`;
	echo "FAIL: `basename $test`"
    fi
    cd ../..
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
