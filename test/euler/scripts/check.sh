#! /bin/sh

if test -z "$1"; then
    echo "usage: check.sh REFERENCES"
    echo "  REFERENCES: reference files. example: \"reference/*\""
    exit 1
fi

if test -e test; then
    :
else
    mkdir test;
fi

PATH=$PATH:scripts
failed=0; all=0;
startdate=`date +%s`
for file in $1; do
    testfile=test/`basename $file`
    sname=`basename $file | awk '{print substr ($1, 1, index ($1, ".") - 1)}'`
    ext=`basename $file | awk '{print substr ($1, index ($1, ".") + 1)}'`
    command=""

    if test "$ext" = "xmgr"; then
	command=`awk '
	    BEGIN {
		com = ""
	    }
	    {
		if ($1 == "@description") {
		    for (i = 2; i <= NF; i++)
			com = com " " $i;
		    print com
		    exit 0;
		}
	    }' < $file`
	command=`echo $command | sed 's/\\\"/"/g' | sed 's/^"//g' | sed 's/"$//g'`
     elif test "$ext" = "tex"; then
	command=`awk '{
		    if ($2 == "command:") {
			for (i = 3; i <= NF; i++)
			    printf ("%s ", $i);
			exit 0;
		    }
		}' < $file`
     fi

     if test -n "$command"; then
	case $command in
	    *stationary.sh* | *order.sh* | *time.sh* )
		command="$command $file"
		;;	    
	esac
	if /bin/sh -c "$command" > $testfile; then
	    all=`expr $all + 1`;
	    echo "PASS: $sname"
	else
	    all=`expr $all + 1`;
	    failed=`expr $failed + 1`;
	    echo "FAIL: $sname"
	fi
	expr `date +%s` - $startdate > timestamp
     fi
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
