#! /bin/sh
# compares two files containing columns of numbers
#
# 1: file1
# 2: file2
# 3: relative tolerance

max=`awk -v file=$2 '
BEGIN {
    max = 0.;
}
{
    for (i = 1; i <= NF; i++)
	val1[i] = $i;
    getline < file
    for (i = 1; i <= NF; i++) {
	v1 = $i - val1[i];
	v1 = v1 > 0. ? v1 : -v1;
	v2 = ($i + val1[i])/2.;
	v2 = v2 > 0. ? v2 : -v2;
	if (v2 > 1e-9 && v1/v2 > max)
	    max = v1/v2;
    }
}
END {print max}' < $1`

echo $max
if test `expr $max "<=" $3` = 1; then
    exit 1
fi

