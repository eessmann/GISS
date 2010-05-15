#!/bin/sh

for f in $1/*.la; do
    lib=`tr -d "'" < $f | grep dlname | cut -d= -f2-3`
    nm -fb $1/.libs/$lib | grep ".* T gfs_.*_class$" | cut -d" " -f3-4
done | sort | uniq | sed -e 's/_class//g' -e 's/^./\U&/' -e 's/_./\U&/g' -e 's/_//g' | \
awk '{ print $0; gsub ("^Gfs", ""); print $0; }'
