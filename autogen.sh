#!/bin/sh

libtoolize --force
aclocal
autoheader
automake -a
autoconf
command=./configure
while [ $# -ne 0 ]; do
        command="$command $1"
        shift
done
$command
