#! /bin/sh

if test -z "$1"; then
    echo "usage: xmgrmerge.sh FILE.xmgr FILE.par"
    exit 1
fi

echo "# ACE/gr parameter file"
echo "#"
awk '{ if (substr ($0, 1, 1) != "#" && $1 != "description") 
          print "@" $0; }' < $2
awk '{ if (substr ($0, 1, 1) != "#") print $0; }' < $1

