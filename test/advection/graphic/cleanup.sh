#! /bin/sh

for file in tests/*; do
    rm -f $file/*.gts $file/*.ps $file/log $file/*~ $file/core
done
rm -f timestamp
