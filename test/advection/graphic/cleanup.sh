#! /bin/sh

for file in tests/*; do
    rm -f $file/*.eps $file/log $file/*~ $file/core
done
rm -f timestamp
