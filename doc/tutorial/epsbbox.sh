#! /bin/bash
for file in $1; do
    if ps2eps -f -q -B -P $file; then
      mv -f $file.eps $file
    fi
done
