#! /bin/bash

for cycle in `seq 1 1 20`; do
    ./simple 8 $cycle 4 4 2>&1 | awk '{
	if ($1 == "div")
	    max = $7;
	else if ($1 == "Time:")
	    print $2 " " max;
    }'
done
