#! /bin/bash
for pos in `seq -0.1 0.001 0.1`; do
  rate=`./shapes 4ellipses | transform --ty $pos | ./poisson -v 6 5 4 -s - 2>&1 | awk -f divergence_rate.awk | awk '{ avg = $3; max = $4;}END{print avg " " max;}'`
  echo $pos $rate
done
