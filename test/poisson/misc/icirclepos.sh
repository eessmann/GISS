#! /bin/bash
for pos in `seq -0.25 0.001 0.25`; do
  rate=`./shapes ellipse | transform --revert --tx $pos | ./poisson -v 6 5 4 -s - 2>&1 | awk -f divergence_rate.awk | awk '{ avg = $3; max = $4;}END{print avg " " max;}'`
  echo $pos $rate
done
