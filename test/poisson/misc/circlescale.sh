#! /bin/bash
for scale in `seq 1 0.01 2`; do
  rate=`./shapes ellipse | transform --scale $scale | ./poisson -v 6 5 4 -s - 2>&1 | awk -f divergence_rate.awk | awk '{ avg = $3; max = $4;}END{print avg " " max;}'`
  echo $scale $rate
done
