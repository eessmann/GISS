#! /bin/bash
for radius in `seq 0. 0.001 0.3`; do
  rate=`./shapes star --dr $radius | ./poisson -v 6 5 4 -s - 2>&1 | awk -f divergence_rate.awk | awk '{ avg = $3; max = $4;}END{print avg " " max;}'`
  echo $radius $rate
done
