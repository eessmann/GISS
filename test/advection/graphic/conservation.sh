#!/bin/sh

if awk '
BEGIN {
  mass0 = 0.;
}
{
  if ($4 == "sum:") {
    if (mass0 == 0.)
      mass0 = $5;
    if (fabs ($5 - mass0) > 1e-10)
      exit 1;
  }
#  else if ($4 == "min:") {
#    if ($5 < -1e-2 || $11 > 1. + 1e-2)
#      exit 1;
#  }
}' < log; then
  exit 0;
else
  exit 1;
fi
