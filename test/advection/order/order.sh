#!/bin/sh

usage()
{
	cat <<EOF
Usage: order LMIN LMAX FILE ORDER

Generates an xmgr file for the order of convergence of simulation FILE.

EOF
	exit $1
}

if test $# -lt 4; then
	usage 1 1>&2
fi

if for level in `seq $1 1 $2`; do
  sed "s/LEVEL/$level/g" < $3 | gerris2D -
done | awk -v order=$4 '
BEGIN { n = 0 }
{
  l[n] = $1
  l1[n] = $2
  l2[n] = $3
  li[n++] = $4
}
END {
  print "@target G0.S3"
  print "@type xy"
  for (i = 0; i < n; i++)
    print l[i] " " l1[i];
  print "@target G0.S4"
  print "@type xy"
  for (i = 0; i < n; i++)
    print l[i] " " l2[i];
  print "@target G0.S5"
  print "@type xy"
  for (i = 0; i < n; i++)
    print l[i] " " li[i];

 print "@target G1.S3"
  print "@type xy"
  for (i = 1; i < n; i++)
    print l[i] " " log(l1[i-1]/l1[i])/log(2);
  print "@target G1.S4"
  print "@type xy"
  for (i = 1; i < n; i++)
    print l[i] " " log(l2[i-1]/l2[i])/log(2);
  print "@target G1.S5"
  print "@type xy"
  for (i = 1; i < n; i++)
    print l[i] " " log(li[i-1]/li[i])/log(2);

  if (log(l1[n-2]/l1[n-1])/log(2) < order || 
      log(l2[n-2]/l2[n-1])/log(2) < order || 
      log(li[n-2]/li[n-1])/log(2) < order)
    exit (1);
  else
    exit (0);  
}'; then
    exit 0
else
    exit 1
fi

