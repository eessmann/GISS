BEGIN {
  n = 0;
}
{
  level[n] = $1;
  totalfirst[n] = $2;
  totalsecond[n] = $3;
  totalinfty[n] = $4;
  reffirst[n] = $5;
  refsecond[n] = $6;
  refinfty[n++] = $7;
}
END {
  for (i = 1; i < n; i++)    
    if (totalfirst[i] > 0. && totalsecond[i] > 0. && totalinfty[i] > 0. && reffirst[i] > 0. && refsecond[i] > 0. && refinfty[i] > 0.)
      print level[i] " " log(totalfirst[i-1]/totalfirst[i])/log(2) " " log(totalsecond[i-1]/totalsecond[i])/log(2) " " log(totalinfty[i-1]/totalinfty[i])/log(2) " " log(reffirst[i-1]/reffirst[i])/log(2) " " log(refsecond[i-1]/refsecond[i])/log(2) " " log(refinfty[i-1]/refinfty[i])/log(2)
}
