{
  if ($1 == "total" && $2 == "err") {
    totalfirst = $6;
    totalsecond = $8;
    totalinfty = $10;
  }
  else if ($1 == "refined" && $2 == "err")
    print $4 " " totalfirst " " totalsecond " " totalinfty " " $6 " " $8 " " $10;
}
