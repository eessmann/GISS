BEGIN {
  ingraph = 0;
  indata = 0;
  sum = 0.;
  n = 0;
}
{
  if ($1 == "@WITH") {
    if ($2 == "G0")
      ingraph = 0;
    else if ($2 == "G1")
      ingraph = 1;
    else if ($2 == "G2")
      ingraph = 0;
    else if ($2 == "G3")
      ingraph = 1;
  }
  else if (ingraph) {
    if ($1 == "@TYPE") {
      if (n > 0) {
	avgorder = sum/n;
	printf ("%g ", avgorder);
      }
      indata = 1;
      sum = 0.;
      n = 0;
    }
    else if ($1 == "&")
      indata = 0;
    else if (indata) {
      sum += $2;
      n++;
    }
  }
}
END {
  if (n > 0) {
    avgorder = sum/n;
    printf ("%g\n", avgorder);
  }
  else
    exit 1;
  exit 0;
}
