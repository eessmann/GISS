BEGIN {
  ingraph = 0;
  indata = 0;
  sum = 0.;
  n = 0;
}
{
  if ($1 == "@target") {
    graph = substr ($2, 1, 2);
    if (graph == "G0")
      ingraph = 0;
    else if (graph == "G1")
      ingraph = 0;
    else if (graph == "G2")
      ingraph = 1;
    else if (graph == "G3")
      ingraph = 0;
  }
  else if (ingraph) {
    if ($1 == "@type") {
      if (n > 0) {
	avgrate = sum/n;
	printf ("%g ", avgrate);
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
    avgrate = sum/n;
    printf ("%g\n", avgrate);
  }
  else
    exit 1;
  exit 0;
}
