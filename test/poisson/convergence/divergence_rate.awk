BEGIN {
  n = 0;
}
{
  if ($1 == "div") {
    iter[n] = $3;
    time[n] = $5;
    first[n] = $7;
    second[n] = $9;
    infty[n++] = $11;
  }
}
END {
  for (i = 0; i < n - 1; i++)
    if (first[i+1] != 0. && second[i+1] != 0. && infty[i+1] != 0.)
      print iter[i+1] " " time[i+1] " " exp(log(first[0]/first[i+1])/(i + 1)) " " exp(log(second[0]/second[i+1])/(i + 1)) " " exp(log(infty[0]/infty[i+1])/(i + 1));
}
