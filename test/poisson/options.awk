BEGIN {
  FS = " |\t|=";
}
{
  n = 0;
  for (i = 1; i <= NF; i++)
    if (substr($i,1,1) == "-")
      o[n++] = $i;
  if (n == 2)
    printf ("s/%s/%s/g ", o[1], o[0]);
}
