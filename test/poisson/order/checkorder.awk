{
  n = NF/2;
  if (n != 6)
    exit 1;
  for (i = 1; i <= n; i++)
    if ($i < $(i + n)*tol)
      exit 1;
  exit 0;
}
