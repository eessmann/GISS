BEGIN{
  FS=" |\"";
}
{
  if ($1 == "@description") {
    for (i = 3; i <= NF; i++) 
      printf ("%s ", $i); 
    exit 0;
  }
}
