awk 'BEGIN {FS = "\"| "} {
  if ($3 == "SYSTEM")
    print $5;
}' < gfs-docs.sgml | sort > /tmp/files
ls sgml/*.sgml > /tmp/files1
diff /tmp/files /tmp/files1 | awk '{
  if ($1 == ">")
    print "not in docs: " $2;
  else if ($1 == "<")
    print "in docs but not in sgml: " $2;
}' | sort




