awk 'BEGIN {FS = ">|<"} {
  if ($2 == "FILE")
    print "tmpl/" $3 ".sgml"; 
}' < gfs-sections.txt | sort > /tmp/files

ls tmpl/*.sgml | sort > /tmp/files1
diff /tmp/files /tmp/files1 | awk '{
  if ($1 == ">")
    print "not in sections: " $2;
}' | sort

