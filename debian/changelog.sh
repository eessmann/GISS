version=`awk '{if ($1 == "Version:") print $2;}' < src/gerris2D.pc`
date=`date "+%a, %e %b %Y %T %z"`

cat <<EOF > debian/changelog
gerris-snapshot ($version) testing; urgency=low

  * gerris-snapshot release (based on Marcelo's official debian)

 -- Stephane Popinet <popinet@users.sf.net>  $date
EOF
