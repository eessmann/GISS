version=$(date +%Y%m%d)
date=$(date +"%a, %e %b %Y %T %z")

cat <<EOF > debian/changelog
gfsview-snapshot ($version) jaunty; urgency=low

EOF

if test -z $1; then
    cat <<EOF >> debian/changelog
  * gfsview-snapshot release

EOF
else
    if html2text < $1 | \
       grep -v "^[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\} " | \
       sed 's/^\([^ ]\{1,\}\)/  \1/g'>> debian/changelog; then :
    else
	exit 1
    fi
fi

cat <<EOF >> debian/changelog
 -- Stephane Popinet <popinet@users.sf.net>  $date
EOF
