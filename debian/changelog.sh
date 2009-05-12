GERRIS_VERSION=$(grep GFS_VERSION config.h | awk 'BEGIN{FS=" |\""}{print $4}')
version=$(grep GFS_BUILD_VERSION src/version.h | awk 'BEGIN{FS=" |\"|-"}{print $4}')
date=$(date +"%a, %e %b %Y %T %z")

if test -z $1; then
    cat <<EOF > debian/changelog
gerris-snapshot ($GERRIS_VERSION-$version) hardy; urgency=low

  * gerris-snapshot release (based on Drew Parsons's official debian)

 -- Stephane Popinet <popinet@users.sf.net>  $date
EOF
else
    cat <<EOF > debian/changelog
gerris-snapshot ($GERRIS_VERSION-$version) hardy; urgency=low

EOF
    if html2text < $1 | grep -v "^[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\} " >> debian/changelog; then :
    else
	exit 1
    fi
cat <<EOF >> debian/changelog
 -- Stephane Popinet <popinet@users.sf.net>  $date
EOF
fi
