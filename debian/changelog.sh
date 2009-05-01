GERRIS_VERSION=$(grep GFS_VERSION config.h | awk 'BEGIN{FS=" |\""}{print $4}')
version=$(grep GFS_BUILD_VERSION src/version.h | awk 'BEGIN{FS=" |\"|-"}{print $4}')
date=$(date +"%a, %e %b %Y %T %z")
cat <<EOF > debian/changelog
gerris-snapshot ($GERRIS_VERSION-$version) testing; urgency=low

  * gerris-snapshot release (based on Drew Parsons's official debian)

 -- Stephane Popinet <popinet@users.sf.net>  $date
EOF
