#!/bin/sh

# do not forget to update ../tools/darcs2dist when changing the way $version is computed
version=`darcs changes --last=1 --xml-output | \
    awk 'BEGIN{RS=" ";FS="="}{if ($1 == "date") print substr($2,4,6) "-" substr($2,10,6);}'`
changes=`darcs whatsnew -s | awk '{
  if ($1 == "M" && substr($2,1,6) != "./doc/") {
    print " + local changes";
    exit 0;
  }
}'`
version="$version$changes"

if test -f version.h ; then
    oldversion=`awk '{if ($2 == "GFS_BUILD_VERSION") print $0;}' < version.h`
fi

if [ "x$oldversion" != "x#define GFS_BUILD_VERSION \"$version\"" ] ; then
    cat <<EOF > version.h
/* version.h
 *
 * This is a generated file.  Please modify 'darcsversion.sh'
 */

#ifndef GFSVERSION_H
#define GFSVERSION_H

#define GFS_BUILD_VERSION "$version"

#endif /* GFSVERSION_H */
EOF
fi
