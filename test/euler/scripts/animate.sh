#! /bin/sh

if test -z "$1"; then
    echo "usage: figures.sh FILES VAR [SOLID] [OPTIONS]"
    echo "  FILES: GFS files. example: \"*.gfs\""
    echo "  VAR: GFS variable"
    echo "  SOLID: a GTS file"
    echo "  OPTIONS: to pass to gts2oogl (default is \"--height --flatten\")"
    exit 1;
fi
var=$2

figures=`mktemp /tmp/figures.XXXXXX`
image=`mktemp /tmp/image.XXXXXX`
solid=`mktemp /tmp/solid.XXXXXX`

if test -z "$3"; then
    :
else
    gts2oogl < $3 > $solid
fi

if test -z "$4"; then
#    options="--height --flatten"
    options=""
else
    options=$4
fi

for file in $1; do
#    ../../tools/gfs2other -g $var < $file | gts2oogl $options > $figures
    ../../tools/gfs2other -S $var $options < $file > $figures
    geomview -nopanels -wins 0 -c - <<EOF
    (camera default { 
			perspective 0 
			fov 1.0
		    })
    (backcolor default 1. 1. 1.)
    (bbox-draw World no)
    (normalization World each)
    (merge-baseap appearance {shading csmooth})
    (load $figures geometry)
    (load $solid geometry)
    (look)
    (snapshot default "$image" ppm 1000 1000)
    (exit)
EOF
    convert -crop 0x0 -geometry 600x600 -comment "$file" $image miff:-
done
rm -f $figures $image $solid
