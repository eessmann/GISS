#! /bin/sh

if test -z "$1"; then
    echo "usage: figures.sh FILES VAR [COLUMNS]"
    echo "  FILES: GFS files. example: \"*.gfs\""
    echo "  VAR: GFS variable"
    echo "  [COLUMNS]: number of columns (default is 3)"
    exit 1;
fi
var=$2

if test -z "$3"; then
    ncol=2
else
    ncol=`expr $3 - 1`
fi

x=0
y=0
col=0
figures=`mktemp /tmp/figures.XXXXXX`
psfig=`mktemp /tmp/psfig.XXXXXX`
echo "LIST {" > $figures
for file in $1; do
    ../../tools/gfs2other -g $var < $file | transform --tx $x --ty $y | gts2oogl --isolines 10 --nosurface --flatten >> $figures
    x1=`awk -v a=$x 'BEGIN {print a - 0.5}'`
    y1=`awk -v a=$y 'BEGIN {print a - 0.5}'`
    x2=`awk -v a=$x 'BEGIN {print a + 0.5}'`
    y2=`awk -v a=$y 'BEGIN {print a + 0.5}'`
    echo "VECT 1 5 0 5 0 $x1 $y1 0 $x1 $y2 0 $x2 $y2 0 $x2 $y1 0 $x1 $y1 0" >> $figures
    if expr $col \< $ncol >/dev/null; then
	col=`expr $col + 1`
	x=`expr $x + 1`
    else
	col=0
        x=0
        y=`expr $y - 1`
    fi
done
echo "}" >> $figures
geomview -nopanels -wins 0 -c - <<EOF
    (camera default { 
			perspective 0 
			fov 1.5
		    })
    (backcolor default 1. 1. 1.)
    (load $figures geometry)
    (bbox-draw World no)
    (normalization World each)
    (look)
    (snapshot default "$psfig" ps 300 300)
    (exit)
EOF
rm -f $figures
sed 's/1 setlinewidth/0.1 setlinewidth/g' < $psfig
rm -f $psfig
