#! /bin/sh

if test -z "$1"; then
    echo "usage: figures.sh FILES [COLUMNS]"
    echo "  FILES: GTS files. example: \"*.gts\""
    echo "  [COLUMNS]: number of columns (default is 3)"
    exit 1;
fi

if test -z "$2"; then
    ncol=2
else
    ncol=`expr $2 - 1`
fi

x=0
y=0
col=0
echo "LIST {" > /tmp/figures.oogl
for file in $1; do
    transform --tx $x --ty $y < $file | gts2oogl --isolines 10 --nosurface | awk '{if (NF == 3) print $1 " " $2 " " 0; else print $0;}' >> /tmp/figures.oogl
    x1=`awk -v a=$x 'BEGIN {print a - 0.5}'`
    y1=`awk -v a=$y 'BEGIN {print a - 0.5}'`
    x2=`awk -v a=$x 'BEGIN {print a + 0.5}'`
    y2=`awk -v a=$y 'BEGIN {print a + 0.5}'`
    echo "VECT 1 5 0 5 0 $x1 $y1 0 $x1 $y2 0 $x2 $y2 0 $x2 $y1 0 $x1 $y1 0" >> /tmp/figures.oogl
    if expr $col \< $ncol >/dev/null; then
	col=`expr $col + 1`
	x=`expr $x + 1`
    else
	col=0
        x=0
        y=`expr $y - 1`
    fi
done
echo "}" >> /tmp/figures.oogl
geomview -nopanels -wins 0 -c - <<EOF
    (camera default { 
			perspective 0 
			fov 1.42217
		    })
    (backcolor default 1. 1. 1.)
    (load /tmp/figures.oogl geometry)
    (bbox-draw World no)
    (normalization World each)
    (look)
    (snapshot default "/tmp/figures.ps" ps 300 300)
    (exit)
EOF
rm -f /tmp/figures.oogl
sed 's/1 setlinewidth/0.1 setlinewidth/g' < /tmp/figures.ps
rm -f /tmp/figures.ps
