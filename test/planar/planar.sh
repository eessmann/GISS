if test x$donotrun != xtrue; then
    awk 'BEGIN{R1 = -0.2; R2=0.3 ; h=(R2-R1)/40. ; for (x = R1; x <= R2; x += h) print 0,x,0;}' > yprofile
    for level in 3 4; do
	if gerris2D -DLEVEL=$level planar.gfs; then :
        mv xplanar xplanar-$level
	else
	    exit 1
	fi
    done
fi



