#! /bin/bash

for file in tutorial/*.html; do
awk '{
    if ($1 == "<PRE>") {
	inpre = 1;
	npre = 0;
    }
    else if ($1 == "</PRE>")
	inpre = 0;
    if (inpre) {
	if (NF > 0 || npre > 1)
	    print $0;
	npre++;
    }
    else
	print $0;
}' < $file | sed 's/<TT> /<TT>/g' > /tmp/`basename $file`
awk 'BEGIN {
    RS = "\"";
}
{
    if (dxscreendone == 0) {
	    if ($1 == "dx-screen")
		tagfound = 1;
    	else if (tagfound && $1 == "HREF=")
        	hreffound = 1;
    	else if (hreffound) {
		print "cp -f dxscreen.png tutorial/" $1;
		system ("cp -f dxscreen.png tutorial/" $1);
        	dxscreendone = 1;
    	}
    }
    if (gfsviewdone == 0) {
    	if ($1 == "gfsview")
		gfsvtagfound = 1;
    	else if (gfsvtagfound && $1 == "HREF=")
        	gfsvhreffound = 1;
    	else if (gfsvhreffound) {
		print "cp -f gfsview.png tutorial/" $1;
		system ("cp -f gfsview.png tutorial/" $1);
        	gfsviewdone = 1;
    	}
    }
    if (gfsviewdone && dxscreendone)
      	exit (0);
}' < $file
mv -f /tmp/`basename $file` $file
done
cp vorticity.png tutorial/img3.png
