# one argument: name = variable name
{
    if ($1 == "#") {
	for (i = 2; i <= NF; i++) {
            if (match($i,"(.*):([A-Z]+)([0-9]+)",a) && a[2] == name) {
		layer[a[3]] = int(a[1]);
		if (int(a[3]) > nl)
		    nl = int(a[3]);
	    }
            else if (match($i,"(.*):([A-Z]+)",a) && a[2] == "Zb")
		zb = int(a[1]);
	}
	nl++
    }
    else {
	dz = $4/nl;
	print $1,$zb,$layer[0]/dz
	for (i = 0; i < nl; i++) {
	    z = dz*(0.5+i)
	    print $1,$zb+z,$layer[i]/dz
	}
	print $1,$zb+dz*nl,$layer[nl - 1]/dz
	print ""
    }
}
