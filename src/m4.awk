function replace_params(s, b,    i)
{
    for (i in b)
	gsub(b[i], "($" i ")", s);
    return s;
}

BEGIN {
    print "m4_changecom()m4_dnl";
}
{
    if ($1 == "GfsDefine" || $1 == "Define") {
	macro = $2;
	delete b;
	if (match(macro, /(.+)\((.+)\)/, a)) {
	    macro = a[1];
	    split(a[2],b,",");
	}
	printf ("m4_define(`%s',`%s", macro, replace_params($3, b));
	for (i = 4; i <= NF; i++)
	    printf (" %s", replace_params($i, b));
	printf ("')\n");
    }
    else if ($1 == "GfsInclude" || $1 == "Include")
	printf ("m4_include(%s)\n", $2);
    else
	print $0;
}
