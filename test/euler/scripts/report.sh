#! /bin/sh
# $1: test files

if test -z "$1"; then
    echo "usage: report.sh REFERENCES"
    echo "  REFERENCES: reference files. example: \"reference/*\""
    exit 1;
fi

if test -e figures; then 
    :
else
    mkdir figures;
fi

texfile=report
cat <<EOF > $texfile.tex
\documentclass[a4paper,11pt]{article}
\usepackage{epsfig}

\oddsidemargin=0mm
\evensidemargin=-5mm
\topmargin=-7mm
\textwidth=16.22cm
\textheight=23.2cm

\renewcommand\floatpagefraction{.9}
\renewcommand\topfraction{.9}
\renewcommand\bottomfraction{.9}
\renewcommand\textfraction{.1}   
\setcounter{totalnumber}{50}
\setcounter{topnumber}{50}
\setcounter{bottomnumber}{50}

\begin{document}
\begin{titlepage}
\mbox{}
\vspace{6cm}
\begin{center}
\Large
{\bf Tests of Euler solver} \\\\
\vspace{5mm}
EOF

gerris2D -V 2>&1 | awk '{print $0 "\\\\"}' >> $texfile.tex
echo "`uname -a` \\\\" >> $texfile.tex
echo "Total running time: `cat timestamp` s\\\\" >> $texfile.tex

cat <<EOF >> $texfile.tex
\today \\\\
\end{center}
\end{titlepage}
\tableofcontents
\listoffigures
\newpage
EOF

for file in $1; do
    testfile=test/`basename $file`
    sname=`basename $file | awk '{print substr ($1, 1, index ($1, ".") - 1)}'`
    ext=`basename $file | awk '{print substr ($1, index ($1, ".") + 1)}'`
    command=""

    if test "$ext" = "xmgr"; then
	command=`awk '
	    BEGIN {
		com = ""
	    }
	    {
		if ($1 == "@description") {
		    for (i = 2; i <= NF; i++)
			com = com " " $i;
		    print com
		    exit 0;
		}
	    }' < $file`
	command=`echo $command | sed 's/\\\"/"/g' | sed 's/^"//g' | sed 's/"$//g'`
	params=`echo $command | awk '{
	    print substr ($1, 1, index ($1, ".sh") - 1);
        }'`
	xmgr -hardcopy -noask -eps -device 2 -printfile figures/$sname.eps -p parameters/$params.par $file $testfile > /dev/null 2>&1
	echo "\\begin{figure}" >> $texfile.tex
	echo "\\begin{center}" >> $texfile.tex
	echo "\\psfig{file=figures/$sname.eps, height=\\hsize, angle=270}" >> $texfile.tex
	echo "\\end{center}" >> $texfile.tex
	esname=`echo $sname | awk 'BEGIN{FS=""}{for (i = 1; i <= NF; i++)if($i=="_")printf("\\\_"); else printf("%s", $i);}'`
	echo "\\caption{$esname: {\tt $command}}" >> $texfile.tex
	echo "\\end{figure}" >> $texfile.tex
	echo "\\clearpage" >> $texfile.tex
     elif test "$ext" = "tex"; then
	command=`awk '{
		    if ($2 == "command:") {
			for (i = 3; i <= NF; i++)
			    printf ("%s ", $i);
			exit 0;
		    }
		}' < $file`
	legend=`awk '{
		    if ($2 == "legend:") {
			for (i = 3; i <= NF; i++)
			    printf ("%s ", $i);
			exit 0;
		    }
		}' < $file`
	echo "\\section{$legend {\\tt $command}}" >> $texfile.tex
	echo "\\subsection*{Reference}" >> $texfile.tex	
	cat $file >> $texfile.tex
	echo "\\subsection*{Test}" >> $texfile.tex
	if test -f $testfile; then
	    cat $testfile >> $texfile.tex
	else
	    echo "This test did not run." >> $texfile.tex
	fi
	echo "\\clearpage" >> $texfile.tex
     fi
done
echo "\\end{document}" >> $texfile.tex
latex -interaction=nonstopmode $texfile.tex > /dev/null 2>&1
latex -interaction=nonstopmode $texfile.tex > /dev/null 2>&1
dvips $texfile.dvi -o $texfile.ps > /dev/null 2>&1
rm -f $texfile.log $texfile.aux $texfile.dvi
