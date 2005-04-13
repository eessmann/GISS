#! /bin/sh
# $1: test files

if test -z "$1"; then
    echo "usage: report.sh REFERENCES"
    echo "  REFERENCES: reference xmgr files. example: \"reference/*.xmgr\""
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
\begin{document}
\begin{titlepage}
\mbox{}
\vspace{6cm}
\begin{center}
\Large
{\bf Tests of spatial order of advection scheme} \\\\
\vspace{5mm}
EOF

../advection -V 2>&1 | awk '{print $0 "\\\\"}' >> $texfile.tex
echo "`uname -a` \\\\" >> $texfile.tex
echo "Total running time: `cat timestamp` s\\\\" >> $texfile.tex

cat <<EOF >> $texfile.tex
\today \\\\
\end{center}
\end{titlepage}
\listoffigures
\newpage
EOF

for file in $1; do
    testfile=test/`basename $file`
    sname=`basename $file | awk '{print substr ($1, 1, index ($1, ".") - 1)}'`
    xmgrace -hardcopy -noask -hdevice EPS -printfile figures/$sname.eps -p order.par $file $testfile
    echo "\\begin{figure}" >> $texfile.tex
    echo "\\begin{center}" >> $texfile.tex
    echo "\\psfig{file=figures/$sname.eps, height=\\hsize, angle=270}" >> $texfile.tex
    echo "\\end{center}" >> $texfile.tex
    esname=`echo $sname | awk 'BEGIN{FS=""}{for (i = 1; i <= NF; i++)if($i=="_")printf("\\\_"); else printf("%s", $i);}'`
    command=`awk '\
	BEGIN{ FS=" |\""; }\
	{\
	    if ($1 == "@description") {\
		for (i = 3; i <= NF; i++)\
		    printf ("%s ", $i);\
		exit 0;\
	    }\
	}\' < $file`
    echo "\\caption{$esname: {\tt $command}}" >> $texfile.tex
    echo "\\end{figure}" >> $texfile.tex
    echo "\\clearpage" >> $texfile.tex
done
echo "\\end{document}" >> $texfile.tex
latex -interaction=nonstopmode $texfile.tex > /dev/null 2>&1
latex -interaction=nonstopmode $texfile.tex > /dev/null 2>&1
dvips $texfile.dvi -o $texfile.ps > /dev/null 2>&1
rm -f $texfile.log $texfile.aux $texfile.dvi
