#! /bin/sh

if test -z "$1"; then
    echo "usage: report.sh DIRECTORIES"
    exit 1;
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
{\bf Graphical tests of advection scheme} \\\\
\vspace{5mm}
EOF

gerris2D -V 2>&1 | awk '{print $0 "\\\\"}' >> $texfile.tex
echo "`uname -a` \\\\" >> $texfile.tex
echo "Total running time: `cat timestamp` s\\\\" >> $texfile.tex

cat <<EOF >> $texfile.tex
\today \\\\
\end{center}
\end{titlepage}
\listoffigures
\newpage
EOF

figure=1
for test in $1; do
    name=`basename $test`
    if test "$name" != "CVS"; then
	cd $test
        esname=`echo $name | awk 'BEGIN{FS=""}{for (i = 1; i <= NF; i++)if($i=="_")printf("\\\_"); else printf("%s", $i);}'`
	caption=`awk '{if ($1 == "caption:") print substr ($0, 10);}' < description.txt`
	printf "Figure %2d: %s\n" $figure $name
	figure=`expr $figure + 1`
	cd ../..
        cat <<EOF >> $texfile.tex
\\begin{figure}
\\begin{center}
\\begin{tabular}{cc}
\\psfig{file=$test/t-0.eps, width=0.45\\hsize} &
\\psfig{file=$test/t-1.eps, width=0.45\\hsize} \\\\
\\psfig{file=$test/t-2.eps, width=0.45\\hsize} &
\\psfig{file=$test/t-3.eps, width=0.45\\hsize}
\\end{tabular}
\\end{center}
\\caption{$esname: $caption}
\\end{figure}
\\clearpage
EOF
    fi
done
echo "\\end{document}" >> $texfile.tex
latex -interaction=nonstopmode $texfile.tex > /dev/null 2>&1
latex -interaction=nonstopmode $texfile.tex > /dev/null 2>&1
dvips $texfile.dvi -o $texfile.ps > /dev/null 2>&1
rm -f $texfile.log $texfile.aux $texfile.dvi $texfile.lof
