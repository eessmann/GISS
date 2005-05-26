#! /bin/sh
# $1: test files

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
{\bf Tests of spatial order of advection scheme} \\\\
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

for file in $1; do
    xmgrace -hardcopy -noask -hdevice EPS -printfile $file/result.eps -p order.par $file/reference.xmgr $file/result.xmgr
    echo "\\begin{figure}" >> $texfile.tex
    echo "\\begin{center}" >> $texfile.tex
    echo "\\psfig{file=$file/result.eps, width=0.8\\hsize}" >> $texfile.tex
    echo "\\end{center}" >> $texfile.tex
    esname=`basename $file | awk 'BEGIN{FS=""}{for (i = 1; i <= NF; i++)if($i=="_")printf("\\\_"); else printf("%s", $i);}'`
    echo "\\caption{$esname: " >> $texfile.tex
    cat $file/description.tex >> $texfile.tex
    echo "}" >> $texfile.tex
    echo "\\end{figure}" >> $texfile.tex
    echo "\\clearpage" >> $texfile.tex
done
echo "\\end{document}" >> $texfile.tex
latex -interaction=nonstopmode $texfile.tex > /dev/null 2>&1
latex -interaction=nonstopmode $texfile.tex > /dev/null 2>&1
dvips $texfile.dvi -o $texfile.ps > /dev/null 2>&1
rm -f $texfile.log $texfile.aux $texfile.dvi $texfile.lof $texfile.tex
