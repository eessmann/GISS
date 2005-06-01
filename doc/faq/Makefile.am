## Process this file with automake to produce Makefile.in

EXTRA_DIST=faq.ps.gz faq.pdf faq.tar.gz

clean-generic:
	$(RM) *.dvi *.aux *.log *.toc *.out

faq.tar.gz: faq.ps.gz
	rm -r -f faq
	latex2html -no_math -html_version 3.2,math -address "" -info "" -split +2 -toc_depth 5 -t "The Gerris FAQ" -local_icons -white faq.tex
	sh pre_fix.sh
	cp -f ../share/darcs.css faq/faq.css
	tar cf faq.tar faq
	gzip -f --best faq.tar

faq.dvi: faq.tex
	latex -interaction=nonstopmode faq.tex > /dev/null 2>&1
	latex -interaction=nonstopmode faq.tex > /dev/null 2>&1
	latex -interaction=nonstopmode faq.tex

faq.ps.gz: faq.dvi
	dvips faq.dvi -o faq.ps
	gzip -f --best faq.ps

faq.pdf: faq.dvi
	dvips -Ppdf -G0 faq.dvi -o faq.ps
	ps2pdf -sPAPERSIZE=a4 -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true faq.ps faq.pdf
	rm -f faq.ps
