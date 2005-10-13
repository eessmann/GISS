## Process this file with automake to produce Makefile.in

EXTRA_DIST=tutorial.ps.gz tutorial.pdf tutorial.tar.gz half-cylinder.gts \
	gfs2D.net

Makefile.deps: Makefile.am tutorial.tex
	awk -f depend.awk file="tutorial" < tutorial.tex > Makefile.deps

clean-generic:
	$(RM) *.dvi *.aux *.log *.toc *.out \
	boundaries.pdf direction.eps direction.pdf dxscreen.eps dxscreen.pdf \
	gfs2oogl.pdf half-cylinder.eps half-cylinder.pdf refined1.pdf refined1_cells.pdf \
	refined2.pdf vorticity.eps vorticity.pdf gfsview.eps tutorial1.tex \
	Makefile.deps

tutorial.tar.gz: tutorial.ps.gz
	rm -r -f tutorial
	sed 's/input{pdf.tex}/usepackage{graphicx}\\newcommand{\\gfx}{eps}/g' < tutorial.tex | sed "s/GFS_VERSION/`$(top_srcdir)/src/gerris2D -V 2>&1 | awk '{ if ($$5 == "version") print $$6}'`/g" > tutorial1.tex
	latex2html -no_math -html_version 3.2,math -address "" -info "" -split +2 -show_section_numbers -toc_depth 5 -t "The Gerris Tutorial" -local_icons -white tutorial1.tex
	mv -f tutorial1 tutorial
	sh pre_fix.sh
	cp -f ../share/darcs.css tutorial/tutorial1.css
	tar cf tutorial.tar tutorial
	gzip -f --best tutorial.tar

tutorial1.dvi: tutorial.tex Makefile.deps
	sed "s/GFS_VERSION/`$(top_srcdir)/src/gerris2D -V 2>&1 | awk '{ if ($$5 == "version") print $$6}'`/g" < tutorial.tex > tutorial1.tex
	latex -interaction=nonstopmode tutorial1.tex > /dev/null 2>&1
	latex -interaction=nonstopmode tutorial1.tex > /dev/null 2>&1
	latex -interaction=nonstopmode tutorial1.tex

tutorial.ps.gz: tutorial1.dvi
	dvips tutorial1.dvi -o tutorial.ps
	gzip -f --best tutorial.ps

tutorial.pdf: tutorial1.dvi
	dvips -Ppdf -G0 tutorial1.dvi -o tutorial1.ps
	ps2pdf -sPAPERSIZE=a4 -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true tutorial1.ps tutorial.pdf
	rm -f tutorial1.ps

%.eps: %.xmgr
	/usr/local/grace/bin/xmgrace -noask -hardcopy -hdevice EPS -printfile $@ $<

%.eps: %.agr
	/usr/local/grace/bin/xmgrace -noask -hardcopy -hdevice EPS -printfile $@ $<

%.eps: %.fig
	fig2dev -L ps $< > $@ ; sh epsbbox.sh $@ 

%.eps: %.jpeg
	convert $< $@.gif
	convert $@.gif EPS:$@
	- rm -f $@.gif

%.eps: %.jpg
	convert $< $@.gif
	convert $@.gif EPS:$@
	- rm -f $@.gif

%.eps: %.gif
	convert $< EPS:$@

%.eps: %.png
	convert $< EPS:$@

%.pdf: %.jpg
	convert $< EPDF:$@

%.pdf: %.png
	convert $< EPDF:$@

%.pdf: %.epsf
	epstopdf $<

%.pdf: %.eps
	epstopdf $<

-include Makefile.deps
