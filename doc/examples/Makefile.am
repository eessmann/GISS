## Process this file with automake to produce Makefile.in

EXAMPLES = \
	cylinder

EXTRA_DIST = \
	$(EXAMPLES) \
	template.tex \
	gfs2tex.py \
	gfs2tex

bin_SCRIPTS = \
	gfs2doc

BUILT_SOURCES= \
	gfs2doc

CLEANFILES = $(BUILT_SOURCES)

dist_pkglib_DATA = gfs2tex.py

gfs2doc: gfs2doc.in

clean-generic:
	$(RM) *.dvi *.aux *.log *.toc *.out examples.tex *.pyc

examples: examples.ps.gz
	latex2html -no_math -html_version 3.2,math -address "" -info "" -split +3 -show_section_numbers -toc_depth 5 -t "Gerris Examples" -local_icons -white examples.tex

examples.dvi: examples.tex
	latex -interaction=nonstopmode examples.tex > /dev/null 2>&1
	latex -interaction=nonstopmode examples.tex > /dev/null 2>&1
	latex -interaction=nonstopmode examples.tex

examples.ps.gz: examples.dvi
	dvips examples.dvi -o examples.ps
	gzip -f --best examples.ps

examples.pdf: examples.dvi
	dvips -Ppdf -G0 examples.dvi -o examples.ps
	ps2pdf -sPAPERSIZE=a4 -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true examples.ps examples.pdf
	rm -f examples.ps

examples.tex: template.tex
	rm -r -f examples
	./gfs2tex
	sed "s/GFS_VERSION/`$(top_srcdir)/src/gerris2D -V 2>&1 | awk '{ if ($$5 == "version") print $$6}'`/g" < template.tex > examples.tex
