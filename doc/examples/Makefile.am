## Process this file with automake to produce Makefile.in

EXAMPLES = \
	cylinder \
	shear \
	rt \
	boussinesq \
	tangaroa

EXTRA_DIST = \
	template.tex \
	gfs2tex.py \
	gfs2tex \
	depend.py \
	test.py \
	l2hconf.pm \
	Makefile.deps

TESTS = test.sh

test.sh: $(EXAMPLES)
	@echo "python test.py $(EXAMPLES)" > test.sh
	@chmod +x test.sh

bin_SCRIPTS = \
	gfs2doc

BUILT_SOURCES= \
	gfs2doc

CLEANFILES = $(BUILT_SOURCES) Makefile.deps

dist_pkglib_DATA = gfs2tex.py

gfs2doc: gfs2doc.in

clean-generic:
	$(RM) *.dvi *.aux *.log *.toc *.out examples.tex *.pyc test.sh gfs2doc

examples: examples.ps.gz l2hconf.pm
	latex2html -no_math -html_version 3.2,math -address "" -info "" -split +3 -show_section_numbers -toc_depth 5 -t "Gerris Examples" -local_icons -white examples.tex
	cp -f ../share/darcs.css examples/examples.css

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

examples.tex: template.tex Makefile.deps
	rm -r -f examples
	python gfs2tex $(EXAMPLES)
	sed "s/GFS_VERSION/`$(top_srcdir)/src/gerris2D -V 2>&1 | awk '{ if ($$5 == "version") print $$6}'`/g" < template.tex > examples.tex

Makefile.deps: Makefile depend.py
	python depend.py $(EXAMPLES) > Makefile.deps

-include Makefile.deps

examples.tar.gz: examples.ps.gz examples.pdf examples $(DOCS)
	tar czf examples.tar.gz examples.ps.gz examples.pdf examples $(DOCS)
