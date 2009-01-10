## Process this file with automake to produce Makefile.in

TESTDIRS = \
	poisson \
	circle \
	dumbell \
	advection \
	shear \
	reynolds \
	periodic \
	merging \
	axi \
	axiadvection \
	lid \
	poiseuille \
	couette \
	kinetic \
	hydrostatic \
	boundaries \
	channel \
	plate \
	spurious \
	capwave \
	oscillation \
	geo \
	waves \
	nz \
	parabola

EXTRA_DIST = \
	template.tex \
	tests.bib \
	depend.py \
	test.py \
	check.py \
	gfs2tex \
	Makefile.deps

TESTS = test.sh

test.sh: $(TESTDIRS)
	@echo "python -u test.py $(TESTDIRS)" > test.sh
	@chmod +x test.sh

clean-generic:
	$(RM) *.dvi *.aux *.log *.toc *.out tests.tex *.pyc test.sh *.bbl *.blg Makefile.deps

DOC = tests

tests: tests.dvi
	hevea -fix $(DOC).tex
	imagen -res 500 -extra "pnmscale 0.24" $(DOC)
	hacha $(DOC).html
	rm -f $(DOC).html
	mv -f $(DOC)[0-9][0-9][0-9].png $(DOC)
##	fixme: the character conversion below is a workaround for a bug in hevea version < 1.09
	for f in *.html; do konwert iso1-utf8 < $$f > $(DOC)/$$f; rm -f $$f; done
	cat $(DOC).css ../doc/share/darcs.css > $(DOC)/$(DOC).css
	sh ../doc/share/fixnav.sh $(DOC)
	cp -f ../doc/share/contents.png ../doc/share/next.png ../doc/share/prev.png $(DOC)
	rm -f *_motif.gif $(DOC).h{tml,aux,ind,toc} $(DOC).image.tex $(DOC).css
	sh ../doc/examples/crossref.sh --url=http://gfs.sourceforge.net/tests/tests $(TESTDIRS)
	mv references tests

tests.dvi: tests.tex tests.bib
	latex -interaction=nonstopmode tests.tex > /dev/null 2>&1
	bibtex tests
	latex -interaction=nonstopmode tests.tex > /dev/null 2>&1
	latex -interaction=nonstopmode tests.tex

tests.pdf: tests.dvi
	dvips -Ppdf -G0 tests.dvi -o tests.ps
	ps2pdf -sPAPERSIZE=a4 -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true tests.ps tests.pdf
	rm -f tests.ps

tests.tex: template.tex Makefile.deps
	rm -r -f tests
	python gfs2tex $(TESTDIRS)
	cp -f template.tex tests.tex

Makefile.deps: Makefile depend.py
	python depend.py $(TESTDIRS) > Makefile.deps

-include Makefile.deps

tests.tar.gz: tests $(DOCS)
	tar czf tests.tar.gz tests $(DOCS)
