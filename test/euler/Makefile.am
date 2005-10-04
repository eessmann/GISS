## Process this file with automake to produce Makefile.in

TESTDIRS = \
	poisson \
	circle \
	advection \
	reynolds \
	periodic \
	merging \
	boundaries \
	channel \
	spurious \
	capwave \
	waves \
	geo

EXTRA_DIST = \
	template.tex \
	tests.bib \
	depend.py \
	test.py \
	gfs2tex \
	l2hconf.pm \
	Makefile.deps

TESTS = test.sh

test.sh: $(TESTDIRS)
	@echo "python test.py $(TESTDIRS)" > test.sh
	@chmod +x test.sh

clean-generic:
	$(RM) *.dvi *.aux *.log *.toc *.out tests.tex *.pyc test.sh *.bbl *.blg

tests: tests.ps.gz l2hconf.pm
	latex2html -no_math -html_version 3.2,math -address "" -info "" -split +3 -show_section_numbers -toc_depth 5 -t "Gerris Tests" -local_icons -white tests.tex
	cp -f ../../doc/share/darcs.css tests/tests.css

tests.dvi: tests.tex tests.bib
	latex -interaction=nonstopmode tests.tex > /dev/null 2>&1
	bibtex tests
	latex -interaction=nonstopmode tests.tex > /dev/null 2>&1
	latex -interaction=nonstopmode tests.tex

tests.ps.gz: tests.dvi
	dvips tests.dvi -o tests.ps
	gzip -f --best tests.ps

tests.pdf: tests.dvi
	dvips -Ppdf -G0 tests.dvi -o tests.ps
	ps2pdf -sPAPERSIZE=a4 -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true tests.ps tests.pdf
	rm -f tests.ps

tests.tex: template.tex Makefile.deps
	rm -r -f tests
	python gfs2tex $(TESTDIRS)
	cp -f template.tex tests.tex

Makefile.deps: Makefile.am depend.py
	python depend.py $(TESTDIRS) > Makefile.deps

-include Makefile.deps

tests.tar.gz: tests $(DOCS)
	tar czf tests.tar.gz tests $(DOCS)
