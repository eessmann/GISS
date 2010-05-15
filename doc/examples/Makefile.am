## Process this file with automake to produce Makefile.in

EXAMPLES = \
	cylinder \
	rt \
	boussinesq \
	tangaroa \
	logo \
	tides \
	ship \
	garden \
	dam \
	hump \
	cyclone

EXTRA_DIST = \
	template.tex \
	gfs2tex.py \
	gfs2tex \
	depend.py \
	test.py \
	Makefile.deps \
	crossref.sh.in \
	gfsxref \
	gfs-mode.el

TESTS = test.sh

test.sh: $(EXAMPLES)
	@echo "python test.py $(EXAMPLES)" > test.sh
	@chmod +x test.sh

bin_SCRIPTS = \
	gfs2doc gfs-highlight gfsxref

BUILT_SOURCES= \
	gfs2doc gfs-highlight gfs.lang gerris.dic gfs-keywords.el

CLEANFILES = $(BUILT_SOURCES) Makefile.deps

pkglib_DATA = gfs2tex.py
pkgdata_DATA = gfs.lang gerris.dic gfs-keywords.el gfs-mode.el

gfs.lang: classes
	$(srcdir)/classes

gerris.dic: classes
	$(srcdir)/classes

gfs-keywords.el: classes
	$(srcdir)/classes

gfs2doc: gfs2doc.in

gfs-highlight: gfs-highlight.in

clean-generic:
	$(RM) *.dvi *.aux *.log *.toc *.out examples.tex *.pyc test.sh gfs2doc
	$(RM) -r examples

DOC = examples

examples: examples.dvi crossref.sh
	hevea -fix $(DOC).tex
	imagen -res 600 -extra "pnmscale 0.24" $(DOC)
	hacha $(DOC).html
	rm -f $(DOC).html
	mv -f $(DOC)[0-9][0-9][0-9].png $(DOC)
##	fixme: the character conversion below is a workaround for a bug in hevea version < 1.09
##	for f in *.html; do konwert iso1-utf8 < $$f > $(DOC)/$$f; rm -f $$f; done
	mv -f *.html $(DOC)
	cat $(DOC).css ../share/darcs.css > $(DOC)/$(DOC).css
	sh ../share/fixnav.sh $(DOC)
	cp -f ../share/contents.png ../share/next.png ../share/prev.png $(DOC)
	rm -f *_motif.gif $(DOC).h{tml,aux,ind,toc} $(DOC).image.tex $(DOC).css
	sh ./crossref.sh --url=http://gfs.sourceforge.net/examples/examples $(EXAMPLES)
	mv references examples

examples.dvi: examples.tex
	latex -interaction=nonstopmode examples.tex > /dev/null 2>&1
	latex -interaction=nonstopmode examples.tex > /dev/null 2>&1
	latex -interaction=nonstopmode examples.tex

examples.pdf: examples.dvi
	dvips -Ppdf -G0 examples.dvi -o examples.ps
	ps2pdf -sPAPERSIZE=a4 -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true examples.ps examples.pdf
	rm -f examples.ps

examples.tex: template.tex Makefile.deps gfs2tex gfs2tex.py
	rm -r -f examples
	python gfs2tex $(EXAMPLES)
	sed "s/GFS_VERSION/`$(top_srcdir)/src/gerris2D -V 2>&1 | awk '{ if ($$5 == "version") print $$6}'`/g" < template.tex | sed 's/\\today/'"`date +\"%B %e, %Y\"`/g" > examples.tex

Makefile.deps: Makefile depend.py
	python depend.py $(EXAMPLES) > Makefile.deps

-include Makefile.deps

examples.tar.gz: examples.pdf examples $(DOCS)
	tar czf examples.tar.gz examples $(DOCS)

INCLUDES = -I$(top_srcdir)/src -I$(includedir) -DG_LOG_DOMAIN=\"Gfs-tools\"\
            $(GTS_CFLAGS)

noinst_PROGRAMS = classes

classes: classes.c $(top_srcdir)/src/init.c
	$(LIBTOOL) --mode=link $(CC) $(AM_CFLAGS) $(INCLUDES) -DFTT_2D=1 \
	classes.c -o classes $(GFS2D_LIBS)
