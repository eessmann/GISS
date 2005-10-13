## Process this file with automake to produce Makefile.in

SUBDIRS = tutorial examples faq

# The name of the module.
DOC_MODULE=gfs

# The top-level SGML file.
DOC_MAIN_SGML_FILE=gfs-docs.sgml

# The directory containing the source code (if it contains documentation).
DOC_SOURCE_DIR=../src

TARGET_DIR=html/$(DOC_MODULE)

EXTRA_DIST=html figures

scan:
	gtkdoc-scan --module=$(DOC_MODULE) --source-dir=$(DOC_SOURCE_DIR) --ignore-headers=""

templates: scan
	gtkdoc-mktmpl --module=$(DOC_MODULE)
	sort -f $(DOC_MODULE)-unused.txt > /tmp/toto; mv -f /tmp/toto $(DOC_MODULE)-unused.txt

sgml: templates
	gtkdoc-mkdb --module=$(DOC_MODULE) --source-dir=$(DOC_SOURCE_DIR)
	sort -f $(DOC_MODULE)-undocumented.txt > /tmp/toto; mv -f /tmp/toto $(DOC_MODULE)-undocumented.txt

html: sgml
	if ! test -d html ; then mkdir html ; fi
	-cd html && gtkdoc-mkhtml $(DOC_MODULE) ../$(DOC_MAIN_SGML_FILE)

clean-local:
	rm -f *~ *.bak *.hierarchy *.signals *.args *-unused.txt

maintainer-clean-local: clean
	rm -rf sgml html $(DOC_MODULE)-decl-list.txt $(DOC_MODULE)-decl.txt 

## install-data-local:
##	install -d -m 0755 $(TARGET_DIR)
##	install -m 0644 html/*.html $(TARGET_DIR)
##	install -m 0644 html/index.sgml $(TARGET_DIR)
##	gtkdoc-fixxref --module=$(DOC_MODULE) --html-dir=$(HTML_DIR)

.PHONY : html sgml templates scan
