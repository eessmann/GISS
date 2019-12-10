## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

# change this to your OSC directory
OSC_DIR=$(HOME)/local/src/osc/home:popinet/$(PACKAGE)-snapshot

# change this to your dput command
DPUT = socksify dput

# change this to your launchpad PPA
# configured in ~/.dput.cf
PPA = gerris-ppa

# change this to the ID of the PPA GPG key
KEY = 07AB22DC

# change this to the html changelog for this launchpad release
LOG = /var/www/gerris/version/gfsview_changes

if HAVE_OSMESA
  BATCH = batch
endif

if HAVE_GTK
  INTERACTIVE = view
endif

SUBDIRS = gl2ps gl $(BATCH) $(INTERACTIVE) debian

EXTRA_DIST = \
	autogen.sh \
	gl2ps/COPYING.GL2PS \
	gl2ps/COPYING.LGPL \
	m4

if DARCS_CONTROLLED
DARCS_CHANGELOG = changelog
else
DARCS_CHANGELOG = 
endif

install-data-local:
	@$(NORMAL_INSTALL)
	if test -d $(srcdir)/pixmaps; then \
	  $(mkinstalldirs) $(DESTDIR)$(pkgdatadir)/pixmaps; \
	  for pixmap in $(srcdir)/pixmaps/*; do \
	    if test -f $$pixmap; then \
	      $(INSTALL_DATA) $$pixmap $(DESTDIR)$(pkgdatadir)/pixmaps; \
	    fi \
	  done \
	fi; \
	if test -d $(srcdir)/fonts; then \
	  $(mkinstalldirs) $(DESTDIR)$(pkgdatadir)/fonts; \
	  for font in $(srcdir)/fonts/*; do \
	    if test -f $$font; then \
	      $(INSTALL_DATA) $$font $(DESTDIR)$(pkgdatadir)/fonts; \
	    fi \
	  done \
	fi

dist-hook:
	if test -d pixmaps; then \
	  mkdir $(distdir)/pixmaps; \
	  for pixmap in pixmaps/*; do \
	    if test -f $$pixmap; then \
	      cp -p $$pixmap $(distdir)/pixmaps; \
	    fi \
	  done \
	fi; \
	if test -d fonts; then \
	  mkdir $(distdir)/fonts; \
	  for font in fonts/*; do \
	    if test -f $$font; then \
	      cp -p $$font $(distdir)/fonts; \
	    fi \
	  done \
	fi

ChangeLog: $(DARCS_CHANGELOG)

changelog:
	darcs changes > ChangeLog

debian-snapshot: $(DARCS_CHANGELOG) config.h
	sh debian/changelog.sh
	dpkg-buildpackage -rfakeroot -b
	rm -f debian/repo/*
	mv -f ../*.changes ../*.deb debian/repo
	cd debian/repo && dpkg-scanpackages . /dev/null | gzip -9c > Packages.gz

launchpad-snapshot: $(DARCS_CHANGELOG) config.h dist
	tar xzf gfsview-?.?.?.tar.gz
	cd gfsview-?.?.? && ./configure && sh debian/changelog.sh $(LOG) && chmod +x debian/rules && dpkg-buildpackage -rfakeroot -k$(KEY) -S -sa
	rm -f gfsview-?.?.?.tar.gz
	rm -r -f gfsview-?.?.?
	$(DPUT) $(PPA) gfsview-snapshot_*_source.changes
	rm -f gfsview-snapshot_*

buildservice-snapshot: $(DARCS_CHANGELOG) config.h dist
	tar xzf $(PACKAGE)-?.?.?.tar.gz
	rm -f $(PACKAGE)-?.?.?.tar.gz
	mv -f $(PACKAGE)-?.?.? $(PACKAGE)-snapshot
	cd $(PACKAGE)-snapshot && ./configure && sh debian/changelog.sh $(LOG) && chmod +x debian/rules && make distclean
	tar czvf $(PACKAGE)-snapshot.tar.gz $(PACKAGE)-snapshot
	dpkg-source -b $(PACKAGE)-snapshot
	rm -r -f $(PACKAGE)-snapshot
	cd $(OSC_DIR) && osc update && rm -f $(OSC_DIR)/*
	cp -f rpm/$(PACKAGE).spec $(OSC_DIR)
	mv -f $(PACKAGE)-snapshot.tar.gz $(PACKAGE)-snapshot_*.tar.gz $(PACKAGE)-snapshot_*.dsc $(OSC_DIR)
	cd $(OSC_DIR) && osc addremove && osc commit -m "snapshot release"
