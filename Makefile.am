## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

# changed this to your dput command
DPUT = socksify dput

# change this to your launchpad PPA
# configured in ~/.dput.cf
PPA = gerris-ppa

# change this to the ID of the PPA GPG key
KEY = 07AB22DC

# change this to the html changelog for this launchpad release
LOG = /var/www/gerris/version/gerris_changes

SUBDIRS = modules src tools test doc desktop debian

EXTRA_DIST = ChangeLog

if DARCS_CONTROLLED
DARCS_CHANGELOG = changelog
else
DARCS_CHANGELOG = 
endif

ChangeLog: $(DARCS_CHANGELOG)

changelog:
	darcs changes | grep -v Ignore-this: > ChangeLog

debian-snapshot: $(DARCS_CHANGELOG) config.h src/version.h
	sh debian/changelog.sh
	sh debian/control.sh
	dpkg-buildpackage -rfakeroot -b
	rm -f debian/repo/*
	mv -f ../*.changes ../*.deb debian/repo
	cd debian/repo && dpkg-scanpackages . /dev/null | gzip -9c > Packages.gz

launchpad-snapshot: $(DARCS_CHANGELOG) config.h src/version.h dist
	tar xzf gerris-?.?.?.tar.gz
	cd gerris-?.?.? && ./configure && sh debian/changelog.sh $(LOG) && dpkg-buildpackage -rfakeroot -k$(KEY) -S -sa
	rm -f gerris-?.?.?.tar.gz
	rm -r -f gerris-?.?.?
	$(DPUT) $(PPA) gerris-snapshot_*_source.changes
	rm -f gerris-snapshot_*
