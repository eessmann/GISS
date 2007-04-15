## Process this file with automake to produce Makefile.in

SUBDIRS = modules src tools test doc desktop

EXTRA_DIST = ChangeLog

if DARCS_CONTROLLED
DARCS_CHANGELOG = changelog
else
DARCS_CHANGELOG = 
endif

ChangeLog: $(DARCS_CHANGELOG)

changelog:
	darcs changes > ChangeLog

debian-snapshot: $(DARCS_CHANGELOG)
	sh debian/control.sh
	dpkg-buildpackage -rfakeroot -b -d
	rm -f debian/repo/*
	mv -f ../*.changes ../*.deb debian/repo
	cd debian/repo && dpkg-scanpackages . /dev/null | gzip -9c > Packages.gz
