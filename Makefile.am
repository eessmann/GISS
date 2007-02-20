## Process this file with automake to produce Makefile.in

SUBDIRS = modules src tools test doc

EXTRA_DIST = ChangeLog

ChangeLog:
	darcs changes > ChangeLog

debian-snapshot:
	sh debian/control.sh
	dpkg-buildpackage -rfakeroot -b -d
	rm -f debian/repo/*
	mv -f ../*.changes ../*.deb debian/repo
	cd debian/repo && dpkg-scanpackages . /dev/null | gzip -9c > Packages.gz
