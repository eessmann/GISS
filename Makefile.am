## Process this file with automake to produce Makefile.in

SUBDIRS = modules src tools test doc

EXTRA_DIST = ChangeLog

dist-hook:
	$(RM) -r `find $(distdir) -name CVS -type d -print`
	$(RM) `find $(distdir) -name ".cvsignore" -type f -print`

ChangeLog:
	darcs changes > ChangeLog

debian-snapshot:
	sh debian/changelog.sh
	sh debian/control.sh
	dpkg-buildpackage -rfakeroot -b -d
	rm -f debian/repo/*
	mv -f ../*.changes ../*.deb debian/repo
	cd debian/repo && dpkg-scanpackages . /dev/null | gzip -9c > Packages.gz
