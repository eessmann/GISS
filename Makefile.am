## Process this file with automake to produce Makefile.in

SUBDIRS = modules src tools test doc

EXTRA_DIST = ChangeLog

dist-hook:
	$(RM) -r `find $(distdir) -name CVS -type d -print`
	$(RM) `find $(distdir) -name ".cvsignore" -type f -print`

ChangeLog:
	darcs changes > ChangeLog
