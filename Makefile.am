## Process this file with automake to produce Makefile.in

SUBDIRS = modules src tools test doc

dist-hook:
	$(RM) -r `find $(distdir) -name CVS -type d -print`
	$(RM) `find $(distdir) -name ".cvsignore" -type f -print`

changelog:
	darcs changes > ChangeLog
