## Process this file with automake to produce Makefile.in

TESTS = test.sh

EXTRA_DIST = \
	check.sh \
	cleanup.sh \
	figures.sh \
	report.sh \
	tests

test.sh:
	@echo "#! /bin/sh" > test.sh
	@echo "rm -f report.ps" >> test.sh
	@echo "for f in tests/*; do if test \"\$$f\" != \"tests/CVS\"; then ./check.sh \"\$$f\"; fi; done" >> test.sh
	@chmod +x test.sh

report.ps: report.sh timestamp
	./report.sh "tests/*"

clean-generic:
	./cleanup.sh
