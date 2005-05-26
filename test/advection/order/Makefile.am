## Process this file with automake to produce Makefile.in

TESTS = test.sh

EXTRA_DIST = \
	order.sh \
	check.sh \
	order.par \
	report.sh \
	tests

test.sh: Makefile.am
	@echo "#! /bin/sh" > test.sh
	@echo "rm -f report.ps" >> test.sh
	@echo "sh check.sh \"tests/*\" 1.9" >> test.sh
	@chmod +x test.sh

report.ps: report.sh timestamp
	sh report.sh "tests/*"

clean-generic:
	$(RM) -f tests/*/result.* tests/*/*~ tests/*/core timestamp test.sh
