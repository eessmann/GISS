## Process this file with automake to produce Makefile.in

TESTS = test.sh

EXTRA_DIST = \
	order.sh \
	check.sh \
	order.par \
	reference \
	avgorder.awk \
	checkorder.awk \
	report.sh

test.sh:
	@echo "#! /bin/sh" > test.sh
	@echo "rm -f report.ps" >> test.sh
	@echo "./check.sh \"reference/*.xmgr\" 0.9" >> test.sh
	@chmod +x test.sh

report.ps: report.sh timestamp
	./report.sh "reference/*.xmgr"

clean-generic:
	$(RM) -r -f figures test
