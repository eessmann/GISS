## Process this file with automake to produce Makefile.in

EXTRA_DIST = \
	order.awk \
	order1.awk \
	report.sh \
	orderfig.sh \
	runcomp.sh \
	runorderfig.sh \
	reference2D \
	reference3D \
	orderfig.par \
	check.sh \
	avgorder.awk \
	checkorder.awk

TESTS = test.sh

REFERENCE = reference2D

test.sh:
	@echo "#! /bin/sh" > test.sh
	@echo "rm -f report.ps" >> test.sh
	@echo "./check.sh \"$(REFERENCE)/*.xmgr\" 0.9" >> test.sh
	@chmod +x test.sh

report.ps: ../poisson report.sh timestamp
	./report.sh "$(REFERENCE)/*.xmgr"

clean-generic:
	$(RM) -r -f figures test
