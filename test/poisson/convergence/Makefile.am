## Process this file with automake to produce Makefile.in

EXTRA_DIST = \
	divergence.awk \
	divergence_rate.awk \
	divfig.sh \
	report.sh \
	runcomp.sh \
	runfig.sh \
	reference2D \
	reference3D \
	divfig.par \
	check.sh \
	checkrate.awk \
	avgrate.awk \
	makeref \
	makeref.sh

TESTS = test.sh

REFERENCE = reference2D

test.sh:
	@echo "#! /bin/sh" > test.sh
	@echo "rm -f report.ps" >> test.sh
	@echo "./check.sh \"$(REFERENCE)/*.xmgr\" 0.9" >> test.sh
	@chmod +x test.sh

report.ps: ../poisson report.sh timestamp
	./report.sh "$(REFERENCE)/*.xmgr"

reference:
	./makeref.sh $(REFERENCE)

clean-generic:
	$(RM) -r -f figures test
