## Process this file with automake to produce Makefile.in

TESTS = test.sh

EXTRA_DIST = \
	check.sh \
	cleanup.sh \
	report.sh \
	tests

test.sh: star.gts
	@echo "#! /bin/sh" > test.sh
	@echo "rm -f report.ps" >> test.sh
	@echo "./check.sh \"tests/*\";" >> test.sh
	@chmod +x test.sh

report.ps: report.sh timestamp
	./report.sh "tests/*"

star.gts:
	../../poisson/shapes star | transform -i -s 0.9 > star.gts

clean-generic:
	./cleanup.sh
