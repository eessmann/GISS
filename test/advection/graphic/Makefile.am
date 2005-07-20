## Process this file with automake to produce Makefile.in

TESTS = test.sh

EXTRA_DIST = \
	check.sh \
	conservation.sh \
	report.sh \
	tests

test.sh: star.gts
	@echo "#! /bin/sh" > test.sh
	@echo "rm -f report.ps" >> test.sh
	@echo "sh check.sh \"tests/*\";" >> test.sh
	@chmod +x test.sh

report.ps: report.sh timestamp
	sh report.sh "tests/*"

star.gts:
	shapes star | transform -i -s 0.9 > star.gts

clean-generic:
	$(RM) -f tests/*/*.eps tests/*/log tests/*/*~ tests/*/core timestamp test.sh star.gts
