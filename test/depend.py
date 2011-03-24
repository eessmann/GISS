#!/usr/bin/python

import sys
import os
import os.path
sys.path.append("../doc/examples")
import gfs2tex

dists = ""
depends = ""
docs = ""
reps = ""
tests = ""
logs = ""
for start in sys.argv[1:]:
    tests += "\\\n\t" + start + ".sh"
    reps += " " + start
    logs += "\\\n\t" + start + ".log"
    f = open(start + ".sh", "w")
    f.write("python -u test.py " + start + "\n")
    f.close()
    os.chmod(start + ".sh",0755)
    for root, dirs, files in os.walk(start,topdown=True):
        if not ".xvpics" in root:
            test = gfs2tex.Example(root)
            name = test.path + "/" + test.name + ".gfs"
            
            docs += "\\\n\t" + name + ".html"
            dists += "\\\n\t" + name
            depends += "\\\n\t" + name
            for f in test.required:
                dists += "\\\n\t" + test.path + "/" + f
            for f in test.generated:
                depends += "\\\n\t" + test.path + "/" + f

print "DOCS = " + docs + dists
print ""
print "EXTRA_DIST += " + dists
print ""
print "TESTS = " + tests + "\\\n\tsummary.sh"
os.chmod("summary.sh",0755)
print ""
print "TESTS_ENVIRONMENT = TESTS=\"" + reps + "\""
print "TEST_EXTENSIONS = .sh"
print "summary.log:" + logs
print ""
print "tests.tex: " + depends
