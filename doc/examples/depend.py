#!/usr/bin/python

import sys
import os
import os.path
import gfs2tex

dists = ""
depends = ""
docs = ""

for start in sys.argv[1:]:
    for root, dirs, files in os.walk(start,topdown=True):
        if not ".xvpics" in root:
            example = gfs2tex.Example(root)
            name = example.path + "/" + example.name + ".gfs"
            docs += "\\\n\t" + name + ".html"
            dists += "\\\n\t" + name
            depends += "\\\n\t" + name
            for f in example.required:
                if os.path.exists("../../_darcs/pristine/doc/examples/" + example.path + "/" + f):
                    dists += "\\\n\t" + example.path + "/" + f
            for f in example.generated:
                depends += "\\\n\t" + example.path + "/" + f
                if f[-4:] == ".mpg":
                    docs += "\\\n\t" + example.path + "/" + f

print "DOCS = " + docs + dists
print ""
print "EXTRA_DIST += " + dists
print ""
print "examples.tex: " + depends
