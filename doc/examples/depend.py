#!/usr/bin/python

import sys
import os
import os.path
import gfs2tex

dists = ""
depends = ""
for start in sys.argv[1:]:
    for root, dirs, files in os.walk(start,topdown=True):
        if not ".xvpics" in root:
            example = gfs2tex.Example(root)
            dists += "\\\n\t" + example.path + "/" + example.name + ".gfs "
            depends += "\\\n\t" + example.path + "/" + example.name + ".gfs "
            for f in example.required:
                dists += "\\\n\t" + example.path + "/" + f + " "
            for f in example.generated:
                depends += "\\\n\t" + example.path + "/" + f + " "
print "EXTRA_DIST += " + dists
print ""
print "examples.tex: " + depends
