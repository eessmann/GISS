#!/usr/bin/python

import sys
import os
import os.path
import gfs2tex

print "EXTRA_DIST +=",
for start in sys.argv[1:]:
    for root, dirs, files in os.walk(start,topdown=True):
        if not ".xvpics" in root:
            example = gfs2tex.Example(root)
            print example.path + "/" + example.name + ".gfs",
            for f in example.required:
                print example.path + "/" + f,
print ""
