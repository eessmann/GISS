import sys
import os
import os.path
import gfs2tex

for start in sys.argv[1:]:
    for root, dirs, files in os.walk(start,topdown=True):
        if not ".xvpics" in root:
            example = gfs2tex.Example(root)
            status,msg = example.test()
            if status != None:
                print root + ": FAIL"
                print " ".join(msg)
            else:
                print root + ": PASS"
