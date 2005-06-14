import sys
import os
import os.path
sys.path.append("../../doc/examples")
import gfs2tex

env = ""
for p in ["gerris2D","gerris3D","gerris2D3"]:
    env += p + "=" + os.getcwd() + "/../../src/" + p + " "
env += " PYTHONPATH=$PYTHONPATH:" + os.getcwd() + " "

n = 0
failed = 0
for start in sys.argv[1:]:
    for root, dirs, files in os.walk(start,topdown=True):
        if not ".xvpics" in root:
            test = gfs2tex.Example(root)
            status,msg = test.run(env)
            if status != None:
                print "FAIL:",root
                if len(msg) > 0:
                    print " ".join(msg)
                print >>open(test.path + "/status",'w'), "{\color{red}FAIL}:"
                failed += 1
            else:
                print "PASS:",root
                print >>open(test.path + "/status",'w'), "{\color{green}PASS}:"
            n += 1
if failed:
    msg = repr(failed) + " of " + repr(n) + " tests failed"
else:
    msg = "All " + repr(n) + " tests passed"

print len(msg)*"="
print msg
print len(msg)*"="

if failed:
    sys.exit(1)
