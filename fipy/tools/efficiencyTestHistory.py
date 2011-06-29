

from datetime import datetime
import os
import pysvn
import subprocess

datafile = open("file", 'r')
datafilelist = datafile.readlines()
revisions = [int(line.lstrip('   revision="').rstrip('">\n')) for line in datafilelist if 'revision="' in line]
revisions.sort()
os.chdir("./trunk/examples")
##pysvn.Client().update(".", revision=pysvn.Revision(pysvn.opt_revision_kind.number, 4581))


for i in revisions[len(revisions)-10:]:
    revisionNumber = pysvn.Client().update(".", revision=pysvn.Revision(pysvn.opt_revision_kind.number, i))
    print "pysvn.Client().info('.')['revision'].number: ", pysvn.Client().info('.')['revision'].number
    os.chdir("../../efficiency_test")
##    print 'hello'
    w, r = os.popen2("python setup.py efficiency_test --uploadToCodespeed --newElements=1000")
    os.wait()
##    print r.readlines()
    os.chdir("../trunk/examples")



