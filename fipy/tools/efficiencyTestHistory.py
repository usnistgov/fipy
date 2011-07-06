##This script runs on using "file", which is created with the command $svn log --xml --incremental > file
##Update file if there is a new revision to /trunk/examples. Otherwise you may ignore this.

from datetime import datetime
import os
import pysvn
import subprocess

def run(startRev):
    datafile = open("file", 'r')
    datafilelist = datafile.readlines()
    revisions = [int(line.lstrip('   revision="').rstrip('">\n')) for line in datafilelist if 'revision="' in line]
    revisions.sort()
    while(startRev not in revisions):
        startRev += 1 
    index = revisions.index(startRev)
    os.chdir("../../../trunk/examples")
    for k in revisions[index:]:
        revisionNumber = pysvn.Client().update(".", revision=pysvn.Revision(pysvn.opt_revision_kind.number, k))
        print "pysvn.Client().info('.')['revision'].number: ", pysvn.Client().info('.')['revision'].number
        os.chdir("../../efficiency_test")
        print 'hello'
        w, r = os.popen2("python setup.py efficiency_test --uploadToCodespeed")
        os.wait()
        os.chdir("../trunk/examples")

if __name__ == '__main__':
    run(4300)



