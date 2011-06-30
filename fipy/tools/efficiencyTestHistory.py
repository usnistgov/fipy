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
    start = 0
    i = 0
    while(start < startRev):
        print 'start', start
        print 'startRev:' , startRev
        start = revisions[i]
        i += 1
    os.chdir("../trunk/examples")
    for k in revisions[i:]:
        revisionNumber = pysvn.Client().update(".", revision=pysvn.Revision(pysvn.opt_revision_kind.number, k))
        print "pysvn.Client().info('.')['revision'].number: ", pysvn.Client().info('.')['revision'].number
        os.chdir("../../efficiency_test")
        w, r = os.popen2("python setup.py efficiency_test")
        print r.readlines()
        os.wait()
        os.chdir("../trunk/examples")

if __name__ == '__main__':
    run(4150)



