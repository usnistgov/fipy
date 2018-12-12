##This script runs on using "file", which is created with the command $svn log --xml --incremental > file
##Update file if there is a new revision to /trunk/examples. Otherwise you may ignore this.


__docformat__ = 'restructuredtext'

__all__ = []

from datetime import datetime
import os
import pysvn
import subprocess
from fipy.tools.performance.efficiency_test import Efficiency_test
from setuptools import setup

def run(startRev):
    dummyCommand = setup(name='dummy', script_name = 'setup.py', script_args = ['test', '--dry-run'])
    test = Efficiency_test(dummyCommand)
    test.initialize_options()
    test.uploadToCodespeed = True
    test.newElements = 100
    datafile = open("file", 'r')
    datafilelist = datafile.readlines()
    revisions = [int(line.lstrip('   revision="').rstrip('">\n')) for line in datafilelist if 'revision="' in line]
    revisions.sort()
    while(startRev not in revisions):
        startRev += 1
    index = revisions.index(startRev)
    os.chdir("../../../trunk/examples")
    for k in revisions[index:3511]:
        revisionNumber = pysvn.Client().update(".", revision=pysvn.Revision(pysvn.opt_revision_kind.number, k))
        print "pysvn.Client().info('.')['revision'].number: ", pysvn.Client().info('.')['revision'].number
        os.chdir("../../efficiency_test")
        print 'hello'
        test.run()
        os.wait()
        os.chdir("../trunk/examples")

if __name__ == '__main__':
    run(3000)
