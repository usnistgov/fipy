##This script runs on using "file", which is created with the command $svn log --xml --incremental > file
##Update file if there is a new revision to /trunk/examples. Otherwise you may ignore this.

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "efficiencyTestHistory.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: Andrew Acquaviva <andrewa@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #
 # ###################################################################
 ##

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
