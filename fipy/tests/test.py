#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/26/03 {3:23:47 PM}
 #                                last update: 11/26/03 {11:11:48 AM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

"""Run all the test cases
"""

import tests.testVariableDiffusion
import tests.testSteadyStateDiffusion
import tests.testExplicitDiffusion
import unittest

class AllTestSuites(unittest.TestSuite):
    """
    Class for adding all the test suites into one suite.
    """
    def __init__(self, suites=()):
        unittest.TestSuite.__init__(self)
        for suite in suites:
            self.addTests(suite._tests)
        
def suite1():
    alltests = AllTestSuites((
        tests.testSteadyStateDiffusion.suite(),
        tests.testVariableDiffusion.suite(),
        tests.testExplicitDiffusion.suite()
        ))
    return alltests

if __name__ == '__main__':
    unittest.main(defaultTest='suite1')

