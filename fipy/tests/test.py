#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/26/03 {3:23:47 PM}
 #                                last update: 11/29/03 {9:08:49 PM} 
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
import tests.testPhase
import unittest

if __name__ == '__main__':
    theSuite = tests.testSteadyStateDiffusion.suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)
    
    theSuite = tests.testExplicitDiffusion.suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)
    
    theSuite = tests.testVariableDiffusion.suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)

    theSuite = tests.testPhase.suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)
    
