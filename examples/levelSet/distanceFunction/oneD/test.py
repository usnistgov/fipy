#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:05:56 PM}
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
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

"""Test steady-state diffusion solutions
"""
 
import unittest
from examples.levelSet.testLevelSetBase import TestLevelSetBase
import fipy.tests.testProgram
import Numeric
import MA
import input

class Test1D(TestLevelSetBase):
    def setUp(self):
        dx = input.dx
        dy = input.dy
        self._testInitialEvaluatedIDs = Numeric.array((4, 5))

        self._testInitialEvaluatedValues = Numeric.array((1., 1., 1., 1., dx / 2., -dx / 2., -1., -1., -1., -1.))

        self._testInitialTrialCellIDs = Numeric.array((3, 6))

        self._testInitialTrialValues = Numeric.array((1., 1., 1., 3. * dx / 2., dx / 2., -dx / 2., -3. * dx / 2., -1., -1., -1.))

        self._testFinalValues = Numeric.array((9. * dx / 2., 7. * dx / 2., 5. * dx / 2., 3. * dx / 2., dx / 2.,
                                               -dx / 2., -3. * dx / 2., -5. * dx / 2., -7. * dx / 2., -9. * dx / 2.))

        self.eqn = input.eqn

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(Test1D))
    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')
