#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:00:57 PM}
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

class TestCircleLevelSet(TestLevelSetBase):
    def setUp(self):
        dx = input.dx
        dy = input.dy
        
        self._testInitialEvaluatedIDs = Numeric.array((26, 27, 28,
                                  36, 37, 38, 39, 40,
                                  46, 47, 51, 52,
                                  57, 58, 62, 63,
                                  68, 69, 73, 74,
                                  80, 81, 82, 83, 84,
                                  92, 93, 94))

        dY = dy / 2.
        dX = dx / 2.
        mm = min (dX, dY)
        
        self._testInitialEvaluatedValues = Numeric.array((-1.,  -1.  , -1.  ,  -1.  ,  -1.  ,  -1.  , -1.  ,  -1.  ,  -1.  , -1., -1.,
                                                          -1.,  -1.  , -1.  ,  -1.  ,  -1.  ,  -1.  , -1.  ,  -1.  ,  -1.  , -1., -1.,
                                                          -1.,  -1.  , -1.  ,  -1.  ,  -dY  ,  -dY  , -dY  ,  -1.  ,  -1.  , -1., -1.,
                                                          -1.,  -1.  , -1.  ,  -mm  ,  mm   ,  dY   , mm   ,  -mm  ,  -1.  , -1., -1.,
                                                          -1.,  -1.  , -dX  ,  mm   ,  1.   ,  1    , 1.   ,  mm   ,  -dX  , -1., -1.,
                                                          -1.,  -1.  , -dX  ,  dX   ,  1.   ,  1    , 1.   ,  dX   ,  -dX  , -1., -1.,
                                                          -1.,  -1.  , -dX  ,  mm   ,  1.   ,  1.   , 1.   ,  mm   ,  -dX  , -1., -1.,
                                                          -1.,  -1.  , -1.  ,  -mm  ,  mm   ,  dY   , mm   ,  -mm  ,  -1.  , -1., -1.,
                                                          -1.,  -1.  , -1.  ,  -1.  ,  -dY  ,  -dY  , -dY  ,  -1.  ,  -1.  , -1., -1.,
                                                          -1.,  -1.  , -1.  ,  -1.  ,  -1.  ,  -1.  , -1.  ,  -1.  ,  -1.  , -1., -1.,
                                                          -1.,  -1.  , -1.  ,  -1.  ,  -1.  ,  -1.  , -1.  ,  -1.  ,  -1.  , -1., -1.))

        v1 = self.evalCell(-dY, -mm, dx, dy)[0] 
        v2 = self.evalCell(-mm, -dX, dx, dy)[0]
        v3 = self.evalCell(mm,  mm,  dx, dy)[1]

        self._testInitialTrialValues = Numeric.array((-1.,  -1.  , -1.  ,  -1.  ,  -1.  ,  -1.  , -1.  ,  -1.  ,  -1.  , -1.  , -1.,
                                                      -1.,  -1.  , -1.  ,  -1.  ,  -3*dY,  -3*dY, -3*dY,  -1.  ,  -1.  , -1.  , -1.,
                                                      -1.,  -1.  , -1.  ,  v1   ,  -dY  ,  -dY  , -dY  ,  v1   ,  -1.  , -1.  , -1.,
                                                      -1.,  -1.  , v2   ,  -mm  ,  mm   ,  dY   , mm   ,  -mm  ,  v2   , -1.  , -1.,
                                                      -1.,  -dX*3, -dX  ,  mm   ,  v3   ,  3*dY , v3   ,  mm   ,  -dX  , -dX*3, -1.,
                                                      -1.,  -dX*3, -dX  ,  dX   ,  3*dX ,  1    , 3*dX ,  dX   ,  -dX  , -dX*3, -1.,
                                                      -1.,  -dX*3, -dX  ,  mm   ,  v3   ,  3*dY , v3   ,  mm   ,  -dX  , -dX*3, -1.,
                                                      -1.,  -1.  , v2   ,  -mm  ,  mm   ,  dY   , mm   ,  -mm  ,  v2   , -1.  , -1.,
                                                      -1.,  -1.  , -1.  ,  v1   ,  -dY  ,  -dY  , -dY  ,  v1   ,  -1.  , -1.  , -1.,
                                                      -1.,  -1.  , -1.  ,  -1.  ,  -3*dY,  -3*dY, -3*dY,  -1.  ,  -1.  , -1.  , -1.,
                                                      -1.,  -1.  , -1.  ,  -1.  ,  -1.  ,  -1.  , -1.  ,  -1.  ,  -1.  , -1.  , -1.))

        
        self._testInitialTrialCellIDs = Numeric.array((15, 16, 17,
                                                  25, 29,
                                                  35, 41,
                                                  45, 48, 49, 50, 53,
                                                  56, 59, 61, 64,
                                                  67, 70, 71, 72, 75,
                                                  79, 85,
                                                  91, 95,
                                                  103, 104, 105))

        self.eqn = input.eqn

    def testFinalValues(self):
        pass
        
def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestCircleLevelSet))
    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')
            
            
