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
import input
import Numeric
import MA

class TestInterior(TestLevelSetBase):
    def setUp(self):
        dx = input.dx
        dy = input.dy
        self._testInitialEvaluatedIDs = Numeric.array((1, 2, 3,
                                                       5, 6, 7, 8, 9,
                                                       10, 11, 13, 14,
                                                       15, 16, 17, 18, 19,
                                                       21, 22, 23))
        dX = dx / 2.
        dY = dy / 2.
        mm = min(dX, dY)
        
        self._testInitialEvaluatedValues = Numeric.array((  1.  ,   dY  ,   dY  ,  dY  ,  1.  ,
                                                            dX  ,   -mm ,   -dY ,  -mm ,  dX  ,
                                                            dX  ,   -dX ,   -1. ,  -dX ,  dX  ,
                                                            dX  ,   -mm ,   -dY ,  -mm ,  dX  ,
                                                            1.  ,   dY  ,   dY  ,  dY  ,  1.  ))

                                                            

        self._testInitialTrialCellIDs = Numeric.array((0, 4, 12, 20, 24))
                                                            
        v1 = self.evalCell(dY, dX, dx, dy)[1]
        v2 = max(-dY*3, -dX*3)   
                                                            
        self._testInitialTrialValues = Numeric.array((  v1  ,   dY  ,   dY  ,  dY  ,  v1  ,
                                                        dX  ,   -mm ,   -dY ,  -mm ,  dX  ,
                                                        dX  ,   -dX ,   -v1 ,  -dX ,  dX  ,
                                                        dX  ,   -mm ,   -dY ,  -mm ,  dX  ,
                                                        v1  ,   dY  ,   dY  ,  dY  ,  v1  ))

        self._testFinalValues = self._testInitialTrialValues

        self.eqn = input.eqn

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestInterior))
    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')

            
            
