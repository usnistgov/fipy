#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 2/13/04 {1:51:48 PM} 
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
import fivol.tests.testProgram
from fivol.tests.testBase import TestBase
from fivol.examples.phase.examples.symmetry.input import SymmetrySystem

class TestSymmetry(TestBase):
    def setUp(self, N = 20, L = 1.):
        self.N = N
        self.L = L
        self.system = SymmetrySystem(N = self.N, L = self.L)
        
    def testResult(self):
        var = self.system.getVar()
        self.system.run()
        dx = self.L / self.N
        for j in range(self.N / 2):
            for i in range(self.N / 2):
                x = dx * (i + 0.5)
                y = dx * (j + 0.5)
                value = x * y
                self.assertWithinTolerance(value, var((self.L - x, y)).getNumericValue(), tol = 1e-10)
                self.assertWithinTolerance(value, var((x, self.L - y)).getNumericValue(), tol = 1e-10)
                self.assertWithinTolerance(value, var((self.L - x, self.L - y)).getNumericValue(), tol = 1e-10)

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestSymmetry))
    return theSuite
    
if __name__ == '__main__':
    fivol.tests.testProgram.main(defaultTest='suite')

            
            
