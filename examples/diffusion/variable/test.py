#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "testVariableDiffusion.py"
 #                                    created: 11/26/03 {3:23:47 PM}
 #                                last update: 12/22/03 {3:56:32 PM} 
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

"""Test steady-state diffusion solutions
"""
 
import unittest
from tests.testBase import TestBase
import Numeric
import input

class TestVariableDiffusion(TestBase):
    """
    Variable diffusion test case. Diffusion
    varies in space. Test is for the steady state profile.
    Test case considers a bar with the following diffusion coefficient:
    
    1    when   0 < x < L/2
    0.1  when   L/2 < x < 3L/4
    1    when   3L/4 < x < L

    where L is the length of the bar. The boundary conditions are a fixed value
    of 0 on the left side and a fixed flux of 1 on the right side. A simple analytical
    answer to this diffusion problem is given in ansFunc(). In order to obtain very
    accurate numerical solutions it is necessary for cell centers to lie at exactly
    L/4 and 3L/4. Thus the test case should give accurate results for meshes of length
    2, 6, 10, 14, 18,...,50,... i.e. 4*n+2 where n is an integer.
    """
    
    def setUp(self):
        parameters = input.getParameters(self.nx, self.ny)
        self.it = parameters['it']
        self.steps = 1
        self.var = parameters['var']
        self.mesh = parameters['mesh']
        self.Lx = parameters['L']
        self.tolerance = 1e-8

    def getTestValues(self):
	x = self.mesh.getCellCenters()[:,0]
	L = self.Lx
	values = Numeric.where(x < 3. * L / 4., 10 * x - 9. * L / 4., x + 18. * L / 4.)
	values = Numeric.where(     x < L / 4.,                    x,           values)
	return values
	    
class TestVariableDiffusion2x1(TestVariableDiffusion):
    """Variable 1D diffusion on a 1x2 mesh
    """
    def setUp(self):
	self.nx = 2
	self.ny = 1
	TestVariableDiffusion.setUp(self)

	    
class TestVariableDiffusion10x1(TestVariableDiffusion):
    """Variable 1D diffusion on a 1x10 mesh
    """
    def setUp(self):
	self.nx = 10
	self.ny = 1
	TestVariableDiffusion.setUp(self)
	
class TestVariableDiffusion50x1(TestVariableDiffusion):
    """Variable 1D diffusion on a 1x50 mesh
    """
    def setUp(self):
	self.nx = 50
	self.ny = 1
	TestVariableDiffusion.setUp(self)

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestVariableDiffusion2x1))
    theSuite.addTest(unittest.makeSuite(TestVariableDiffusion10x1))
    theSuite.addTest(unittest.makeSuite(TestVariableDiffusion50x1))
    return theSuite
    
if __name__ == '__main__':
    theSuite = suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)

            
            
