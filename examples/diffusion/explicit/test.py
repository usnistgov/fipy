#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "testExplicitDiffusion.py"
 #                                    created: 11/27/03 {3:23:47 PM}
 #                                last update: 1/13/04 {1:01:01 PM} 
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
from tests.testBase import TestBase
import input

class TestExplicitDiffusion(TestBase):
    """Generic steady-state diffusion class
    Same as TestSteadyStateDiffusion biut for explcit case
    Constructs a mesh, variable, equation, and iterator based
    on the mesh dimensions specified by the child class
    """
    def setUp(self):
        parameters = input.getParameters(self.nx, self.ny)
	self.steps = parameters['steps']
	self.timeStep = parameters['timeStep']
	self.tolerance = parameters['tolerance']
        self.mesh = parameters['mesh']
        self.var = parameters['var']
        self.it = parameters['it']
        self.valueLeft = parameters['valueLeft']
        self.valueRight = parameters['valueRight']

    def getTestValues(self):
	(lx,ly) = self.mesh.getPhysicalShape()
	vl = self.valueLeft
	vr = self.valueRight
	x = self.mesh.getCellCenters()[:,0]
	return vl + (vr - vl) * x / lx
	
class  TestExplicitDiffusion10(TestExplicitDiffusion):
    """Steady-state 1D diffusion on a 10x1 mesh
    """
    def setUp(self):
	self.nx = 10
	self.ny = 1
	TestExplicitDiffusion.setUp(self)

class  TestExplicitDiffusion50(TestExplicitDiffusion):
    """Steady-state 1D diffusion on a 50x1 mesh
    """
    def setUp(self):
	self.nx = 50
	self.ny = 1
	TestExplicitDiffusion.setUp(self)

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestExplicitDiffusion10))
    theSuite.addTest(unittest.makeSuite(TestExplicitDiffusion50))
    return theSuite
    
if __name__ == '__main__':
    theSuite = suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)
    
            
            
