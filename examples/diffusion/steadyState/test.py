#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:05:59 PM} 
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

from fipy.tests.testBase import TestBase
import fipy.tests.testProgram

from input import getParameters

class Test(TestBase):
    """Generic steady-state diffusion class
    
    	Constructs a mesh, variable, equation, and iterator based
	on the mesh dimensions specified by the child class
    """
    def setUp(self):
        parameters = getParameters(self.nx, self.ny)
        
	self.steps = parameters['steps']
##	self.timeStep = parameters['timeStep']
	self.tolerance = parameters['tolerance']
        self.valueLeft = parameters['valueLeft']
        self.valueRight = parameters['valueRight']
        self.mesh = parameters['mesh']
        self.var = parameters['var']
##        self.eq = parameters['eq']
        self.it = parameters['it']

    def getTestValues(self):
	(lx,ly) = self.mesh.getPhysicalShape()
	vl = self.valueLeft
	vr = self.valueRight
	x = self.mesh.getCellCenters()[:,0]
	return vl + (vr - vl) * x / lx
	
class Test20x20(Test):
    """Steady-state 1D diffusion on a 20x20 mesh
    """
    def setUp(self):
	self.nx = 20
	self.ny = 20
	Test.setUp(self)	    
	
class Test50x50(Test):
    """Steady-state 1D diffusion on a 50x50 mesh
    """
    def setUp(self):
	self.nx = 50
	self.ny = 50
	Test.setUp(self)	    

class  Test1D(Test):
    """Steady-state 1D diffusion on a 50x1 mesh
    """
    def setUp(self):
	self.nx = 50
	self.ny = 1
	Test.setUp(self)

def suite():
    theSuite = unittest.TestSuite()
    
    theSuite.addTest(unittest.makeSuite(Test1D))
    theSuite.addTest(unittest.makeSuite(Test20x20))
    theSuite.addTest(unittest.makeSuite(Test50x50))
    
    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')

            
            
