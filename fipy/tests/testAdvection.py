#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "testAdvection.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 12/5/03 {4:53:50 PM} 
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

"""Test steady convection diffusion solutions
"""

import unittest
from testBase import TestBase
from meshes.grid2D import Grid2D
from equations.advectionEquation import AdvectionEquation
from solvers.linearCGSSolver import LinearCGSSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator
from variables.variable import Variable
from terms.exponentialConvectionTerm import ExponentialConvectionTerm
from terms.powerLawConvectionTerm import PowerLawConvectionTerm
import Numeric

class TestAdvection(TestBase):
    """Generic advection test class
    
    	Constructs a mesh, variable, equation, and iterator based
	on the mesh dimensions specified by the child class
    """
    def setUp(self):
	self.steps = 1
	self.timeStep = 1.
	
        self.valueLeft = 0.
        self.valueRight = 1.

        self.mesh = Grid2D(self.L/self.nx,1.,self.nx,self.ny)
        
        self.var = Variable(
            name = "concentration",
            mesh = self.mesh,
	    value = self.valueLeft,
            viewer = 'None')

	self.eq = AdvectionEquation(
	    var = self.var,
	    transientCoeff = self.tranCoeff,
	    convectionCoeff = self.convCoeff,
	    solver = LinearCGSSolver(
		tolerance = 1.e-15, 
		steps = 1000
	    ),
	    convectionScheme = self.convectionScheme,
	    boundaryConditions=(
		FixedValue(faces = self.mesh.getFacesLeft(),value = self.valueLeft),
		FixedValue(faces = self.mesh.getFacesRight(),value = self.valueRight),
		FixedFlux(faces = self.mesh.getFacesTop(),value = 0.),
		FixedFlux(faces = self.mesh.getFacesBottom(),value = 0.)
	    )
	)

        self.it = Iterator((self.eq,))

    def getTestValue(self, coords):
# 	(lx,ly) = self.mesh.getPhysicalShape()
# 	vl = self.valueLeft
# 	vr = self.valueRight
	x = coords[0]
	val = (1. - Numeric.exp(-self.convCoeff[0] * x / self.diffCoeff)) / (1. - Numeric.exp(-self.convCoeff[0] * self.L / self.diffCoeff))
	
	return val
	
            
	    
class  TestSteadyConvectionDiffusion1DExponential(TestSteadyConvectionDiffusion):
    """Steady-state 1D diffusion on a 100x1 mesh, with exponentional convection scheme
    """
    def setUp(self):
	self.L = 10.
	self.nx = 1000
	self.ny = 1
	self.diffCoeff = 1.
	self.convCoeff = (10.,0)
	self.tolerance = 1e-10
	self.convectionScheme = ExponentialConvectionTerm
	TestSteadyConvectionDiffusion.setUp(self)

class  TestSteadyConvectionDiffusion2DExponential(TestSteadyConvectionDiffusion):
    """Steady-state 1D diffusion on a 10x10 mesh, with exponentional convection scheme
    """
    def setUp(self):
	self.L = 10.
	self.nx = 10
	self.ny = 10
	self.diffCoeff = 1.
	self.convCoeff = (10.,0)
	self.tolerance = 1e-10
	self.convectionScheme = ExponentialConvectionTerm
	TestSteadyConvectionDiffusion.setUp(self)
	
class  TestSteadyConvectionDiffusion1DExponentialBackwards(TestSteadyConvectionDiffusion):
    """Steady-state 1D diffusion on a 100x1 mesh, with exponentional convection scheme
    """
    def setUp(self):
	self.L = 10.
	self.nx = 1000
	self.ny = 1
	self.diffCoeff = 1.
	self.convCoeff = (-10.,0)
	self.tolerance = 1e-10
	self.convectionScheme = ExponentialConvectionTerm
	TestSteadyConvectionDiffusion.setUp(self)
	
class  TestSteadyConvectionDiffusion1DPowerLaw(TestSteadyConvectionDiffusion):
    """Steady-state 1D diffusion on a 100x1 mesh, with power law convection scheme
    """
    def setUp(self):
	self.L = 10.
	self.nx = 1000
	self.ny = 1
	self.diffCoeff = 1.
	self.convCoeff = (10.,0)
	self.tolerance = 1e-2
	self.convectionScheme = PowerLawConvectionTerm
	TestSteadyConvectionDiffusion.setUp(self)
	
def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestSteadyConvectionDiffusion1DExponential))
    theSuite.addTest(unittest.makeSuite(TestSteadyConvectionDiffusion1DPowerLaw))
    theSuite.addTest(unittest.makeSuite(TestSteadyConvectionDiffusion1DExponentialBackwards))
    theSuite.addTest(unittest.makeSuite(TestSteadyConvectionDiffusion2DExponential))
    return theSuite
    
if __name__ == '__main__':
    theSuite = suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)

            
            
