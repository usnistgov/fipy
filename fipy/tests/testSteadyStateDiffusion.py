#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "testSteadyStateDiffusion.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 12/10/03 {2:28:16 PM} 
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
from testBase import TestBase
from meshes.grid2D import Grid2D
from equations.diffusionEquation import DiffusionEquation
from solvers.linearPCGSolver import LinearPCGSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator
from variables.cellVariable import CellVariable

class TestSteadyStateDiffusion(TestBase):
    """Generic steady-state diffusion class
    
    	Constructs a mesh, variable, equation, and iterator based
	on the mesh dimensions specified by the child class
    """
    def setUp(self):
	self.steps = 1
	self.timeStep = 1.
	self.tolerance = 1e-8

        self.valueLeft = 0.
        self.valueRight = 1.

        self.mesh = Grid2D(1.,1.,self.nx,self.ny)
        
        self.var = CellVariable(
            name = "concentration",
            mesh = self.mesh,
	    value = self.valueLeft,
            viewer = None)

        self.eq = DiffusionEquation(
            self.var,
            transientCoeff = 0., 
            diffusionCoeff = 1.,
            solver = LinearPCGSolver(
            tolerance = 1.e-15, 
            steps = 1000
            ),
            boundaryConditions=(
            FixedValue(self.mesh.getFacesLeft(),self.valueLeft),
            FixedValue(self.mesh.getFacesRight(),self.valueRight),
            FixedFlux(self.mesh.getFacesTop(),0.),
            FixedFlux(self.mesh.getFacesBottom(),0.)
            )
            )

        self.it = Iterator((self.eq,))

    def getTestValues(self):
	(lx,ly) = self.mesh.getPhysicalShape()
	vl = self.valueLeft
	vr = self.valueRight
	x = self.mesh.getCellCenters()[:,0]
	return vl + (vr - vl) * x / lx
	
class TestSteadyStateDiffusion20x20(TestSteadyStateDiffusion):
    """Steady-state 1D diffusion on a 20x20 mesh
    """
    def setUp(self):
	self.nx = 20
	self.ny = 20
	TestSteadyStateDiffusion.setUp(self)	    
	
class TestSteadyStateDiffusion50x50(TestSteadyStateDiffusion):
    """Steady-state 1D diffusion on a 50x50 mesh
    """
    def setUp(self):
	self.nx = 50
	self.ny = 50
	TestSteadyStateDiffusion.setUp(self)	    

class  TestSteadyStateDiffusion1D(TestSteadyStateDiffusion):
    """Steady-state 1D diffusion on a 50x1 mesh
    """
    def setUp(self):
	self.nx = 50
	self.ny = 1
	TestSteadyStateDiffusion.setUp(self)

class  TestSteadyStateDiffusion1D(TestSteadyStateDiffusion):
    """Steady-state 1D diffusion on a 50x1 mesh
    """
    def setUp(self):
	self.nx = 50
	self.ny = 1
	TestSteadyStateDiffusion.setUp(self)
	
def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestSteadyStateDiffusion1D))
    theSuite.addTest(unittest.makeSuite(TestSteadyStateDiffusion20x20))
    theSuite.addTest(unittest.makeSuite(TestSteadyStateDiffusion50x50))
    return theSuite
    
if __name__ == '__main__':
    theSuite = suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)

            
            
