#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 12/22/03 {4:04:04 PM} 
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
from meshes.grid2D import Grid2D
import elphf
from componentVariable import ComponentVariable
from phaseVariable import PhaseVariable

from solvers.linearPCGSolver import LinearPCGSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator
from examples.phase.modularVariable import ModularVariable
from variables.cellVariable import CellVariable

class TestElPhF(TestBase):
    """
    Simple test case for the phase field equation.
    """
    def setUp(self):
	self.tolerance = 1e-7

	self.L = self.nx * self.dx

	self.mesh = Grid2D(
	    dx = self.dx,
	    dy = self.dy,
	    nx = self.nx,
	    ny = self.ny)
        
	self.phase = PhaseVariable(
	    name = "phase",
	    mesh = self.mesh,
	    value = 1.
	    )
	    
	self.var1 = ComponentVariable(
	    name = "c1",
	    mesh = self.mesh,
	    standardPotential = 1.,
	    barrierHeight = 1.,
	    value = self.valueLeft
	    )

	self.var2 = ComponentVariable(
	    name = "c2",
	    mesh = self.mesh,
	    standardPotential = 1.,
	    barrierHeight = 1.,
	    value = self.valueRight
	    )
	    
	setCells = self.mesh.getCells(self.func)

	self.var1.setValue(self.valueRight,setCells)
	self.var2.setValue(self.valueLeft,setCells)
	    
	fields = {
	    'phase': self.phase,
	    'substitutionals': (self.var1,self.var2)
	}
	
	parameters = {
	    'diffusivity': 1.,
	    'time step duration': self.timeStep
	}

	self.it = elphf.makeIterator(mesh = self.mesh, fields = fields, parameters = parameters)
	
    def testResult(self):
	self.it.iterate(steps = self.steps)
	
	array = self.var1.getValue()
	self.assertArrayWithinTolerance(array, self.val1, self.tolerance)	
	array = self.var2.getValue()
	self.assertArrayWithinTolerance(array, self.val2, self.tolerance)	

class TestElPhF1D(TestElPhF):
    def setUp(self):
	self.steps = 40
	self.timeStep = 10000.
        self.nx = 40
        self.ny = 1
	self.dx = self.dy = 1.
	self.valueLeft=0.3
	self.valueRight=0.6
# 	self.valueLeft="0.3 mol/l"
# 	self.valueRight="0.6 mol/l"
	self.val1 = self.val2 = 0.45
        self.L = self.nx * self.dx
	self.func = lambda center, L = self.L: center[0] > L/2
        TestElPhF.setUp(self)

class TestElPhF2D(TestElPhF):
    def setUp(self):
	self.steps = 40
	self.timeStep = 10000.
	self.nx = 40
	self.ny = 40
	self.dx = self.dy = 1.
	self.valueLeft=0.3
	self.valueRight=0.6
# 	self.valueLeft="0.3 mol/l"
# 	self.valueRight="0.6 mol/l"
	self.val1 = self.val2 = 0.45
        self.L = self.nx * self.dx
	self.func = lambda center, L = self.L: center[0] > L/2
	TestElPhF.setUp(self)

class TestElPhF2DCorner(TestElPhF):
    def setUp(self):
	self.steps = 40
	self.timeStep = 10000.
	self.nx = 40
	self.ny = 40
	self.dx = self.dy = 1.
	self.valueLeft=0.3
	self.valueRight=0.6
# 	self.valueLeft="0.3 mol/l"
# 	self.valueRight="0.6 mol/l"
	self.val1 = 0.375
	self.val2 = 0.525
        self.L = self.nx * self.dx
	self.func = lambda center, L = self.L: (center[0] > L/2) and (center[1] > L/2) 
	TestElPhF.setUp(self)
	
def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestElPhF1D))
    theSuite.addTest(unittest.makeSuite(TestElPhF2D))
    theSuite.addTest(unittest.makeSuite(TestElPhF2DCorner))
    return theSuite
    
if __name__ == '__main__':
    theSuite = suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)

            
            
