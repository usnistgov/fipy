#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "testSteadyStateDiffusion.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 12/9/03 {2:42:58 PM} 
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
import os
import cPickle
import tests
from meshes.grid2D import Grid2D
from phase.phaseEquation import PhaseEquation
from solvers.linearPCGSolver import LinearPCGSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator
from phase.modularVariable import ModularVariable
from variables.cellVariable import CellVariable
import Numeric

class TestPhase(TestBase):
    """
    Simple test case for the phase field equation.
    """
    def setUp(self):
	self.steps = 100
	self.timeStep = 0.02
	self.tolerance = 1e-7

        phaseParameters={
            'tau' :        0.1,
            'epsilon' :    0.008,
            's' :          0.01,
            'alpha' :      0.015,
            'c2':          0.0,
            'anisotropy':  0.,
            'symmetry':    4.
            }
        
        valueLeft=1.
        valueRight=1.
        
        nx = self.nx
        ny = self.ny
        dx = self.L/nx
        dy = self.L/ny
        
        self.mesh = Grid2D(dx,dy,nx,ny)
        
        self.var = CellVariable(
            name = 'PhaseField',
            mesh = self.mesh,
            value = 1.
            )
	    
	phaseParameters['phi'] = self.var
        
        theta = ModularVariable(
            name = 'Theta',
            mesh = self.mesh,
            value = 1.,
            hasOld = 0
            )
	    
        func = self.func

        rightCells = self.mesh.getCells(func)
        
        theta.setValue(0.,rightCells)

	phaseParameters['theta'] = theta
	
	phaseParameters['temperature'] = 1.
	
        eq = PhaseEquation(
            self.var,
            solver = LinearPCGSolver(
		tolerance = 1.e-15, 
		steps = 1000
            ),
            boundaryConditions=(
		FixedValue(self.mesh.getFacesLeft(),valueLeft),
		FixedValue(self.mesh.getFacesRight(),valueRight)
	    ),
            parameters = phaseParameters
	)
        
        self.it = Iterator((eq,))
	
    def getTestValues(self):
	filestream=os.popen('gunzip --fast -c < %s/%s'%(tests.__path__[0],self.testFile),'r')
	
	testData = cPickle.load(filestream)
	filestream.close()

	return testData

class TestPhase1D(TestPhase):
    def setUp(self):
        self.nx = 100
        self.ny = 1
        L = self.L = 1.5
        def func(x,L=L):
            if x[0] > L / 2.:
                return 1
            else:
                return 0
        self.func = func
	self.testFile = 'testPhaseData.gz'
        TestPhase.setUp(self)

class TestPhaseCircle(TestPhase):
    def setUp(self):
        self.nx = 100
        self.ny = 100
        L = self.L = 1.5
        def func(x,L=L):
            r = L / 4.
            c = (L / 2., L / 2.)
            if (x[0] - c[0])**2 + (x[1] - c[1])**2 < r**2:
                return 1
            else:
                return 0
        self.func = func
	self.testFile = 'testCirclePhaseData.gz'
        TestPhase.setUp(self)

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestPhase1D))
    theSuite.addTest(unittest.makeSuite(TestPhaseCircle))
    return theSuite
    
if __name__ == '__main__':
    theSuite = suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)

            
            
