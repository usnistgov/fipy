#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 12/23/03 {6:26:24 PM} 
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
 
from __future__ import nested_scopes

import Numeric

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

	self.mesh = Grid2D(
	    dx = self.dx,
	    dy = self.dy,
	    nx = self.nx,
	    ny = self.ny)
	    
	setCells = self.mesh.getCells(self.func)
		
	fields = {}
        
	fields['phase'] = PhaseVariable(
	    name = "phase",
	    mesh = self.mesh,
	    value = self.phaseValue['left']
	    )
	    
	fields['phase'].setValue(self.phaseValue['right'],setCells)
# 	fields['phase'].setValue(0.9,(self.mesh.cells[17],))
# 	fields['phase'].setValue(0.8,(self.mesh.cells[18],))
# 	fields['phase'].setValue(0.51,(self.mesh.cells[19],))
# 	fields['phase'].setValue(0.49,(self.mesh.cells[20],))
# 	fields['phase'].setValue(0.2,(self.mesh.cells[21],))
# 	fields['phase'].setValue(0.1,(self.mesh.cells[22],))
	
	fields['substitutionals'] = ()
	
	for component in self.components:
	    component['var'] = ComponentVariable(
		name = component['name'],
		mesh = self.mesh,
		standardPotential = component['standard potential'],
		barrierHeight = component['barrier height'],
		value = component['left value']
		)
		
	    component['var'].setValue(component['right value'],setCells)
	    
	    fields['substitutionals'] += (component['var'],)
	    
	self.fields = fields
	
	parameters = {
	    'diffusivity': 1.,
	    'time step duration': self.timeStep,
	    'solvent standard potential': self.solvent['standard potential'],
	    'solvent barrier height': self.solvent['barrier height'],
	    'phase mobility': self.phaseMobility,
	    'phase gradient energy': self.phaseGradientEnergy
	}

	self.it = elphf.makeIterator(mesh = self.mesh, fields = fields, parameters = parameters)
	
    def testResult(self):
# 	self.fields['phase'].plot()
# 	raw_input()
	
	self.it.iterate(steps = self.steps)
	
# 	print "dx,dy: ", self.mesh.getCellDistances()
# 	print "vol: ", self.mesh.getCellVolumes()	
# # 	total = Numeric.array(self.fields['solvent'][:])
# # 	
# 	for component in self.components:
# 	    print "Cj: ", component['var'][:]
# 	    print "Cj avg.: ", Numeric.sum(component['var'][:])/(self.nx * self.ny)
# # 	    print "substitutionalSum: ", component['var'].substitutionalSum[:]
# # 	    print "weightedDiffusivity: ", component['var'].weightedDiffusivity[:]
# # 	    print "subsConvCoeff: ", component['var'].subsConvCoeff[:]
# # 	    print "pConvCoeff: ", component['var'].pConvCoeff[:]
# # 	    print "gConvCoeff: ", component['var'].gConvCoeff[:]
# # 	    print "convCoeff: ", component['var'].subsConvCoeff[:] + component['var'].pConvCoeff[:] + component['var'].gConvCoeff[:]
# # 
# # 	    total += component['var'][:]
# # 	    
# 	print "Cn: ", self.fields['solvent'][:]
# 	print "Cn avg.: ", Numeric.sum(self.fields['solvent'][:])/(self.nx * self.ny)
# # 	
# # 	print total

	self.fields['phase'].plot()
# 	print "phase: ", self.fields['phase']
# 	print "p: ", self.fields['phase'].get_p()[:]
# 	print "g: ", self.fields['phase'].get_g()[:]
# 	print "gFace: ", self.fields['phase'].get_gFace()[:]
# 	print "grad-p: ", self.fields['phase'].get_p().getFaceGrad()
# 	print "p' grad-phase: ", (30 * self.fields['phase'].get_gFace().transpose() * self.fields['phase'].getFaceGrad())[:]
	    
	for component in self.components:
	    component['var'].plot()
	    
	raw_input()
	    
	for component in self.components:
	    component['var'].plot()
	    array = component['var'].getValue()
	    self.assertArrayWithinTolerance(array, component['final value'], self.tolerance)	

class TestElPhF1D(TestElPhF):
    def setUp(self):
	self.steps = 40
	self.timeStep = 10000.
        self.nx = 40
        self.ny = 1
	self.dx = self.dy = 1.
	L = self.nx * self.dx
	
	self.phaseValue = {
	    'left': 1.,
	    'right': 1.
	}
	self.components = (
	    {
		'name': "c1",
		'standard potential': 1.,
		'barrier height': 1.,
		'left value': 0.3,
		'right value': 0.6,
		'final value': 0.45
	    },
	    {
		'name': "c2",
		'standard potential': 1.,
		'barrier height': 1.,
		'left value': 0.6,
		'right value': 0.3,
		'final value': 0.45
	    }
	)

	self.func = lambda center: center[0] > L/2
        TestElPhF.setUp(self)

class TestElPhF2D(TestElPhF):
    def setUp(self):
	self.steps = 40
	self.timeStep = 10000.
	self.nx = 40
	self.ny = 40
	self.dx = self.dy = 1.
	L = self.nx * self.dx
	
	self.phaseValue = {
	    'left': 1.,
	    'right': 1.
	}
	self.components = (
	    {
		'name': "c1",
		'standard potential': 1.,
		'barrier height': 1.,
		'left value': 0.3,
		'right value': 0.6,
		'final value': 0.45
	    },
	    {
		'name': "c2",
		'standard potential': 1.,
		'barrier height': 1.,
		'left value': 0.6,
		'right value': 0.3,
		'final value': 0.45
	    }
	)
	
	self.func = lambda center: center[0] > L/2
	TestElPhF.setUp(self)

class TestElPhF2DCorner(TestElPhF):
    def setUp(self):
	self.steps = 40
	self.timeStep = 10000.
	self.nx = 40
	self.ny = 40
	self.dx = self.dy = 1.
	L = self.nx * self.dx
	
	self.phaseValue = {
	    'left': 1.,
	    'right': 1.
	}
	self.components = (
	    {
		'name': "c1",
		'standard potential': 1.,
		'barrier height': 1.,
		'left value': 0.3,
		'right value': 0.6,
		'final value': 0.375
	    },
	    {
		'name': "c2",
		'standard potential': 1.,
		'barrier height': 1.,
		'left value': 0.6,
		'right value': 0.3,
		'final value': 0.525
	    }
	)
	
	self.func = lambda center: (center[0] > L/2) and (center[1] > L/2) 
	TestElPhF.setUp(self)
	
class TestElPhF1DphaseSeparate(TestElPhF):
    def setUp(self):
	self.steps = 40
	self.timeStep = 100000.
	self.nx = 1000
	self.ny = 1
	self.dx = self.dy = 0.01
	L = self.nx * self.dx
	
	self.phaseValue = {
	    'left': 1.,
	    'right': 0.
	}
	
	self.phaseMobility = 1.
	self.phaseGradientEnergy = 0.025
	
	self.solvent = {
	    'standard potential': Numeric.log(.1/.2),
	    'barrier height': 1.
	}
	
	self.components = (
	    {
		'name': "c1",
		'standard potential': Numeric.log(.3/.4) - self.solvent['standard potential'],
		'barrier height': 0.0,
		'left value': 0.35,
		'right value': 0.35,
		'final left': 0.4,
		'final right': 0.3
	    },
	    {
		'name': "c2",
		'standard potential': Numeric.log(.4/.3) - self.solvent['standard potential'],
		'barrier height': 0.0,
		'left value': 0.35,
		'right value': 0.35,
		'final left': 0.3,
		'final right': 0.4
	    },
	    {
		'name': "c3",
		'standard potential': Numeric.log(.2/.1) - self.solvent['standard potential'],
		'barrier height': 0.0,
		'left value': 0.15,
		'right value': 0.15,
		'final left': 0.1,
		'final right': 0.2
	    }
	)
	
	
	for component in self.components:
	    component['final value'] = Numeric.zeros((self.nx,self.ny),'d')
	    component['final value'][:self.nx/2] = component['final left']
	    component['final value'][self.nx/2:] = component['final right']
	
	self.func = lambda center: (center[0] > L/2)
	TestElPhF.setUp(self)
		
class TestElPhF1DphaseSeparateBoring(TestElPhF):
    def setUp(self):
	self.steps = 40
	self.timeStep = 100000.
	self.nx = 400
	self.ny = 1
	self.dx = self.dy = 0.01
	L = self.nx * self.dx
	
	self.phaseValue = {
	    'left': 1.,
	    'right': 0.
	}
	
	self.phaseMobility = 1.
	self.phaseGradientEnergy = 0.1
	
	self.solvent = {
	    'standard potential': Numeric.log(.7/.3),
	    'barrier height': 1.
	}
	
	self.components = ()
	
	self.components = (
	    {
		'name': "c1",
		'standard potential': Numeric.log(.3/.7) - self.solvent['standard potential'],
		'barrier height': 0., #100. - self.solvent['barrier height'],
		'left value': 0.5,
		'right value': 0.5,
		'final left': 0.7,
		'final right': 0.3
	    },
	)
	
	for component in self.components:
	    component['final value'] = Numeric.zeros((self.nx,self.ny),'d')
	    component['final value'][:self.nx/2] = component['final left']
	    component['final value'][self.nx/2:] = component['final right']
	
	self.func = lambda center: (center[0] > L/2)
	TestElPhF.setUp(self)
		
def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestElPhF1D))
#    theSuite.addTest(unittest.makeSuite(TestElPhF2D))
#    theSuite.addTest(unittest.makeSuite(TestElPhF2DCorner))
#    theSuite.addTest(unittest.makeSuite(TestElPhF1DphaseSeparate))
#    theSuite.addTest(unittest.makeSuite(TestElPhF1DphaseSeparateBoring))
    return theSuite
    
if __name__ == '__main__':
    theSuite = suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)

            
            
