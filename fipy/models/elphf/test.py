#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 4/1/05 {2:50:32 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
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

import Numeric

import fipy.tests.testProgram

from fipy.tests.testBase import TestBase

class TestElPhF(TestBase):
    """
    Simple test case for the phase field equation.
    """
    def setUp(self, input):
	self.mesh = input.mesh
	self.it = input.it
	self.fields = input.fields
	self.parameters = input.parameters
	
	self.tolerance = 1e-7
	self.steps = 40

	self.final = {
	    'phase': [1],
	    'potential': [0],
	    'substitutionals': [],
	    'interstitials': []
	}
	
    def assertFieldWithinTolerance(self, field, final):
	self.assertWithinTolerance(field[0], final[0], self.tolerance)	
	self.assertWithinTolerance(field[-1], final[-1], self.tolerance)	

    def testResult(self):
	for step in range(self.steps):
	    self.it.timestep()

	self.assertFieldWithinTolerance(self.fields['phase'], self.final['phase'])	
	self.assertFieldWithinTolerance(self.fields['potential'], self.final['potential'])	
	
	for i in range(len(self.fields['substitutionals'])):
	    self.assertFieldWithinTolerance(self.fields['substitutionals'][i], self.final['substitutionals'][i])	
	
	for i in range(len(self.fields['interstitials'])):
	    self.assertFieldWithinTolerance(self.fields['interstitials'][i], self.final['interstitials'][i])	
		
class TestElPhF1D(TestElPhF):
    def setUp(self):
	import input1D
	TestElPhF.setUp(self, input1D)
	
	self.final['substitutionals'] = [[0.45],[0.45]]

class TestElPhF2D(TestElPhF):
    def setUp(self):
	import input2D
	TestElPhF.setUp(self, input2D)
	
	self.final['substitutionals'] = [[0.45],[0.45]]
    
class TestElPhF2Dcorner(TestElPhF):
    def setUp(self):
	import input2Dcorner
	TestElPhF.setUp(self, input2Dcorner)
	
	self.final['substitutionals'] = [[0.375],[0.525]]

class TestElPhF1Dphase(TestElPhF):
    def setUp(self):
	import input1Dphase
	TestElPhF.setUp(self, input1Dphase)
	
	self.tolerance = 1e-4
	self.steps = 10
	
    def testResult(self):
	for step in range(self.steps):
	    self.it.timestep()

	field  = self.parameters['phase']
    
	x = Numeric.arange(float(self.parameters['mesh']['nx'])) 
	x -= (self.parameters['mesh']['nx'] - 1.) / 2.
	x *= self.parameters['mesh']['dx']
	d = Numeric.sqrt(field['gradient energy'] / (self.parameters['solvent']['barrier height']))
	final = (1. - Numeric.tanh(x/(2.*d))) / 2.
	
	self.assertArrayWithinTolerance(field['var'].getNumericValue(), final, self.tolerance)
	
class TestElPhF1DphaseBinary(TestElPhF):
    def setUp(self):
	import input1DphaseBinary
	TestElPhF.setUp(self, input1DphaseBinary)
	
	self.tolerance = 2e-3
	
	self.final['phase'] = [1.,0.]
	self.final['substitutionals'] = [[0.7,0.3]]

class TestElPhF1DphaseQuaternary(TestElPhF):
    def setUp(self):
	import input1DphaseQuaternary
	TestElPhF.setUp(self, input1DphaseQuaternary)
	
	self.tolerance = 2e-3
	
	self.final['phase'] = [1.,0.]
	self.final['substitutionals'] = [[0.4,0.3],[0.3,0.4],[0.1,0.2]]
	
class TestElPhF1DphaseTernaryAndElectrons(TestElPhF):
    def setUp(self):
	import input1DphaseTernAndElectrons
	TestElPhF.setUp(self, input1DphaseTernAndElectrons)
	
	self.tolerance = 2e-3
	
	self.final['phase'] = [1.,0.]
	self.final['interstitials'] = [[0.4,0.3]]
	self.final['substitutionals'] = [[0.3,0.4],[0.1,0.2]]
	
class TestElPhF1DpoissonAllCharge(TestElPhF):
    def setUp(self):
	import input1DpoissonAllCharge
	TestElPhF.setUp(self, input1DpoissonAllCharge)
	
	self.tolerance = 2e-5
	
    def testResult(self):
	for step in range(self.steps):
	    self.it.timestep()

	x = Numeric.arange(float(self.parameters['mesh']['nx'])) 
	x += 0.5
	x *= self.parameters['mesh']['dx']
	final = (x**2)/2 - 2*x
	
	self.assertArrayWithinTolerance(self.fields['potential'].getNumericValue(), final, self.tolerance)
		
class TestElPhF1DpoissonLeftCharge(TestElPhF):
    def setUp(self):
	import input1DpoissonLeftCharge
	TestElPhF.setUp(self, input1DpoissonLeftCharge)
	
	self.tolerance = 2e-5
	
    def testResult(self):
	for step in range(self.steps):
	    self.it.timestep()

	x = Numeric.arange(float(self.parameters['mesh']['nx'])) 
	x += 0.5
	x *= self.parameters['mesh']['dx']
	final = Numeric.where(x < 1, (x**2)/2 - x, -0.5)
	
	self.assertArrayWithinTolerance(self.fields['potential'].getNumericValue(), final, self.tolerance)
    
class TestElPhF1DpoissonRightCharge(TestElPhF):
    def setUp(self):
	import input1DpoissonRightCharge
	TestElPhF.setUp(self, input1DpoissonRightCharge)
	
	self.tolerance = 2e-5
	
    def testResult(self):
	for step in range(self.steps):
	    self.it.timestep()

	x = Numeric.arange(float(self.parameters['mesh']['nx'])) 
	x += 0.5
	x *= self.parameters['mesh']['dx']
	final = Numeric.where(x < 1, -x, ((x-1)**2)/2 - x)
	
	self.assertArrayWithinTolerance(self.fields['potential'].getNumericValue(), final, self.tolerance)

def _suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestElPhF1D))
    theSuite.addTest(unittest.makeSuite(TestElPhF2D))
    theSuite.addTest(unittest.makeSuite(TestElPhF2Dcorner))    
    theSuite.addTest(unittest.makeSuite(TestElPhF1Dphase))
    theSuite.addTest(unittest.makeSuite(TestElPhF1DphaseBinary))
    theSuite.addTest(unittest.makeSuite(TestElPhF1DphaseQuaternary))
    theSuite.addTest(unittest.makeSuite(TestElPhF1DphaseTernaryAndElectrons))
    theSuite.addTest(unittest.makeSuite(TestElPhF1DpoissonAllCharge))
    theSuite.addTest(unittest.makeSuite(TestElPhF1DpoissonLeftCharge))
    theSuite.addTest(unittest.makeSuite(TestElPhF1DpoissonRightCharge))

    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')

            
            
