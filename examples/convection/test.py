#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "testStdyConvectionDiffusion.py"
 #                                    created: 12/16/03 {3:23:47 PM}
 #                                last update: 1/16/04 {11:26:15 AM} 
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
 
import unittest

import Numeric

from fivol.tests.testBase import TestBase

from exponential1DSystem import Exponential1DSystem
from exponential2DSystem import Exponential2DSystem
from exponential1DBackSystem import Exponential1DBackSystem
from exponential1DUpSystem import Exponential1DUpSystem
from power1DSystem import Power1DSystem
from exponential1DScSystem import Exponential1DScSystem

class TestSteadyConvectionDiffusionSc(TestBase):
    """steady-state convection-diffusion-source 
    """
    def setUp(self):
        parameters = self.system.getParameters()
	self.steps = parameters['steps']
        self.var = parameters['var']
        self.it = parameters['it']
	self.mesh = parameters['mesh']
        self.convCoeff = parameters['convection coeff']
        self.diffCoeff = parameters['diffusion coeff']
        self.sourceCoeff = parameters['source coeff']
        self.L = parameters['L']
        
    def getTestValues(self):
	if self.convCoeff[0] != 0.:
	    axis = 0
	else:
	    axis = 1
	x = self.mesh.getCellCenters()[:,axis]
        AA = -self.sourceCoeff * x / self.convCoeff[axis]
        BB = 1. + self.sourceCoeff * self.L / self.convCoeff[axis]
        CC = 1. - Numeric.exp(-self.convCoeff[axis] * x / self.diffCoeff)
        DD = 1. - Numeric.exp(-self.convCoeff[axis] * self.L / self.diffCoeff)
	return AA + BB * CC / DD
	
	    
class  TestSteadyConvectionDiffusion1DExponential(TestSteadyConvectionDiffusionSc):
    """Steady-state 1D diffusion on a 100x1 mesh, with exponentional convection scheme
    """
    def setUp(self):
        self.system = Exponential1DSystem()
        self.tolerance = 1e-10
	TestSteadyConvectionDiffusionSc.setUp(self)

class  TestSteadyConvectionDiffusion2DExponential(TestSteadyConvectionDiffusionSc):
    """Steady-state 1D diffusion on a 10x10 mesh, with exponentional convection scheme
    """
    def setUp(self):
        self.system = Exponential2DSystem()
        self.tolerance = 1e-10
	TestSteadyConvectionDiffusionSc.setUp(self)
	
class  TestSteadyConvectionDiffusion1DExponentialBackwards(TestSteadyConvectionDiffusionSc):
    """Steady-state 1D diffusion on a 100x1 mesh, with exponentional convection scheme
    """
    def setUp(self):
        self.system = Exponential1DBackSystem()
        self.tolerance = 1e-10
	TestSteadyConvectionDiffusionSc.setUp(self)
	
class  TestSteadyConvectionDiffusion1DExponentialUp(TestSteadyConvectionDiffusionSc):
    """Steady-state 1D diffusion on a 100x1 mesh, with exponentional convection scheme
    """
    def setUp(self):
        self.system = Exponential1DBackSystem()
        self.tolerance = 1e-10
	TestSteadyConvectionDiffusionSc.setUp(self)
	
class  TestSteadyConvectionDiffusion1DPowerLaw(TestSteadyConvectionDiffusionSc):
    """Steady-state 1D diffusion on a 100x1 mesh, with power law convection scheme
    """
    def setUp(self):
        self.system = Power1DSystem()
	self.tolerance = 1e-2
	TestSteadyConvectionDiffusionSc.setUp(self)
        
class  TestSteadyConvectionDiffusion1DExponentialSc(TestSteadyConvectionDiffusionSc):
    """Steady-state 1D diffusion on a 100x1 mesh, with exponentional convection scheme
    """
    def setUp(self):
        self.system = Exponential1DScSystem()
	self.tolerance = 1e-6
	TestSteadyConvectionDiffusionSc.setUp(self)

    def getBoundaryConditions(self):
	return (
	    FixedFlux(faces = self.mesh.getFacesLeft(),value = 0.),
	    FixedFlux(faces = self.mesh.getFacesRight(),value = 0.),
	    FixedValue(faces = self.mesh.getFacesTop(),value = self.valueRight),
	    FixedValue(faces = self.mesh.getFacesBottom(),value = self.valueLeft)
	)
	

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestSteadyConvectionDiffusion1DExponential))
    theSuite.addTest(unittest.makeSuite(TestSteadyConvectionDiffusion1DExponentialUp))
    theSuite.addTest(unittest.makeSuite(TestSteadyConvectionDiffusion1DPowerLaw))
    theSuite.addTest(unittest.makeSuite(TestSteadyConvectionDiffusion1DExponentialBackwards))
    theSuite.addTest(unittest.makeSuite(TestSteadyConvectionDiffusion2DExponential))
    theSuite.addTest(unittest.makeSuite(TestSteadyConvectionDiffusion1DExponentialSc))
    return theSuite
    
if __name__ == '__main__':
    theSuite = suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)

            
            
