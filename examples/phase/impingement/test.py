#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 01/21/04 {11:53:49 AM}
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
import os
import cPickle

from fivol.tests.testBase import TestBase
import fivol.examples.phase.examples.impingement

from input import ImpingementSystem

class TestImpingement(TestBase):
    def setUp(self):
        self.system = ImpingementSystem(nx = self.nx, ny = self.ny, initialConditions = self.initialConditions)
        self.testFile = 'testImpingement.gz'
        parameters = self.system.getParameters()

	self.steps = parameters['steps']
	self.tolerance = 1e-10

        self.it = parameters['it']
        self.var = parameters['theta']
        
    def getTestValues(self):
	filestream=os.popen('gunzip --fast -c < %s/%s'%(fivol.examples.phase.examples.impingement.__path__[0],self.testFile),'r')
	
	testData = cPickle.load(filestream)
	filestream.close()

	return testData

class Test1D(TestImpingement):
    def setUp(self):
        self.nx = 40
        self.ny = 1
        
        def getRightCells(cell, Lx = 1., Ly = 1.):
            if cell.getCenter()[0] > Lx / 2.:
                return 1

        def getAllCells(cell, Lx = 1., Ly = 1.):
            return 1.
    
        self.initialConditions = (
        { 'phase value' : 1., 'theta value' : 1., 'func' : getAllCells },
        { 'phase value' : 1., 'theta value' : 0., 'func' : getRightCells }        
        )

        TestImpingement.setUp(self)
        
def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(Test1D))
    return theSuite
    
if __name__ == '__main__':
    theSuite = suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)

            
            
