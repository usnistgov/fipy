#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:00:11 PM}
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

from fipy.tests.testBase import TestBase
import fipy.tests.testProgram

import fipy.examples.phase.examples.impingement
from fipy.examples.phase.examples.impingement.input1D import System1D
from fipy.examples.phase.examples.impingement.input4Particles import System4Particles
import fipy.tools.dump as dump

class TestImpingement(TestBase):
    def setUp(self):

        parameters = self.system.getParameters()

	self.steps = parameters['steps']
	self.tolerance = 1e-10

        self.it = parameters['it']
        self.var = parameters['theta']
        
    def getTestValues(self):
	filestream=os.popen('gunzip --fast -c < %s/%s'%(fipy.examples.phase.examples.impingement.__path__[0],self.testFile),'r')
	
	testData = cPickle.load(filestream)
	filestream.close()

	return testData

class Test1D(TestImpingement):
    def setUp(self):
        self.system = System1D(nx = 40, ny = 1)
        self.testFile = 'testImpingement.gz'
        TestImpingement.setUp(self)

class Test4Particles(TestImpingement):
    def setUp(self):
        self.system = System4Particles(nx = 20, ny = 20, drivingForce = 10.)
        self.testFile = '4ParticleData.gz'
        TestImpingement.setUp(self)

class Test4ParticlesRestart(TestImpingement):
    def setUp(self):
        self.testFile = '4ParticleData.gz'
        self.system = System4Particles(nx = 20, ny = 20, steps = 5, drivingForce = 10.)
        TestImpingement.setUp(self)
        for step in range(self.steps):
	    self.it.timestep()

        ## dump data
        dump.write(self.system.getDumpData(), 'restartData')

        ## run another bunch of steps
        restartData = dump.read('restartData')
        system = System4Particles(steps = 5, drivingForce = 10., restartData = restartData)
        TestImpingement.setUp(self)

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(Test1D))
    theSuite.addTest(unittest.makeSuite(Test4Particles))
    theSuite.addTest(unittest.makeSuite(Test4ParticlesRestart))
    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')

            
            
