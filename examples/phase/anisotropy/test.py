#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:02:10 PM} 
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

import examples.phase.anisotropy

from inputTemperature import AnisotropySystem

class TestAnisotropy(TestBase):
    def setUp(self):
        self.system = AnisotropySystem()
        self.testFile = 'testAnisotropy.gz'
        parameters = self.system.getParameters()

	self.steps = parameters['steps']
	self.tolerance = 1e-10

        self.it = parameters['it']
        self.var = parameters['var']
        
    def getTestValues(self):
	filestream=os.popen('gunzip --fast -c < %s/%s'%(examples.phase.anisotropy.__path__[0],self.testFile),'r')
	
	testData = cPickle.load(filestream)
	filestream.close()

	return testData

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestAnisotropy))
    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')

            
            
