#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "testVariable.py"
 #                                    created: 2/20/04 {11:19:30 AM} 
 #                                last update: 2/20/04 {12:21:35 PM} 
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

import unittest

from Numeric import array

import fivol.tests.testProgram
from fivol.tests.testBase import TestBase
import Numeric
from fivol.meshes.grid2D import Grid2D
from fivol.variables.cellVariable import CellVariable
import fivol.tools.dump as dump

class TestVariablePickle(TestBase):
    def setUp(self):
        mesh = Grid2D(dx = 1., dy = 1., nx = 10, ny = 10)

	self.var = CellVariable(mesh = mesh, value = 1., hasOld = 1, name = 'test')

        self.var.setValue(mesh.getCellCenters()[:,0] * mesh.getCellCenters()[:,1])

        pickledVar = dump.write(self.var, 'pickledVar')
        self.unPickledVar = dump.read('pickledVar')
        
    def testResult(self):
        pass

    def testValue(self):
	self.assertWithinTolerance(self.var.getValue(), self.unPickledVar.getValue(), 1e-10)

    def testOldValue(self):
        self.assertWithinTolerance(self.var.getOld().getValue(), self.unPickledVar.getOld().getValue(), 1e-10)
	
def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestVariablePickle))
    return theSuite
    
if __name__ == '__main__':
    fivol.tests.testProgram.main(defaultTest='suite')

