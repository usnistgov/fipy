#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based phase field solver
 # 
 #  FILE: "testVariable.py"
 #                                    created: 2/20/04 {11:19:30 AM} 
 #                                last update: 12/9/04 {6:14:23 PM} 
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

import unittest

from Numeric import array

import fipy.tests.testProgram
from fipy.tests.testBase import TestBase
import Numeric
from fipy.meshes.grid2D import Grid2D
from fipy.variables.cellVariable import CellVariable
import fipy.tools.dump as dump

class TestVariablePickle(TestBase):
    def setUp(self):
        mesh = Grid2D(dx = 1., dy = 1., nx = 10, ny = 10)

	self.var = CellVariable(mesh = mesh, value = 1., hasOld = 1, name = 'test')

        self.var.setValue(mesh.getCellCenters()[:,0] * mesh.getCellCenters()[:,1])

        import tempfile
        import os
        tmp = tempfile.gettempdir()
        fileName = os.path.join(tmp, 'data')
        pickledVar = dump.write(self.var, fileName)
        self.unPickledVar = dump.read(fileName)
        
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
    fipy.tests.testProgram.main(defaultTest='suite')

