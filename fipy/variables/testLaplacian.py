#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "testLaplacian.py"
 #                                    created: 4/30/04 {11:19:30 AM} 
 #                                last update: 5/6/04 {4:21:05 PM} 
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

##import sets

import unittest

import fipy.tests.testProgram
from fipy.tests.testBase import TestBase
from fipy.meshes.grid2D import Grid2D
from fipy.variables.cellVariable import CellVariable
import fipy.tools.array as array

class TestLaplacian(TestBase):
    def setUp(self, nx = 5, ny = 5, dx = 1., dy = 1.):
	self.mesh = Grid2D(nx = nx, ny = ny, dx = dx, dy = dx)
	self.var = CellVariable(
	    mesh = self.mesh,
	    value = self.getValue())
	    
    def testResult(self):
	result = self.var.getLaplacian(self.order)

	# boundary cells will be incorrect, so only look at interior
	# we need to strip off as many exterior layers of cells as
	# the order/2 of our Laplacian

        try:
            import sets
	
            cellFaceIDs = self.mesh.getCellFaceIDs()
            exteriorCellIDs = sets.Set(self.mesh.getExteriorCellIDs())
            interiorCellIDs = sets.Set(self.mesh.getInteriorCellIDs())
	
            for order in range(self.order/2 - 1):
                exteriorFaceIDs = sets.Set(array.take(self.mesh.getCellFaceIDs(), list(exteriorCellIDs)).flat)
                neighboringCellIDs = sets.Set([cellID for cellID in interiorCellIDs if len(sets.Set(cellFaceIDs[cellID]) & exteriorFaceIDs) != 0])
                exteriorCellIDs |= neighboringCellIDs
                interiorCellIDs -= neighboringCellIDs

        except:
            exteriorCellIDs = list(self.mesh.getExteriorCellIDs())
            interiorCellIDs = list(self.mesh.getInteriorCellIDs())
            cellToCellIDs = list(self.mesh.getCellToCellIDs())
        
            for order in range(self.order / 2 - 1):

                newExteriorCellIDs = []
                newInteriorCellIDs = []
                
                for cellID in interiorCellIDs:
                    inNewExteriorCellIDs = 0

                    for adjCellID in cellToCellIDs[cellID]:
                        if adjCellID in exteriorCellIDs:
                            inNewExteriorCellIDs = 1
                
                    if inNewExteriorCellIDs == 1:
                        newExteriorCellIDs.append(cellID)
                    else:
                        newInteriorCellIDs.append(cellID)

                exteriorCellIDs = newExteriorCellIDs
                interiorCellIDs = newInteriorCellIDs

        result = array.take(result, list(interiorCellIDs))
                
        self.assertArrayWithinTolerance(result, self.answer, 1e-10)

class TestLaplacian2(TestLaplacian):
    def setUp(self):
	self.answer = 4
	self.order = 2
	TestLaplacian.setUp(self)
	
    def getValue(self):
	centers = self.mesh.getCellCenters()
	return centers[:,0]**2 + centers[:,1]**2

class TestLaplacian4(TestLaplacian):
    def setUp(self):
	self.answer = 48
	self.order = 4
	TestLaplacian.setUp(self)
	
    def getValue(self):
	centers = self.mesh.getCellCenters()
	return centers[:,0]**4 + centers[:,1]**4

class TestLaplacian6(TestLaplacian):
    def setUp(self):
	self.answer = 1440
	self.order = 6
	TestLaplacian.setUp(self, nx = 7, ny = 7)
	
    def getValue(self):
	centers = self.mesh.getCellCenters()
	return centers[:,0]**6 + centers[:,1]**6

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestLaplacian2))
    theSuite.addTest(unittest.makeSuite(TestLaplacian4))
    theSuite.addTest(unittest.makeSuite(TestLaplacian6))
    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')

