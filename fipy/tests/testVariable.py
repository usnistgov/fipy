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

from fivol.meshes.grid2D import Grid2D
##from fivol.meshes.vertex import Vertex
from fivol.variables.cellVariable import CellVariable
from fivol.variables.arithmeticCellToFaceVariable import ArithmeticCellToFaceVariable
from fivol.variables.harmonicCellToFaceVariable import HarmonicCellToFaceVariable

class TestMesh(Grid2D):
    def __init__(self, dx, dy, nx, ny, factor):
	self.factor = factor
	Grid2D.__init__(self, dx, dy, nx, ny)
	
##    def createVertices(self):
##	"""Return list of Vertices
##	"""
##	vertices = ()
##	ny=self.ny
##	nx=self.nx
##	dx=self.dx
##	dy=self.dy
##	for j in range(ny+1):
##	    dx = self.dx
##	    for	i in range(nx+1):
##		vertices += (Vertex(array([i * dx, j * dy],'d')),)
##		dx = dx * self.factor
##	return vertices	

class TestMean(TestBase):
    def setUp(self, value, dx = 1., dy = 1., factor = 1):
	self.mesh = TestMesh(dx, dy, 2, 1, factor = factor)
	self.var = CellVariable(
	    mesh = self.mesh,
	    value = value)
	self.alpha = 1. / (2. * factor)

    def testResult(self):
	self.assertWithinTolerance(self.result[0], self.answer, 1e-10)
	
class TestArithmeticMean(TestMean):
    def setUp(self, value, dx = 1., dy = 1., factor = 1):
	TestMean.setUp(self, array(value,'d'), dx, dy, factor)
	self.mean = ArithmeticCellToFaceVariable(self.var)
	self.result = self.mean.getNumericValue()
	self.answer = (value[0] - value[1]) * self.alpha + value[1]
	
class TestArithmeticMean1(TestArithmeticMean):
    def setUp(self):
	TestArithmeticMean.setUp(self, value = (1,2), factor = 1)

class TestArithmeticMean2(TestArithmeticMean):
    def setUp(self):
	TestArithmeticMean.setUp(self, value = (1,2), factor = 2)
	
class TestArithmeticMean3(TestArithmeticMean):
    def setUp(self):
	TestArithmeticMean.setUp(self, value = (1,2), factor = 10)

class TestHarmonicMean(TestMean):
    def setUp(self, value, dx = 1., dy = 1., factor = 1):
	TestMean.setUp(self, value, dx, dy, factor)
	self.mean = HarmonicCellToFaceVariable(self.var)
	self.result = self.mean.getNumericValue()
	self.answer = value[0] * value[1] / ((value[1] - value[0]) * self.alpha + value[0])

class TestHarmonicMean1(TestHarmonicMean):
    def setUp(self):
	TestHarmonicMean.setUp(self, value = (1,2), factor = 1)

class TestHarmonicMean2(TestHarmonicMean):
    def setUp(self):
	TestHarmonicMean.setUp(self, value = (1,2), factor = 2)
	
class TestHarmonicMean3(TestHarmonicMean):
    def setUp(self):
	TestHarmonicMean.setUp(self, value = (1,2), factor = 10)

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestArithmeticMean1))
    theSuite.addTest(unittest.makeSuite(TestArithmeticMean2))
    theSuite.addTest(unittest.makeSuite(TestArithmeticMean3))
    theSuite.addTest(unittest.makeSuite(TestHarmonicMean1))
    theSuite.addTest(unittest.makeSuite(TestHarmonicMean2))
    theSuite.addTest(unittest.makeSuite(TestHarmonicMean3))
    return theSuite
    
if __name__ == '__main__':
    fivol.tests.testProgram.main(defaultTest='suite')

