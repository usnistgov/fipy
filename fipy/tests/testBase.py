#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "testBase.py"
 #                                    created: 12/5/03 {4:34:49 PM} 
 #                                last update: 1/25/04 {8:49:43 PM} 
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

import Numeric

class TestBase(unittest.TestCase):
    def assertWithinTolerance(self, first, second, tol = 1e-10, msg=None):
	"""Fail if the two objects are unequal by more than tol.
	"""
	if abs(first - second) > tol:
	    raise self.failureException, (msg or '%s !~ %s' % (first, second))

    def assertArrayWithinTolerance(self, first, second, atol = 1e-10, rtol = 1e-10, msg=None):
	"""Fail if the two objects are unequal by more than tol.
	"""
        
	if not Numeric.allclose(first, second, rtol, atol):
	    raise self.failureException, (msg or '\n%s\nis not\n%s' % (first, second))
        
    def getTestValue(self, cell):
	pass
	
    def getTestValues(self):
	values = self.var.getNumericValue().copy()
	for cell in self.mesh.getCells():
	    id = cell.getId()
	    values[id] = self.getTestValue(cell)
	return values
	
    def testResult(self):
	for step in range(self.steps):
	    self.it.timestep() #, maxSweeps = 10)
	array = self.var.getNumericValue()
	values = self.getTestValues()
	values = Numeric.reshape(values, Numeric.shape(array))
	self.assertArrayWithinTolerance(array, values, self.tolerance)
