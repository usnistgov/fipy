## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "testBase.py"
 #                                    created: 12/5/03 {4:34:49 PM} 
 #                                last update: 12/5/03 {5:12:02 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
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

class TestBase(unittest.TestCase):
    def getTestValue(self, cell):
	pass
	
    def assertWithinTolerance(self, first, second, tol = 1e-10, msg=None):
	"""Fail if the two objects are unequal by more than tol.
	"""
	if abs(first - second) > tol:
	    raise self.failureException, (msg or '%s !~ %s' % (first, second))
	    
    def testResult(self):
	self.it.iterate(steps = self.steps, timeStep = self.timeStep)
	array = self.var.getArray()

	for cell in self.mesh.getCells():
	    id = cell.getId()
	    val = self.getTestValue(cell.getCenter())
	    norm = abs(array[id] - val)
	    self.assertWithinTolerance(norm, 0.0, self.tolerance,("cell(%g)'s value of %g differs from %g by %g" % (id,array[id],val,norm)))

