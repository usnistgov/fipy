#!/usr/bin/env python

"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "test.py"
                                   created: 11/10/03 {3:23:47 PM}
                               last update: 11/20/03 {4:57:42 PM} 
 Author: Jonathan Guyer
 E-mail: guyer@nist.gov
   mail: NIST
    www: http://ctcms.nist.gov
 
========================================================================
This software was developed at the National Institute of Standards
and Technology by employees of the Federal Government in the course
of their official duties.  Pursuant to title 17 Section 105 of the
United States Code this software is not subject to copyright
protection and is in the public domain.  PFM is an experimental
system.  NIST assumes no responsibility whatsoever for its use by
other parties, and makes no guarantees, expressed or implied, about
its quality, reliability, or any other characteristic.  We would
appreciate acknowledgement if the software is used.

This software can be redistributed and/or modified freely
provided that any derivative works bear some notice that they are
derived from it, and any modified versions bear some notice that
they have been modified.
========================================================================
 
 Description: 

 History

 modified   by  rev reason
 ---------- --- --- -----------
 2003-11-10 JEG 1.0 original
###################################################################
"""

import unittest
from meshes.grid2D import Grid2D
from equations.diffusionEquation import DiffusionEquation
from solvers.linearPCGSolver import LinearPCGSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator
from variables.variable import Variable

class TestSteadyStateDiffusion(unittest.TestCase):
    def setUp(self):

        self.valueLeft = 0.
        self.valueRight = 1.

        self.mesh = Grid2D(1.,1.,self.nx,self.ny)
        
        self.var = Variable(
            name = "concentration",
            mesh = self.mesh,
	    value = self.valueLeft,
            viewer = 'None')

        self.eq = DiffusionEquation(
            self.var,
            name = "concentration",
            transientCoeff = 0., 
            diffusionCoeff = 1.,
            solver = LinearPCGSolver(
            tolerance = 1.e-15, 
            steps = 1000
            ),
            boundaryConditions=(
            FixedValue(self.mesh.getFacesLeft(),self.valueLeft),
            FixedValue(self.mesh.getFacesRight(),self.valueRight),
            FixedFlux(self.mesh.getFacesTop(),0.),
            FixedFlux(self.mesh.getFacesBottom(),0.)
            )
            )

        self.it = Iterator((self.eq,))

    def assertWithinTolerance(self, first, second, tol = 1e-10, msg=None):
        """Fail if the two objects are unequal by more than tol.
        """
        if abs(first - second) > tol:
            raise self.failureException, (msg or '%s !~ %s' % (first, second))
        
    def testResult(self):
        self.it.iterate(1,1.)
        array = self.var.getArray()
        (lx,ly) = self.mesh.getPhysicalShape()
        vl = self.valueLeft
        vr = self.valueRight

        for cell in self.mesh.getCells():
            coords = cell.getCenter()
            id = cell.getId()
            val = vl + (vr - vl) * coords[0] / lx
            norm = abs(array[id] - val)        
            self.assertWithinTolerance(norm, 0.0, 1e-8,("cell(%g)'s value of %g differs from %g by %g" % (id,array[id],val,norm)))
            
	    
class TestSteadyStateDiffusion20x20(TestSteadyStateDiffusion):
    def setUp(self):
	self.nx = 20
	self.ny = 20
	TestSteadyStateDiffusion.setUp(self)	    
	
class TestSteadyStateDiffusion50x50(TestSteadyStateDiffusion):
    def setUp(self):
	self.nx = 50
	self.ny = 50
	TestSteadyStateDiffusion.setUp(self)	    

class  TestSteadyStateDiffusion1D(TestSteadyStateDiffusion):
    def setUp(self):
	self.nx = 50
	self.ny = 1
	TestSteadyStateDiffusion.setUp(self)	    
	
def suite():
    suite1D = unittest.makeSuite(TestSteadyStateDiffusion1D, 'test')
    suite20 = unittest.makeSuite(TestSteadyStateDiffusion20x20, 'test')
    suite50 = unittest.makeSuite(TestSteadyStateDiffusion50x50, 'test')
    alltests = unittest.TestSuite((suite1D,suite20,suite50))
    return alltests
    
if __name__ == '__main__':
    unittest.main(defaultTest='suite')

            
            
