#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "testExplicitDiffusion.py"
 #                                    created: 11/27/03 {3:23:47 PM}
 #                                last update: 11/27/03 {11:11:48 AM} 
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
from meshes.grid2D import Grid2D
from equations.explicitDiffusionEquation import ExplicitDiffusionEquation
from solvers.linearPCGSolver import LinearPCGSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator
from variables.variable import Variable

class TestExplicitDiffusion(unittest.TestCase):
    """Generic steady-state diffusion class
    Same as TestSteadyStateDiffusion biut for explcit case
    Constructs a mesh, variable, equation, and iterator based
    on the mesh dimensions specified by the child class
    """
    def setUp(self):

        self.valueLeft = 0.
        self.valueRight = 1.

        self.mesh = Grid2D(1.,1.,self.nx,self.ny)
        
        self.var = Variable(
            name = "concentration",
            mesh = self.mesh,
	    value = self.valueLeft,
            viewer = 'None')

        self.eq = ExplicitDiffusionEquation(
            self.var,
            name = "concentration",
            transientCoeff = 1., 
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
        self.it.iterate(10000,0.2)
        array = self.var.getArray()
        (lx,ly) = self.mesh.getPhysicalShape()
        vl = self.valueLeft
        vr = self.valueRight

        for cell in self.mesh.getCells():
            coords = cell.getCenter()
            id = cell.getId()
            val = vl + (vr - vl) * coords[0] / lx
            norm = abs(array[id] - val)
            self.assertWithinTolerance(norm, 0.0, 1e-3,("cell(%g)'s value of %g differs from %g by %g" % (id,array[id],val,norm)))
            
	    
class  TestExplicitDiffusion10(TestExplicitDiffusion):
    """Steady-state 1D diffusion on a 10x1 mesh
    """
    def setUp(self):
	self.nx = 10
	self.ny = 1
	TestExplicitDiffusion.setUp(self)

class  TestExplicitDiffusion50(TestExplicitDiffusion):
    """Steady-state 1D diffusion on a 50x1 mesh
    """
    def setUp(self):
	self.nx = 50
	self.ny = 1
	TestExplicitDiffusion.setUp(self)
	
def suite():
    suite10 = unittest.makeSuite(TestExplicitDiffusion10, 'test')
    suite50 = unittest.makeSuite(TestExplicitDiffusion50, 'test')
    alltests = unittest.TestSuite((suite10,suite50))
    return alltests
    
if __name__ == '__main__':
    unittest.main(defaultTest='suite')

            
            
