#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "testVariableDiffusion.py"
 #                                    created: 11/26/03 {3:23:47 PM}
 #                                last update: 11/26/03 {11:11:48 AM} 
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
from equations.diffusionEquation import DiffusionEquation
from solvers.linearPCGSolver import LinearPCGSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator
from variables.variable import Variable
import Numeric

class TestVariableDiffusion(unittest.TestCase):
    """
    Variable diffusion test case. Diffusion
    varies in space. Test is for the steady state profile.
    Test case considers a bar with the following diffusion coefficient:
    
    1    when   0 < x < L/2
    0.1  when   L/2 < x < 3L/4
    1    when   3L/4 < x < L

    where L is the length of the bar. The boundary conditions are a fixed value
    of 0 on the left side and a fixed flux of 1 on the right side. A simple analytical
    answer to this diffusion problem is given in ansFunc(). In order to obtain very
    accurate numerical solutions it is necessary for cell centers to lie at exactly
    L/4 and 3L/4. Thus the test case should give accurate results for meshes of length
    2, 6, 10, 14, 18,...,50,... i.e. 4*n+2 where n is an integer.
    """
    
    def setUp(self):

        self.valueLeft = 0.
        self.fluxRight = 1.

        Lx = 1.2345
        Ly = 2.3456

        dx = Lx / self.nx
        dy = Ly

        self.mesh = Grid2D(dx,dy,self.nx,self.ny)
        
        self.var = Variable(
            name = "diffusion",
            mesh = self.mesh,
	    value = self.valueLeft,
            viewer = 'None')

        def diffFunc(x,L = Lx):
            if x[0] < L / 4.:
                return 1.
            elif x[0] < 3.* L / 4.:
                return 0.1
            else:
                return 1.

        def ansFunc(x,L = Lx):
            if x[0] < L / 4.:
                return x[0]
            elif x[0] < 3. * L / 4.:
                return 10 * x[0] - 9. * L / 4.
            else:
                return x[0] + 18. * L / 4.
            
        self.ansFunc = ansFunc

        faces = self.mesh.getFaces()
        diffCoeff = Numeric.zeros((len(faces)),'d')
        for face in faces:
            diffCoeff[face.getId()] = diffFunc(face.getCenter())        

        self.eq = DiffusionEquation(
            self.var,
            name = "diffusion",
            transientCoeff = 0., 
            diffusionCoeff = diffCoeff,
            solver = LinearPCGSolver(
            tolerance = 1.e-15, 
            steps = 1000
            ),
            boundaryConditions=(
            FixedValue(self.mesh.getFacesLeft(),self.valueLeft),
            FixedFlux(self.mesh.getFacesRight(),self.fluxRight),
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

        for cell in self.mesh.getCells():
            coords = cell.getCenter()
            id = cell.getId()
            val = self.ansFunc(coords)
            norm = abs(array[id] - val)        
            self.assertWithinTolerance(norm, 0.0, 1e-8,("cell(%g)'s value of %g differs from %g by %g" % (id,array[id],val,norm)))
            



class TestVariableDiffusion2x1(TestVariableDiffusion):
    """Variable 1D diffusion on a 1x2 mesh
    """
    def setUp(self):
	self.nx = 2
	self.ny = 1
	TestVariableDiffusion.setUp(self)

	    
class TestVariableDiffusion10x1(TestVariableDiffusion):
    """Variable 1D diffusion on a 1x10 mesh
    """
    def setUp(self):
	self.nx = 10
	self.ny = 1
	TestVariableDiffusion.setUp(self)
	
class TestVariableDiffusion50x1(TestVariableDiffusion):
    """Variable 1D diffusion on a 1x50 mesh
    """
    def setUp(self):
	self.nx = 50
	self.ny = 1
	TestVariableDiffusion.setUp(self)

def suite():
    suite2 = unittest.makeSuite(TestVariableDiffusion2x1, 'test')
    suite10 = unittest.makeSuite(TestVariableDiffusion10x1, 'test')
    suite50 = unittest.makeSuite(TestVariableDiffusion50x1, 'test')
    alltests = unittest.TestSuite((suite2,suite10,suite50))
    return alltests
    
if __name__ == '__main__':
    unittest.main(defaultTest='suite')

            
            
