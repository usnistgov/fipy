#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "advectionEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 9/3/04 {10:37:34 PM} 
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
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

"""

The `AdvectionTerm` object constructs the b vector contribution for
the advection term given by

.. raw:: latex

    $$ u | \\nabla \\phi | $$

from the advection equation given by:

.. raw:: latex

    $$ \\frac{\\partial \\phi}{\\partial t} + u | \\nabla \\phi | = 0$$

The construction of the gradient magnitude term requires upwinding.
The formula used here is given by:

.. raw:: latex

    $$ u_P | \\nabla \\phi |_P = \\max \\left( u_P , 0 \\right) \\left[  \\sum_A \\min \\left( \\frac{ \\phi_A - \\phi_P } { d_{AP}}, 0 \\right)^2 \\right]^{1/2} +  \\min \\left( u_P , 0 \\right) \\left[  \\sum_A \\max \\left( \\frac{ \\phi_A - \\phi_P } { d_{AP}}, 0 \\right)^2 \\right]^{1/2} $$

Here are some simple test cases for this problem:

   >>> from fipy.meshes.grid2D import Grid2D
   >>> mesh = Grid2D(dx = 1., dy = 1., nx = 3, ny = 1) 
   
Trivial test:
   >>> L, b = AdvectionTerm(0., mesh = mesh).buildMatrix(Numeric.zeros(3, 'd'))
   >>> Numeric.allclose(b, Numeric.zeros(3, 'd'), atol = 1e-10)
   1
   
Less trivial test:

   >>> L, b = AdvectionTerm(1., mesh = mesh).buildMatrix(Numeric.arange(3))
   >>> Numeric.allclose(b, Numeric.array((0., -1., -1.)), atol = 1e-10)
   1

Even less trivial

   >>> L, b = AdvectionTerm(-1., mesh = mesh).buildMatrix(Numeric.arange(3))
   >>> Numeric.allclose(b, Numeric.array((1., 1., 0.)), atol = 1e-10)
   1

Another trivial test case (more trivial than a trivial test case
standing on a harpsichord singing 'trivial test cases are here again')

   >>> vel = Numeric.array((-1, 2, -3))
   >>> L, b = AdvectionTerm(vel, mesh = mesh).buildMatrix(Numeric.array((4,6,1)))
   >>> Numeric.allclose(b, -vel * Numeric.array((2, Numeric.sqrt(5**2 + 2**2), 5)), atol = 1e-10)
   1

Somewhat less trivial test case:

   >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 2)
   >>> vel = Numeric.array((3, -5, -6, -3)) 
   >>> L, b = AdvectionTerm(vel, mesh = mesh).buildMatrix(Numeric.array((3 , 1, 6, 7)))
   >>> answer = -vel * Numeric.array((2, Numeric.sqrt(2**2 + 6**2), 1, 0))
   >>> Numeric.allclose(b, answer, atol = 1e-10)
   1
   
"""
__docformat__ = 'restructuredtext'

import Numeric
import MA

from fipy.terms.term import Term
from fipy.tools.sparseMatrix import SparseMatrix

class AdvectionTerm(Term):

    def __init__(self, coeff = None, mesh = None):

        self.coeff = coeff
        Term.__init__(self, weight = None, mesh = mesh)
        
        NCells = self.mesh.getNumberOfCells()
        self.L = SparseMatrix(size = NCells)
        
    def buildMatrix(self, oldArray, coeffScale = None, varScale = None, dt = None):

        NCells = self.mesh.getNumberOfCells()
        NCellFaces = self.mesh.getMaxFacesPerCell()

        cellValues = Numeric.repeat(oldArray[:,Numeric.NewAxis], NCellFaces, axis = 1)
        
        cellIDs = Numeric.repeat(Numeric.arange(NCells)[:,Numeric.NewAxis], NCellFaces, axis = 1)
        cellToCellIDs = self.mesh.getCellToCellIDs()

        cellToCellIDs = MA.where(cellToCellIDs.mask(), cellIDs, cellToCellIDs) 

        adjacentValues = Numeric.take(oldArray, cellToCellIDs)

        differences = self.getDifferences(adjacentValues, cellValues, oldArray, cellToCellIDs)
        differences = differences.filled(fill_value = 0)
        
        minsq = Numeric.sqrt(Numeric.sum(Numeric.minimum(differences, Numeric.zeros((NCells, NCellFaces)))**2, axis = 1))
        maxsq = Numeric.sqrt(Numeric.sum(Numeric.maximum(differences, Numeric.zeros((NCells, NCellFaces)))**2, axis = 1))

        coeff = Numeric.array(self.coeff)
        coeffXdiffereneces = coeff * ((coeff > 0.) * minsq + (coeff < 0.) * maxsq)

        return (self.L, -coeffXdiffereneces * self.mesh.getCellVolumes())
        
    def getDifferences(self, adjacentValues, cellValues, oldArray, cellToCellIDs):
        return (adjacentValues - cellValues) / self.mesh.getCellToCellDistances()

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
        
        

    


        

    

    
      

        
        
        
        
