#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "advectionEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 9/3/04 {10:29:55 PM} 
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

The `HigherOrderAdvectionTerm` object constructs the `b` vector contribution for
the advection term given by

.. raw:: latex

    $$ u | \\nabla \\phi | $$

from the advection equation given by:

.. raw:: latex

    $$ \\frac{\\partial \\phi}{\\partial t} + u | \\nabla \\phi | = 0$$

The construction of the gradient magnitude term requires upwinding as in the standard
`AdvectionTerm`. The higher order terms are incorperated as follows.
The formula used here is given by:

.. raw:: latex

    $$ u_P | \\nabla \\phi |_P = \\max \\left( u_P , 0 \\right) \\left[  \\sum_A \\min \\left( D_{AP}, 0 \\right)^2 \\right]^{1/2} +  \\min \\left( u_P , 0 \\right) \\left[  \\sum_A \\max \\left( D_{AP}, 0 \\right)^2 \\right]^{1/2} $$

where,

.. raw:: latex

    $$ D_{AP} = \\frac{ \\phi_A - \\phi_P } { d_{AP}} - \\frac{ d_{AP} } {2} m \\left(L_A, L_P \\right) $$

and

.. raw:: latex

    $$ m\\left(x, y\\right) = x \;\; \\text{if} \;\; |x| \le |y| \;\; xy \ge 0 $$
    $$ m\\left(x, y\\right) = y \;\; \\text{if} \;\; |x| > |y|   \;\; xy \ge 0 $$
    $$ m\\left(x, y\\right) = 0 \;\; \\text{if} \;\; xy < 0 $$
    
also,

.. raw:: latex

    $$ L_A = \\frac{\\phi_{AA} + \\phi_P - 2 \\phi_A}{d_{AP}^2} $$
    $$ L_P = \\frac{\\phi_{A} + \\phi_{PP} - 2 \\phi_P}{d_{AP}^2} $$

Here are some simple test cases for this problem:

   >>> from fipy.meshes.grid2D import Grid2D
   >>> mesh = Grid2D(dx = 1., dy = 1., nx = 3, ny = 1) 
   
Trivial test:

   >>> from fipy.variables.cellVariable import CellVariable
   >>> coeff = CellVariable(mesh = mesh, value = Numeric.zeros(3, 'd'))
   >>> L, b = HigherOrderAdvectionTerm(0., mesh = mesh).buildMatrix(coeff)
   >>> Numeric.allclose(b, Numeric.zeros(3, 'd'), atol = 1e-10)
   1
   
Less trivial test:

   >>> coeff = CellVariable(mesh = mesh, value = Numeric.arange(3))
   >>> L, b = HigherOrderAdvectionTerm(1., mesh = mesh).buildMatrix(coeff)
   >>> Numeric.allclose(b, Numeric.array((0., -1., -1.)), atol = 1e-10)
   1

Even less trivial

   >>> coeff = CellVariable(mesh = mesh, value = Numeric.arange(3)) 
   >>> L, b = HigherOrderAdvectionTerm(-1., mesh = mesh).buildMatrix(coeff)
   >>> Numeric.allclose(b, Numeric.array((1., 1., 0.)), atol = 1e-10)
   1

Another trivial test case (more trivial than a trivial test case
standing on a harpsichord singing 'trivial test cases are here again')

   >>> vel = Numeric.array((-1, 2, -3))
   >>> coeff = CellVariable(mesh = mesh, value = Numeric.array((4,6,1))) 
   >>> L, b = HigherOrderAdvectionTerm(vel, mesh = mesh).buildMatrix(coeff)
   >>> Numeric.allclose(b, -vel * Numeric.array((2, Numeric.sqrt(5**2 + 2**2), 5)), atol = 1e-10)
   1

Somewhat less trivial test case:

   >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 2)
   >>> vel = Numeric.array((3, -5, -6, -3))
   >>> coeff = CellVariable(mesh = mesh, value = Numeric.array((3 , 1, 6, 7)))
   >>> L, b = HigherOrderAdvectionTerm(vel, mesh = mesh).buildMatrix(coeff)
   >>> answer = -vel * Numeric.array((2, Numeric.sqrt(2**2 + 6**2), 1, 0))
   >>> Numeric.allclose(b, answer, atol = 1e-10)
   1

For the above test cases the `HigherOrderAdvectionTerm` gives the
same result as the `AdvectionTerm`. The following test imposes a quadratic
field. The higher order term can resolve this field correctly.

.. raw:: latex

    $$ \\phi = x^2 $$

The returned vector `b` should have the value:

.. raw:: latex

    $$ -|\\nabla \\phi| = -\\left|\\frac{\\partial \\phi}{\\partial x}\\right| = - 2 |x| $$

Build the test case in the following way,

   >>> mesh = Grid2D(dx = 1., dy = 1., nx = 5, ny = 1)
   >>> vel = 1.
   >>> coeff = CellVariable(mesh = mesh, value = mesh.getCellCenters()[:,0]**2)
   >>> L, b = AdvectionTerm(vel, mesh = mesh).buildMatrix(coeff)

The first order term is not accurate. The first and last element are ignored because they
don't have any neighbors for higher order evaluation

   >>> Numeric.allclose(b[1:-1], -2 * mesh.getCellCenters()[:,0][1:-1])
   0

The higher order term is spot on.

   >>> L, b = HigherOrderAdvectionTerm(vel, mesh = mesh).buildMatrix(coeff)
   >>> Numeric.allclose(b[1:-1], -2 * mesh.getCellCenters()[:,0][1:-1])
   1

The `HigherOrderAdvectionTerm` will also resolve a circular field with
more accuracy,

.. raw:: latex

    $$ \\phi = \\left( x^2 + y^2 \\right)^{1/2} $$

Build the test case in the following way,

   >>> mesh = Grid2D(dx = 1., dy = 1., nx = 10, ny = 10)
   >>> vel = 1.
   >>> x = mesh.getCellCenters()[:,0]
   >>> y = mesh.getCellCenters()[:,1]
   >>> r = Numeric.sqrt(x**2 + y**2)
   >>> coeff = CellVariable(mesh = mesh, value = r)
   >>> L, b = AdvectionTerm(1., mesh = mesh).buildMatrix(coeff)
   >>> error = Numeric.reshape(Numeric.reshape(b, (10,10))[2:-2,2:-2] + 1, (36,))
   >>> print error[Numeric.argmax(error)]
   0.123105625618

The maximum error is large (about 12%) for the first order advection.

   >>> L, b = HigherOrderAdvectionTerm(1., mesh = mesh).buildMatrix(coeff)
   >>> error = Numeric.reshape(Numeric.reshape(b, (10,10))[2:-2,2:-2] + 1, (36,))
   >>> print error[Numeric.argmax(error)]
   0.0201715476597

The maximum error is 2% when using a higher order contribution.

"""
__docformat__ = 'restructuredtext'

import MA
import Numeric

from advectionTerm import AdvectionTerm
import fipy.tools.array as array

class HigherOrderAdvectionTerm(AdvectionTerm):

    def getDifferences(self, adjacentValues, cellValues, oldArray, cellToCellIDs):
        
        dAP = self.mesh.getCellToCellDistances()
        
##        adjacentGradient = Numeric.take(oldArray.getGrad(), cellToCellIDs)
        from fipy.meshes.numMesh.mesh import MAtake
        adjacentGradient = MAtake(oldArray.getGrad(), self.mesh.getCellToCellIDs())
        adjacentNormalGradient = array.dot(adjacentGradient, self.mesh.getCellNormals(), axis = 2)
        adjacentUpValues = cellValues + 2 * dAP * adjacentNormalGradient

        cellIDs = Numeric.reshape(Numeric.repeat(Numeric.arange(self.mesh.getNumberOfCells()), self.mesh.getMaxFacesPerCell()), cellToCellIDs.shape)
        cellIDs = MA.masked_array(cellIDs, mask = self.mesh.getCellToCellIDs().mask())
        cellGradient = MAtake(oldArray.getGrad(), cellIDs)
        cellNormalGradient = array.dot(cellGradient, self.mesh.getCellNormals(), axis = 2)
        cellUpValues = adjacentValues - 2 * dAP * cellNormalGradient
        
        cellLaplacian = (cellUpValues + adjacentValues - 2 * cellValues) / dAP**2

        adjacentLaplacian = (adjacentUpValues + cellValues - 2 * adjacentValues) / dAP**2
        adjacentLaplacian = adjacentLaplacian.filled(0)
        cellLaplacian = cellLaplacian.filled(0)

##        print 'adjacentLaplacian',adjacentLaplacian
##        print 'cellLaplacian',cellLaplacian

        mm = Numeric.where(cellLaplacian * adjacentLaplacian < 0.,
                           0.,
                           Numeric.where(abs(cellLaplacian) > abs(adjacentLaplacian),
                                         adjacentLaplacian,
                                         cellLaplacian))
        

        return AdvectionTerm.getDifferences(self, adjacentValues, cellValues, oldArray, cellToCellIDs) -  mm * dAP / 2.

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
        
        

    


        

    

    
      

        
        
        
        
