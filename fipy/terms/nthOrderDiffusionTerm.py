#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "nthOrderDiffusionTerm.py"
 #                                    created: 5/10/04 {11:24:01 AM} 
 #                                last update: 11/30/04 {6:41:08 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  FiPy
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2004-05-10 JEG 1.0 original
 # ###################################################################
 ##


r"""
This term represents a higher order diffusion term. The order of the term is determined
by the number of `coeffs`, such that::

    NthOrderDiffusionTerm(D1, mesh, bcs)
	    
represents a typical 2nd-order diffusion term of the form

.. raw:: latex
        
    $$ \nabla\cdot\left(D_1 \nabla \phi\right) $$
              
and::

    NthOrderDiffusionTerm((D1,D2), mesh, bcs)
            
represents a 4th-order Cahn-Hilliard term of the form
        
.. raw:: latex
        
    $$ \nabla\cdot\left[ D_1 \nabla\cdot\left(D_2 \nabla \phi\right) \right] $$

and so on.

Test, 2nd order, 1 dimension, fixed flux of zero both ends.

   >>> from fipy.meshes.grid2D import Grid2D
   >>> mesh = Grid2D(1., 1., 2, 1)
   >>> term = NthOrderDiffusionTerm(coeffs = (1,))
   >>> coeff = term.getCoeff(mesh)
   >>> print term.getCoefficientMatrix(mesh, coeff)
   -1.000000   1.000000  
    1.000000  -1.000000  
   >>> from fipy.variables.cellVariable import CellVariable
   >>> L,b = term.buildMatrix(var = CellVariable(mesh = mesh))
   >>> print L
   -1.000000   1.000000  
    1.000000  -1.000000  
   >>> print b
   [ 0., 0.,]

Test, 2nd order, 1 dimension, fixed flux 3, fixed value of 4

   >>> from fipy.boundaryConditions.fixedFlux import FixedFlux
   >>> from fipy.boundaryConditions.fixedValue import FixedValue
   >>> bcLeft = FixedFlux(mesh.getFacesLeft(), 3.)
   >>> bcRight = FixedValue(mesh.getFacesRight(), 4.)
   >>> term = NthOrderDiffusionTerm(coeffs = (1.,))
   >>> coeff = term.getCoeff(mesh)
   >>> print term.getCoefficientMatrix(mesh, coeff)
   -1.000000   1.000000  
    1.000000  -1.000000  
   >>> L,b = term.buildMatrix(var = CellVariable(mesh = mesh), boundaryConditions = (bcLeft, bcRight))
   >>> print L
   -1.000000   1.000000  
    1.000000  -3.000000  
   >>> print b
   [-3.,-8.,]

Test, 4th order, 1 dimension, x = 0; fixed flux 3, fixed curvatures 0,
x = 2, fixed value 1, fixed curvature 0

   >>> bcLeft1 = FixedFlux(mesh.getFacesLeft(), 3.)
   >>> from fipy.boundaryConditions.nthOrderBoundaryCondition import NthOrderBoundaryCondition
   >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.getFacesLeft(), 0., 2)
   >>> bcRight1 = FixedValue(mesh.getFacesRight(), 4.)
   >>> bcRight2 =  NthOrderBoundaryCondition(mesh.getFacesRight(), 0., 2)
   >>> term = NthOrderDiffusionTerm(coeffs = (1., 1.))
   >>> coeff = term.getCoeff(mesh)
   >>> print term.getCoefficientMatrix(mesh, coeff)
   -1.000000   1.000000  
    1.000000  -1.000000  
   >>> L,b = term.buildMatrix(var = CellVariable(mesh = mesh), 
   ...                        boundaryConditions = (bcLeft1, bcLeft2, bcRight1, bcRight2))
   >>> print L
    4.000000  -6.000000  
   -4.000000  10.000000  
   >>> print b
   [  1., 21.,]

Test, 4th order, 1 dimension, x = 0; fixed flux 3, fixed curvature 2,
x = 2, fixed value 4, fixed 3rd order -1

   >>> bcLeft1 = FixedFlux(mesh.getFacesLeft(), 3.)
   >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.getFacesLeft(), 2., 2)
   >>> bcRight1 = FixedValue(mesh.getFacesRight(), 4.)
   >>> bcRight2 =  NthOrderBoundaryCondition(mesh.getFacesRight(), -1., 3)
   >>> term = NthOrderDiffusionTerm(coeffs = (-1., 1.))
   >>> coeff = term.getCoeff(mesh)
   >>> print term.getCoefficientMatrix(mesh, coeff)
    1.000000  -1.000000  
   -1.000000   1.000000  
   >>> L,b = term.buildMatrix(var = CellVariable(mesh = mesh),
   ...                        boundaryConditions = (bcLeft1, bcLeft2, bcRight1, bcRight2))
   >>> print L
   -4.000000   6.000000  
    2.000000  -4.000000  
   >>> print b
   [ 3.,-4.,]

Test when dx = 0.5.

   >>> mesh = Grid2D(dx = .5, dy = 1., nx = 2, ny = 1)
   >>> bcLeft1 = FixedValue(mesh.getFacesLeft(), 0.)
   >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.getFacesLeft(), 1., 2)
   >>> bcRight1 = FixedFlux(mesh.getFacesRight(), 1.)
   >>> bcRight2 =  NthOrderBoundaryCondition(mesh.getFacesRight(), 0., 3)
   >>> term = NthOrderDiffusionTerm(coeffs = (1., 1.))
   >>> coeff = term.getCoeff(mesh)
   >>> print term.getCoefficientMatrix(mesh, coeff)
   -2.000000   2.000000  
    2.000000  -2.000000  
   >>> L,b = term.buildMatrix(var = CellVariable(mesh = mesh), 
   ...                        boundaryConditions = (bcLeft1, bcLeft2, bcRight1, bcRight2))
   >>> ans = Numeric.array(((8e+01, -32),(-32, 16)))
   >>> print Numeric.allclose(Numeric.array(L), ans)
   1
   >>> print b
   [-8., 4.,]

Test when dx = 0.25.

   >>> mesh = Grid2D(dx = .25, dy = 1., nx = 2, ny = 1)
   >>> bcLeft1 = FixedValue(mesh.getFacesLeft(), 0.)
   >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.getFacesLeft(), 1., 2)
   >>> bcRight1 = FixedFlux(mesh.getFacesRight(), 1.)
   >>> bcRight2 =  NthOrderBoundaryCondition(mesh.getFacesRight(), 0., 3)
   >>> term = NthOrderDiffusionTerm(coeffs = (1., 1.))
   >>> coeff = term.getCoeff(mesh)
   >>> print term.getCoefficientMatrix(mesh, coeff)
   -4.000000   4.000000  
    4.000000  -4.000000  
   >>> L,b = term.buildMatrix(var = CellVariable(mesh = mesh), 
   ...                        boundaryConditions = (bcLeft1, bcLeft2, bcRight1, bcRight2))
   >>> ans = Numeric.array(((6.4e+2, -2.56e+2), (-2.56e+2, 1.28e+2)))
   >>> print Numeric.allclose(Numeric.array(L), ans)
   1
   >>> print b
   [-24., 16.,]
   
   
"""

__docformat__ = 'restructuredtext'

import Numeric

from fipy.terms.term import Term
from fipy.tools.sparseMatrix import SparseMatrix
from fipy.tools.sparseMatrix import SparseIdentityMatrix
from fipy.variables.addOverFacesVariable import AddOverFacesVariable
import fipy.tools.array as array
import fipy.tools.vector

class NthOrderDiffusionTerm(Term):
    def __init__(self, coeffs):

	self.order = len(coeffs) * 2

	if len(coeffs) > 0:
            self.nthCoeff = coeffs[0]
	else:
	    self.nthCoeff = None

	Term.__init__(self)
	
	if self.order > 0:
	    self.lowerOrderDiffusionTerm = NthOrderDiffusionTerm(coeffs = coeffs[1:])


    def getBoundaryConditions(self, boundaryConditions):
        higherOrderBCs = []
        lowerOrderBCs = []
        for bc in boundaryConditions:
            bcDeriv = bc.getDerivative(self.order - 2)
            if bcDeriv:
                higherOrderBCs.append(bcDeriv)
            else:
                lowerOrderBCs.append(bc)
                
        return higherOrderBCs, lowerOrderBCs
                
    def calcCoeff(self, mesh):
        if self.nthCoeff is not None:
            self.coeff = -self.nthCoeff * mesh.getFaceAreas() / mesh.getCellDistances()
        else:
            self.coeff = None
        
    def getCoefficientMatrix(self, mesh, coeff):
	coefficientMatrix = SparseMatrix(size = mesh.getNumberOfCells(), bandwidth = mesh.getMaxFacesPerCell())

	interiorCoeff = Numeric.array(coeff)
        
	array.put(interiorCoeff, mesh.getExteriorFaceIDs(), 0)
        from fipy.variables.faceVariable import FaceVariable
        interiorCoeff = FaceVariable(mesh = mesh, value = interiorCoeff)
        
	contributions = array.take(interiorCoeff, mesh.getCellFaceIDs())

        interiorCoeff0 = -array.take(Numeric.array(coeff), mesh.getInteriorFaceIDs())
        interiorCoeff1 = -array.take(Numeric.array(coeff), mesh.getInteriorFaceIDs())
        interiorFaceCellIDs = array.take(mesh.getFaceCellIDs(), mesh.getInteriorFaceIDs())

	contributions = array.sum(contributions, 1)	
	coefficientMatrix.addAtDiagonal(contributions)
        
	coefficientMatrix.addAt(interiorCoeff0, interiorFaceCellIDs[:,0], interiorFaceCellIDs[:,1])
	coefficientMatrix.addAt(interiorCoeff1, interiorFaceCellIDs[:,1], interiorFaceCellIDs[:,0])
	
	return coefficientMatrix
	
    def buildMatrix(self, var, boundaryConditions = (), dt = 1.):
        mesh = var.getMesh()
        
        N = mesh.getNumberOfCells()
        volumes = mesh.getCellVolumes()
        if self.order > 0:

            coeff = self.getCoeff(mesh)
            
            coefficientMatrix = self.getCoefficientMatrix(mesh, coeff)

            boundaryB = Numeric.zeros(N,'d')
                
            higherOrderBCs, lowerOrderBCs = self.getBoundaryConditions(boundaryConditions)
            M = mesh.getMaxFacesPerCell()

            coeffs = {
		'cell 1 diag':     self.coeff,
		'cell 1 offdiag': -self.coeff
	    }
	    
	    coeffs['cell 2 offdiag'] = coeffs['cell 1 offdiag']
	    coeffs['cell 2 diag'] = coeffs['cell 1 diag']

            for boundaryCondition in higherOrderBCs:
                LL, bb = boundaryCondition.buildMatrix(N, M, coeffs)
                
                coefficientMatrix += LL
                boundaryB += bb
                
            lowerOrderL, lowerOrderb = self.lowerOrderDiffusionTerm.buildMatrix(var = var, 
                                                                                boundaryConditions = lowerOrderBCs, 
                                                                                dt = dt)
            lowerOrderb = lowerOrderb / volumes
            volMatrix = SparseMatrix(size = N)
            volMatrix.addAtDiagonal(1. / volumes )
            lowerOrderL = volMatrix * lowerOrderL
    
            L = coefficientMatrix * lowerOrderL

            b = coefficientMatrix * lowerOrderb + boundaryB

        else:
            N = mesh.getNumberOfCells()
            L = SparseMatrix(size = N)
            L.addAtDiagonal(volumes)
        
            b = Numeric.zeros((N),'d')
            
        return (L, b)
        
def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
