#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "nthOrderDiffusionTerm.py"
 #                                    created: 5/10/04 {11:24:01 AM} 
 #                                last update: 11/19/04 {5:21:41 PM} 
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
   >>> term = NthOrderDiffusionTerm((1,), mesh, ())
   >>> print term.getCoefficientMatrix()
    1.000000  -1.000000  
   -1.000000   1.000000  
   >>> L,b = term.buildMatrix((0,0), 1., 1.)
   >>> print L
    1.000000  -1.000000  
   -1.000000   1.000000  
   >>> print b
   [ 0., 0.,]

Test, 2nd order, 1 dimension, fixed flux 3, fixed value of 4

   >>> from fipy.boundaryConditions.fixedFlux import FixedFlux
   >>> from fipy.boundaryConditions.fixedValue import FixedValue
   >>> bcLeft = FixedFlux(mesh.getFacesLeft(), 3.)
   >>> bcRight = FixedValue(mesh.getFacesRight(), 4.)
   >>> term = NthOrderDiffusionTerm((1.,), mesh, (bcLeft, bcRight))
   >>> print term.getCoefficientMatrix()
    1.000000  -1.000000  
   -1.000000   1.000000  
   >>> L,b = term.buildMatrix((0.,0.), 1., 1.)
   >>> print L
    1.000000  -1.000000  
   -1.000000   3.000000  
   >>> print b
   [ 3., 8.,]

Test, 4th order, 1 dimension, x = 0; fixed flux 3, fixed curvatures 0,
x = 2, fixed value 1, fixed curvature 0

   >>> bcLeft1 = FixedFlux(mesh.getFacesLeft(), 3.)
   >>> from fipy.boundaryConditions.nthOrderBoundaryCondition import NthOrderBoundaryCondition
   >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.getFacesLeft(), 0., 2)
   >>> bcRight1 = FixedValue(mesh.getFacesRight(), 4.)
   >>> bcRight2 =  NthOrderBoundaryCondition(mesh.getFacesRight(), 0., 2)
   >>> term = NthOrderDiffusionTerm((1., 1.), mesh, (bcLeft1, bcLeft2, bcRight1, bcRight2))
   >>> print term.getCoefficientMatrix()
   -1.000000   1.000000  
    1.000000  -1.000000  
   >>> L,b = term.buildMatrix((0.,0.), 1., 1.)
   >>> print L
   -4.000000   6.000000  
    4.000000  -10.000000 
   >>> print b
   [ -1.,-21.,]

Test, 4th order, 1 dimension, x = 0; fixed flux 3, fixed curvature 2,
x = 2, fixed value 4, fixed 3rd order -1

   >>> bcLeft1 = FixedFlux(mesh.getFacesLeft(), 3.)
   >>> from fipy.boundaryConditions.nthOrderBoundaryCondition import NthOrderBoundaryCondition
   >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.getFacesLeft(), 2., 2)
   >>> bcRight1 = FixedValue(mesh.getFacesRight(), 4.)
   >>> bcRight2 =  NthOrderBoundaryCondition(mesh.getFacesRight(), -1., 3)
   >>> term = NthOrderDiffusionTerm((1., 1.), mesh, (bcLeft1, bcLeft2, bcRight1, bcRight2))
   >>> print term.getCoefficientMatrix()
   -1.000000   1.000000  
    1.000000  -1.000000  
   >>> L,b = term.buildMatrix((0.,0.), 1., 1.)
   >>> print L
   -4.000000   6.000000  
    2.000000  -4.000000  
   >>> print b
   [ 3.,-4.,]

Test when dx = 0.5.

   >>> mesh = Grid2D(dx = .5, dy = 1., nx = 2, ny = 1)
   >>> bcLeft1 = FixedValue(mesh.getFacesLeft(), 0.)
   >>> from fipy.boundaryConditions.nthOrderBoundaryCondition import NthOrderBoundaryCondition
   >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.getFacesLeft(), 1., 2)
   >>> bcRight1 = FixedFlux(mesh.getFacesRight(), 1.)
   >>> bcRight2 =  NthOrderBoundaryCondition(mesh.getFacesRight(), 0., 3)
   >>> term = NthOrderDiffusionTerm((1., 1.), mesh, (bcLeft1, bcLeft2, bcRight1, bcRight2))
   >>> print term.getCoefficientMatrix()
   -2.000000   2.000000  
    2.000000  -2.000000  
   >>> L,b = term.buildMatrix((0.,0.), 1., 1.)
   >>> print L
   -8.00e+01  32.000000  
   32.000000  -16.000000 
   >>> print b
   [ 8.,-4.,]

Test when dx = 0.25.

   >>> mesh = Grid2D(dx = .25, dy = 1., nx = 2, ny = 1)
   >>> bcLeft1 = FixedValue(mesh.getFacesLeft(), 0.)
   >>> from fipy.boundaryConditions.nthOrderBoundaryCondition import NthOrderBoundaryCondition
   >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.getFacesLeft(), 1., 2)
   >>> bcRight1 = FixedFlux(mesh.getFacesRight(), 1.)
   >>> bcRight2 =  NthOrderBoundaryCondition(mesh.getFacesRight(), 0., 3)
   >>> term = NthOrderDiffusionTerm((1., 1.), mesh, (bcLeft1, bcLeft2, bcRight1, bcRight2))
   >>> print term.getCoefficientMatrix()
   -4.000000   4.000000  
    4.000000  -4.000000  
   >>> L,b = term.buildMatrix((0.,0.), 1., 1.)
   >>> print L
   -6.40e+02   2.56e+02  
    2.56e+02  -1.28e+02  
   >>> print b
   [ 24.,-16.,]
   
   
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
    def __init__(self, coeffs, mesh, boundaryConditions):

	self.order = len(coeffs) * 2

##        if globalOrder is None:
##            self.globalOrder = self.order
##        else:
##            self.globalOrder = globalOrder
        
	if len(coeffs) > 0:
	    self.coeff = coeffs[0] * mesh.getFaceAreas() / mesh.getCellDistances()
            ## Added to change the sign so that all terms are added to
            ## the right hand side, without this the 4th order term is
            ## added to the left side of the equation.
            if self.order % 4 == 0:
                self.coeff = -self.coeff
	else:
	    self.coeff = None

	self.boundaryConditions = []
        lowerBoundaryConditions = []
        for bc in boundaryConditions:
            bcDeriv = bc.getDerivative(self.order - 2)
	    if bcDeriv:
		self.boundaryConditions.append(bcDeriv)
	    else:
		lowerBoundaryConditions.append(bc)
            
	Term.__init__(self, weight = None, mesh = mesh)
	
	N = mesh.getNumberOfCells()
	if self.order > 0:
	    self.lowerOrderDiffusionTerm = NthOrderDiffusionTerm(coeffs = coeffs[1:], mesh = mesh, boundaryConditions = lowerBoundaryConditions)

    def getCoefficientMatrix(self):
	mesh = self.getMesh()
	
	coefficientMatrix = SparseMatrix(size = mesh.getNumberOfCells(), bandwidth = mesh.getMaxFacesPerCell())

	interiorCoeff = Numeric.array(self.coeff)
        
	array.put(interiorCoeff, mesh.getExteriorFaceIDs(), 0)
        from fipy.variables.faceVariable import FaceVariable
        interiorCoeff = FaceVariable(mesh = mesh, value = interiorCoeff)
        
	contributions = array.take(interiorCoeff, mesh.getCellFaceIDs())

        ## divide through by the volume for each cell

        interiorCoeff0 = -array.take(Numeric.array(self.coeff), mesh.getInteriorFaceIDs())
        interiorCoeff1 = -array.take(Numeric.array(self.coeff), mesh.getInteriorFaceIDs())
        interiorFaceCellIDs = array.take(mesh.getFaceCellIDs(), mesh.getInteriorFaceIDs())

##        if self.globalOrder != self.order:
##            vols = Numeric.array(mesh.getCellVolumes())
##            contributions = contributions / vols[:,Numeric.NewAxis]
##            interiorCoeff0 = interiorCoeff0 / Numeric.take(vols, interiorFaceCellIDs[:,0])
##            interiorCoeff1 = interiorCoeff1 / Numeric.take(vols, interiorFaceCellIDs[:,1])

	contributions = array.sum(contributions, 1)	
	coefficientMatrix.addAtDiagonal(contributions)
        
##	interiorCoeff1 = -array.take(Numeric.array(self.coeff), mesh.getInteriorFaceIDs())
##	interiorFaceCellIDs = array.take(mesh.getFaceCellIDs(), mesh.getInteriorFaceIDs())
	
	coefficientMatrix.addAt(interiorCoeff0, interiorFaceCellIDs[:,0], interiorFaceCellIDs[:,1])
	coefficientMatrix.addAt(interiorCoeff1, interiorFaceCellIDs[:,1], interiorFaceCellIDs[:,0])
	
	return coefficientMatrix
	
    def buildMatrix(self, oldArray, coeffScale, varScale, dt = 1.):
        N = self.getMesh().getNumberOfCells()
        volumes = self.mesh.getCellVolumes()
        if self.order > 0:

	    coefficientMatrix = self.getCoefficientMatrix()
	    boundaryB = 0
	    
##             boundaryB = Numeric.zeros(N,'d')
                
	    M = self.getMesh().getMaxFacesPerCell()
	    
            for boundaryCondition in self.boundaryConditions:
		LL, bb = boundaryCondition.buildMatrix(N, M, self.coeff,-self.coeff, coeffScale)
		
		coefficientMatrix += LL
		boundaryB += bb / varScale
            
            lowerOrderL, lowerOrderb = self.lowerOrderDiffusionTerm.buildMatrix(oldArray, coeffScale, varScale)

            lowerOrderb = lowerOrderb / volumes
            volMatrix = SparseMatrix(size = N)
            volMatrix.addAtDiagonal(1. / volumes )
            lowerOrderL = volMatrix * lowerOrderL
    
            L = coefficientMatrix * lowerOrderL

            iseven = not ((self.order / 2) % 2)
            if iseven:
                boundaryB = -boundaryB

##            print "order",self.order
##            print "coefficientMatrix",coefficientMatrix
##            print "L",L
##            print "lowerOrderL",lowerOrderL

            b = coefficientMatrix * lowerOrderb + boundaryB

##            print "order",self.order
##            print "coefficientMatrix",coefficientMatrix
##            print "L",L
##            print "lowerOrderL",lowerOrderL

##            if self.globalOrder != self.order:
##                volumes = self.mesh.getCellVolumes()
##                b = b / volumes
##                L = SparseMatrix(size = N).addAtDiagonal(1. / volumes ) * L
    
        else:
            N = self.getMesh().getNumberOfCells()
            L = SparseMatrix(size = N)
            L.addAtDiagonal(volumes)
        
##            L = SparseIdentityMatrix(size = N)
            b = Numeric.zeros((N),'d')
            
##        for boundaryCondition in self.boundaryConditions:
##            LL,bb,ids = boundaryCondition.getContribution(self.coeff,-self.coeff)
	    
##            self.L.addAt(LL / coeffScale,ids,ids)
##            fipy.tools.vector.putAdd(self.b, ids, bb/(coeffScale * varScale))

        return (L, b)
        
def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
