#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "diffusionTerm.py"
 #                                    created: 11/13/03 {11:39:03 AM} 
 #                                last update: 1/17/06 {12:00:58 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
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
 #  2003-11-13 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import Numeric

from fipy.terms.term import Term
from fipy.tools.sparseMatrix import _SparseMatrix
from fipy.tools import numerix

class DiffusionTerm(Term):

    r"""
    This term represents a higher order diffusion term. The order of the term is determined
    by the number of `coeffs`, such that::

        DiffusionTerm(D1, mesh, bcs)

    represents a typical 2nd-order diffusion term of the form

    .. raw:: latex

        $$ \nabla\cdot\left(D_1 \nabla \phi\right) $$

    and::

        DiffusionTerm((D1,D2), mesh, bcs)

    represents a 4th-order Cahn-Hilliard term of the form

    .. raw:: latex

        $$ \nabla \cdot \left\{ D_1 \nabla \left[ \nabla\cdot\left( D_2 \nabla \phi\right) \right] \right\} $$

    and so on.

    """

    def __init__(self, coeff = (1.,)):
        """
        Create a `DiffusionTerm`.

        :Parameters:
          - `coeff`: `Tuple` or `list` of `FaceVariables` or numbers.
	  
        """
        if type(coeff) not in (type(()), type([])):
            coeff = (coeff,)

        self.order = len(coeff) * 2

        if len(coeff) > 0:
            self.nthCoeff = coeff[0]

            from fipy.variables.variable import Variable
            if not isinstance(self.nthCoeff, Variable):
                self.nthCoeff = Variable(value = self.nthCoeff)

	    from fipy.variables.cellVariable import CellVariable
	    from fipy.variables.vectorCellVariable import VectorCellVariable
	    if isinstance(self.nthCoeff, VectorCellVariable) or isinstance(self.nthCoeff, CellVariable):
		self.nthCoeff = self.nthCoeff.getArithmeticFaceValue()

	    from fipy.variables.vectorFaceVariable import VectorFaceVariable
	    if isinstance(self.nthCoeff, VectorFaceVariable):
		self.nthCoeff = self.nthCoeff.dot(self.nthCoeff.getMesh()._getFaceNormals()**2)

	    from fipy.variables.faceVariable import FaceVariable
	    if not isinstance(self.nthCoeff, FaceVariable):
		if self.nthCoeff.getShape() != ():
		    raise TypeError, "The coefficient must be a FaceVariable, CellVariable, VectorFaceVariable, VectorCellVariable, or a scalar value."

        else:
            self.nthCoeff = None

        Term.__init__(self, coeff = coeff)
        
        if self.order > 0:
            self.lowerOrderDiffusionTerm = DiffusionTerm(coeff = coeff[1:])
        
    def __neg__(self):
        """
        Negate the term.

           >>> -DiffusionTerm(coeff = [1.])
           DiffusionTerm(coeff = [-1.0])

           >>> -DiffusionTerm()
           DiffusionTerm(coeff = [-1.0])
           
        """
        negatedCoeff = list(self.coeff)
        negatedCoeff[0] = -negatedCoeff[0]
        return self.__class__(coeff = negatedCoeff)
            
    def _getBoundaryConditions(self, boundaryConditions):
        higherOrderBCs = []
        lowerOrderBCs = []

        for bc in boundaryConditions:
            bcDeriv = bc._getDerivative(self.order - 2)
            if bcDeriv:
                higherOrderBCs.append(bcDeriv)
            else:
                lowerOrderBCs.append(bc)
                
        return higherOrderBCs, lowerOrderBCs
                
    def _calcGeomCoeff(self, mesh):
        if self.nthCoeff is not None:
            return self.nthCoeff * mesh._getFaceAreas() / mesh._getCellDistances()
        else:
            return None
        
    def _getCoefficientMatrix(self, mesh, coeff):
        coefficientMatrix = _SparseMatrix(size = mesh.getNumberOfCells(), bandwidth = mesh._getMaxFacesPerCell())

        interiorCoeff = Numeric.array(coeff)
        
        numerix.put(interiorCoeff, mesh.getExteriorFaceIDs(), 0)
        from fipy.variables.faceVariable import FaceVariable
        interiorCoeff = FaceVariable(mesh = mesh, value = interiorCoeff)
        
        contributions = numerix.take(interiorCoeff, mesh._getCellFaceIDs())

        interiorCoeff0 = -numerix.take(Numeric.array(coeff), mesh.getInteriorFaceIDs())
        interiorCoeff1 = -numerix.take(Numeric.array(coeff), mesh.getInteriorFaceIDs())
        interiorFaceCellIDs = numerix.take(mesh.getFaceCellIDs(), mesh.getInteriorFaceIDs())

        contributions = numerix.sum(contributions, 1)	
        coefficientMatrix.addAtDiagonal(contributions)
        
        coefficientMatrix.addAt(interiorCoeff0, interiorFaceCellIDs[:,0], interiorFaceCellIDs[:,1])
        coefficientMatrix.addAt(interiorCoeff1, interiorFaceCellIDs[:,1], interiorFaceCellIDs[:,0])
        
        return coefficientMatrix
        
    def _bcAdd(self, coefficientMatrix, boundaryB, LLbb):
        coefficientMatrix += LLbb[0]
        boundaryB += LLbb[1]
        
    def _doBCs(self, higherOrderBCs, N, M, coeffs, coefficientMatrix, boundaryB):
        [self._bcAdd(coefficientMatrix, boundaryB, boundaryCondition._buildMatrix(N, M, coeffs)) for boundaryCondition in higherOrderBCs]
            
        return coefficientMatrix, boundaryB

##    def _buildMatrix(self, var, boundaryConditions = (), dt = 1., coefficientMatrix = None):
    def _buildMatrix(self, var, boundaryConditions = (), dt = 1.):
        mesh = var.getMesh()
        
        N = mesh.getNumberOfCells()
        volumes = mesh.getCellVolumes()
        if self.order > 0:

            coeff = self._getGeomCoeff(mesh)
            
##            if coefficientMatrix is None:
            coefficientMatrix = self._getCoefficientMatrix(mesh, -coeff)

            boundaryB = Numeric.zeros(N,'d')
                
            higherOrderBCs, lowerOrderBCs = self._getBoundaryConditions(boundaryConditions)
            M = mesh._getMaxFacesPerCell()

            coeffs = {
                'cell 1 diag':    -coeff,
                'cell 1 offdiag':  coeff
            }

            coeffs['cell 2 offdiag'] = coeffs['cell 1 offdiag']
            coeffs['cell 2 diag'] = coeffs['cell 1 diag']

            coefficientMatrix, boundaryB = self._doBCs(higherOrderBCs, N, M, coeffs, coefficientMatrix, boundaryB)
##             for boundaryCondition in higherOrderBCs:
##                 LL, bb = boundaryCondition._buildMatrix(N, M, coeffs)
##                 
##                 coefficientMatrix += LL
##                 boundaryB += bb

            lowerOrderL, lowerOrderb = self.lowerOrderDiffusionTerm._buildMatrix(var = var, 
                                                                                 boundaryConditions = lowerOrderBCs, 
                                                                                 dt = dt)
##                                                                              coefficientMatrix = coefficientMatrix)

            lowerOrderb = lowerOrderb / volumes
            volMatrix = _SparseMatrix(size = N)
            volMatrix.addAtDiagonal(1. / volumes )
            lowerOrderL = volMatrix * lowerOrderL
            
            L = coefficientMatrix * lowerOrderL

            b = coefficientMatrix * lowerOrderb + boundaryB

        else:
            N = mesh.getNumberOfCells()
            L = _SparseMatrix(size = N)
            L.addAtDiagonal(volumes)
        
            b = Numeric.zeros((N),'d')
            
        return (L, b)
        
    def _test(self):
        r"""
        Test, 2nd order, 1 dimension, fixed flux of zero both ends.

           >>> from fipy.meshes.grid1D import Grid1D
           >>> mesh = Grid1D(dx = 1., nx = 2)
           >>> term = DiffusionTerm(coeff = (1,))
           >>> coeff = term._getGeomCoeff(mesh)
           >>> print term._getCoefficientMatrix(mesh, coeff)
            1.000000  -1.000000  
           -1.000000   1.000000  
           >>> from fipy.variables.cellVariable import CellVariable
           >>> L,b = term._buildMatrix(var = CellVariable(mesh = mesh))
           >>> print L
           -1.000000   1.000000  
            1.000000  -1.000000  
           >>> print b
           [ 0., 0.,]
           
        The coefficient must be a `FaceVariable`, a `CellVariable` (which will
        be interpolated to a `FaceVariable`), or a scalar value 
        
           >>> from fipy.variables.faceVariable import FaceVariable
           >>> term = DiffusionTerm(coeff = FaceVariable(mesh = mesh, value = 1))
           >>> coeff = term._getGeomCoeff(mesh)
           >>> print term._getCoefficientMatrix(mesh, coeff)
            1.000000  -1.000000  
           -1.000000   1.000000  
           >>> L,b = term._buildMatrix(var = CellVariable(mesh = mesh))
           >>> print L
           -1.000000   1.000000  
            1.000000  -1.000000  
           >>> print b
           [ 0., 0.,]

           >>> term = DiffusionTerm(coeff = CellVariable(mesh = mesh, value = 1))
           >>> coeff = term._getGeomCoeff(mesh)
           >>> print term._getCoefficientMatrix(mesh, coeff)
            1.000000  -1.000000  
           -1.000000   1.000000  
           >>> L,b = term._buildMatrix(var = CellVariable(mesh = mesh))
           >>> print L
           -1.000000   1.000000  
            1.000000  -1.000000  
           >>> print b
           [ 0., 0.,]

           >>> from fipy.variables.variable import Variable
           >>> term = DiffusionTerm(coeff = Variable(value = 1))
           >>> coeff = term._getGeomCoeff(mesh)
           >>> print term._getCoefficientMatrix(mesh, coeff)
            1.000000  -1.000000  
           -1.000000   1.000000  
           >>> L,b = term._buildMatrix(var = CellVariable(mesh = mesh))
           >>> print L
           -1.000000   1.000000  
            1.000000  -1.000000  
           >>> print b
           [ 0., 0.,]
           
           >>> term = DiffusionTerm(coeff = ((1,2),))
           Traceback (most recent call last):
               ...
           TypeError: The coefficient must be a FaceVariable, CellVariable, VectorFaceVariable, VectorCellVariable, or a scalar value.
           >>> from fipy.variables.vectorFaceVariable import VectorFaceVariable
           >>> term = DiffusionTerm(coeff = VectorFaceVariable(mesh = mesh, value = (1,)))
           >>> from fipy.variables.vectorCellVariable import VectorCellVariable
           >>> term = DiffusionTerm(coeff = VectorCellVariable(mesh = mesh, value = (1,)))

        Test, 2nd order, 1 dimension, fixed flux 3, fixed value of 4

           >>> from fipy.boundaryConditions.fixedFlux import FixedFlux
           >>> from fipy.boundaryConditions.fixedValue import FixedValue
           >>> bcLeft = FixedFlux(mesh.getFacesLeft(), 3.)
           >>> bcRight = FixedValue(mesh.getFacesRight(), 4.)
           >>> term = DiffusionTerm(coeff = (1.,))
           >>> coeff = term._getGeomCoeff(mesh)
           >>> print term._getCoefficientMatrix(mesh, coeff)
            1.000000  -1.000000  
           -1.000000   1.000000  
           >>> L,b = term._buildMatrix(var = CellVariable(mesh = mesh), 
           ...                         boundaryConditions = (bcLeft, bcRight))
           >>> print L
           -1.000000   1.000000  
            1.000000  -3.000000  
           >>> print b
           [-3.,-8.,]

        Test, 4th order, 1 dimension, x = 0; fixed flux 3, fixed curvatures 0,
        x = 2, fixed value 1, fixed curvature 0

           >>> bcLeft1 = FixedFlux(mesh.getFacesLeft(), 3.)
           >>> from fipy.boundaryConditions.nthOrderBoundaryCondition \
           ...     import NthOrderBoundaryCondition
           >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.getFacesLeft(), 0., 2)
           >>> bcRight1 = FixedValue(mesh.getFacesRight(), 4.)
           >>> bcRight2 =  NthOrderBoundaryCondition(mesh.getFacesRight(), 0., 2)
           >>> term = DiffusionTerm(coeff = (1., 1.))
           >>> coeff = term._getGeomCoeff(mesh)
           >>> print term._getCoefficientMatrix(mesh, coeff)
            1.000000  -1.000000  
           -1.000000   1.000000  
           >>> L,b = term._buildMatrix(var = CellVariable(mesh = mesh), 
           ...                         boundaryConditions = (bcLeft1, bcLeft2, 
           ...                                               bcRight1, bcRight2))
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
           >>> term = DiffusionTerm(coeff = (-1., 1.))
           >>> coeff = term._getGeomCoeff(mesh)
           >>> print term._getCoefficientMatrix(mesh, coeff)
           -1.000000   1.000000  
            1.000000  -1.000000  
           >>> L,b = term._buildMatrix(var = CellVariable(mesh = mesh),
           ...                         boundaryConditions = (bcLeft1, bcLeft2, 
           ...                                               bcRight1, bcRight2))
           >>> print L
           -4.000000   6.000000  
            2.000000  -4.000000  
           >>> print b
           [ 3.,-4.,]

        Test when dx = 0.5.

           >>> mesh = Grid1D(dx = .5, nx = 2)
           >>> bcLeft1 = FixedValue(mesh.getFacesLeft(), 0.)
           >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.getFacesLeft(), 1., 2)
           >>> bcRight1 = FixedFlux(mesh.getFacesRight(), 1.)
           >>> bcRight2 =  NthOrderBoundaryCondition(mesh.getFacesRight(), 0., 3)
           >>> term = DiffusionTerm(coeff = (1., 1.))
           >>> coeff = term._getGeomCoeff(mesh)
           >>> print term._getCoefficientMatrix(mesh, coeff)
            2.000000  -2.000000  
           -2.000000   2.000000  
           >>> L,b = term._buildMatrix(var = CellVariable(mesh = mesh), 
           ...                         boundaryConditions = (bcLeft1, bcLeft2, 
           ...                                               bcRight1, bcRight2))
           >>> ans = Numeric.array(((8e+01, -32),(-32, 16)))
           >>> print Numeric.allclose(Numeric.array(L), ans)
           1
           >>> print b
           [-8., 4.,]

        Test when dx = 0.25.

           >>> mesh = Grid1D(dx = .25, nx = 2)
           >>> bcLeft1 = FixedValue(mesh.getFacesLeft(), 0.)
           >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.getFacesLeft(), 1., 2)
           >>> bcRight1 = FixedFlux(mesh.getFacesRight(), 1.)
           >>> bcRight2 =  NthOrderBoundaryCondition(mesh.getFacesRight(), 0., 3)
           >>> term = DiffusionTerm(coeff = (1., 1.))
           >>> coeff = term._getGeomCoeff(mesh)
           >>> print term._getCoefficientMatrix(mesh, coeff)
            4.000000  -4.000000  
           -4.000000   4.000000  
           >>> L,b = term._buildMatrix(var = CellVariable(mesh = mesh), 
           ...                         boundaryConditions = (bcLeft1, bcLeft2, 
           ...                                               bcRight1, bcRight2))
           >>> ans = Numeric.array(((6.4e+2, -2.56e+2), (-2.56e+2, 1.28e+2)))
           >>> print Numeric.allclose(Numeric.array(L), ans)
           1
           >>> print b
           [-24., 16.,]

        The following tests are to check that DiffusionTerm can take any of the four
        main Variable types.

	   >>> from fipy.meshes.tri2D import Tri2D
	   >>> mesh = Tri2D(nx = 1, ny = 1)
	   >>> print DiffusionTerm(CellVariable(value = 1, mesh = mesh)).nthCoeff
	   [ 1., 1., 1., 1., 1., 1., 1., 1.,]
	   >>> print DiffusionTerm(FaceVariable(value = 1, mesh = mesh)).nthCoeff
	   [ 1., 1., 1., 1., 1., 1., 1., 1.,]
	   >>> print DiffusionTerm(VectorCellVariable(value = (0.5,1), mesh = mesh)).nthCoeff
	   [ 1.  , 1.  , 0.5 , 0.5 , 0.75, 0.75, 0.75, 0.75,]
	   >>> print DiffusionTerm(VectorFaceVariable(value = (0.5, 1), mesh = mesh)).nthCoeff
	   [ 1.  , 1.  , 0.5 , 0.5 , 0.75, 0.75, 0.75, 0.75,]
	   >>> mesh = Tri2D(nx = 1, ny = 1, dy = 0.1)
	   >>> val = DiffusionTerm(VectorFaceVariable(value = (0.5, 1), mesh = mesh)).nthCoeff
	   >>> print numerix.allclose((1, 1, 0.5, 0.5, 0.9950495, 0.9950495, 0.9950495, 0.9950495), val)
	   1

        """
        pass

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
