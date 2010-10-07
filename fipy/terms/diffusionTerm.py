#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "diffusionTerm.py"
 #
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
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.terms.term import Term
from fipy.tools import numerix

class DiffusionTerm(Term):

    r"""
    This term represents a higher order diffusion term. The order of the term is determined
    by the number of `coeffs`, such that::

        DiffusionTerm(D1)

    represents a typical 2nd-order diffusion term of the form

    .. math::

       \nabla\cdot\left(D_1 \nabla \phi\right)

    and::

        DiffusionTerm((D1,D2))

    represents a 4th-order Cahn-Hilliard term of the form

    .. math::

       \nabla \cdot \left\{ D_1 \nabla \left[ \nabla\cdot\left( D_2 \nabla \phi\right) \right] \right\}

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
            if isinstance(self.nthCoeff, CellVariable):
                self.nthCoeff = self.nthCoeff.getArithmeticFaceValue()

        else:
            self.nthCoeff = None

        Term.__init__(self, coeff = coeff)
        
        if self.order > 0:
            self.lowerOrderDiffusionTerm = DiffusionTerm(coeff = coeff[1:])
        
    def __neg__(self):
        """
        Negate the term.

        >>> -DiffusionTerm(coeff=[1.])
        DiffusionTerm(coeff=[-1.0])

        >>> -DiffusionTerm()
        DiffusionTerm(coeff=[-1.0])
           
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

    def _getNormals(self, mesh):
        return mesh._getFaceCellToCellNormals()

    def _getRotationTensor(self, mesh):
        if not hasattr(self, 'rotationTensor'):

            from fipy.variables.faceVariable import FaceVariable
            rotationTensor = FaceVariable(mesh=mesh, rank=2)
            
            rotationTensor[:, 0] = self._getNormals(mesh)

            if mesh.getDim() == 2:
                rotationTensor[:,1] = rotationTensor[:,0].dot((((0, 1), (-1, 0))))
            elif mesh.getDim() ==3:
                epsilon = 1e-20

                div = numerix.sqrt(1 - rotationTensor[2,0]**2)
                flag = numerix.resize(div > epsilon, (mesh.getDim(), mesh._getNumberOfFaces()))

                rotationTensor[0, 1] = 1
                rotationTensor[:, 1] = numerix.where(flag,
                                                     rotationTensor[:,0].dot((((0, 1, 0), (-1, 0, 0), (0, 0, 0)))) / div,
                                                     rotationTensor[:, 1])

                rotationTensor[1, 2] = 1
                rotationTensor[:, 2] = numerix.where(flag,
                                                     rotationTensor[:,0] * rotationTensor[2,0] / div,
                                                     rotationTensor[:, 2])
                rotationTensor[2, 2] = -div
                
            self.rotationTensor = rotationTensor

        return self.rotationTensor
    
    def _treatMeshAsOrthogonal(self, mesh):
        return mesh._isOrthogonal()

    def _calcAnisotropySource(self, coeff, mesh, var):

        if not hasattr(self, 'anisotropySource'):
            if len(coeff) > 1:
                gradients = var.getGrad().getHarmonicFaceValue().dot(self._getRotationTensor(mesh))
                from fipy.variables.addOverFacesVariable import _AddOverFacesVariable
                self.anisotropySource = _AddOverFacesVariable(gradients[1:].dot(coeff[1:])) * mesh.getCellVolumes()

    def _calcGeomCoeff(self, mesh):
        if self.nthCoeff is not None:
          
            coeff = self.nthCoeff
            shape = numerix.getShape(coeff)

            from fipy.variables.faceVariable import FaceVariable
            if isinstance(coeff, FaceVariable):
                rank = coeff.getRank()
            else:
                rank = len(shape)

            if rank == 0 and self._treatMeshAsOrthogonal(mesh):
                tmpBop = (coeff * mesh._getFaceAreas() / mesh._getCellDistances())[numerix.newaxis, :]

            else:

                if rank == 1 or rank == 0:
                    coeff = coeff * numerix.identity(mesh.getDim())

                if rank > 0:
                    shape = numerix.getShape(coeff)
                    if mesh.getDim() != shape[0] or mesh.getDim() != shape[1]:
                        raise IndexError, 'diffusion coefficent tensor is not an appropriate shape for this mesh'          

                faceNormals = FaceVariable(mesh=mesh, rank=1, value=mesh._getFaceNormals())
                rotationTensor = self._getRotationTensor(mesh)
                rotationTensor[:,0] = rotationTensor[:,0] / mesh._getCellDistances()
                
                tmpBop = faceNormals.dot(coeff).dot(rotationTensor) * mesh._getFaceAreas()

            return tmpBop

        else:

            return None

    def _getCoefficientMatrix(self, SparseMatrix, mesh, coeff):
        interiorCoeff = numerix.array(coeff)
        
        interiorCoeff[mesh.getExteriorFaces().getValue()] = 0
        
        interiorCoeff = numerix.take(interiorCoeff, mesh._getCellFaceIDs())

        coefficientMatrix = SparseMatrix(mesh=mesh, bandwidth = mesh._getMaxFacesPerCell() + 1)
        coefficientMatrix.addAtDiagonal(numerix.sum(interiorCoeff, 0))
        del interiorCoeff

        interiorFaces = mesh.getInteriorFaceIDs()
        interiorFaceCellIDs = mesh.getInteriorFaceCellIDs()

        interiorCoeff = -numerix.take(coeff, interiorFaces, axis=-1)
        coefficientMatrix.addAt(interiorCoeff, interiorFaceCellIDs[0], interiorFaceCellIDs[1])
        interiorCoeff = -numerix.take(coeff, interiorFaces, axis=-1)
        coefficientMatrix.addAt(interiorCoeff, interiorFaceCellIDs[1], interiorFaceCellIDs[0])
        
        return coefficientMatrix
        
    def _bcAdd(self, coefficientMatrix, boundaryB, LLbb):
        coefficientMatrix += LLbb[0]
        boundaryB += LLbb[1]
        
    def _doBCs(self, SparseMatrix, higherOrderBCs, N, M, coeffs, coefficientMatrix, boundaryB):
        for boundaryCondition in higherOrderBCs:
            self._bcAdd(coefficientMatrix, boundaryB, boundaryCondition._buildMatrix(SparseMatrix, N, M, coeffs))
            
        return coefficientMatrix, boundaryB

    def __add__(self, other):
        if isinstance(other, DiffusionTerm):
            from fipy.terms.collectedDiffusionTerm import _CollectedDiffusionTerm
            if isinstance(other, _CollectedDiffusionTerm):
                return other + self
            elif other.order == self.order and self.order <= 2:
                if self.order == 0:
                    return self
                elif self.order == 2:
                    return self.__class__(coeff=self.coeff[0] + other.coeff[0])
            else:
                term = _CollectedDiffusionTerm()
                term += self
                term += other
                return term
        else:
            return Term.__add__(self, other)

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions = (), dt = 1., equation=None):
        mesh = var.getMesh()
        
        N = mesh.getNumberOfCells()
        M = mesh._getMaxFacesPerCell()

        if self.order > 2:

            higherOrderBCs, lowerOrderBCs = self._getBoundaryConditions(boundaryConditions)
            
            lowerOrderL, lowerOrderb = self.lowerOrderDiffusionTerm._buildMatrix(var = var, SparseMatrix=SparseMatrix,
                                                                                 boundaryConditions = lowerOrderBCs, 
                                                                                 dt = dt,
                                                                                 equation=equation)
            del lowerOrderBCs
            
            lowerOrderb = lowerOrderb / mesh.getCellVolumes()
            volMatrix = SparseMatrix(mesh=var.getMesh(), bandwidth = 1)
            
            volMatrix.addAtDiagonal(1. / mesh.getCellVolumes() )
            lowerOrderL = volMatrix * lowerOrderL
            del volMatrix

            if not hasattr(self, 'coeffDict'):

                coeff = self._getGeomCoeff(mesh)[0]
                minusCoeff = -coeff
                
                coeff.dontCacheMe()
                minusCoeff.dontCacheMe()

                self.coeffDict = {
                    'cell 1 diag':     minusCoeff,
                    'cell 1 offdiag':  coeff
                    }
                del coeff
                del minusCoeff

                self.coeffDict['cell 2 offdiag'] = self.coeffDict['cell 1 offdiag']
                self.coeffDict['cell 2 diag'] = self.coeffDict['cell 1 diag']


            mm = self._getCoefficientMatrix(SparseMatrix, mesh, self.coeffDict['cell 1 diag'])
            L, b = self._doBCs(SparseMatrix, higherOrderBCs, N, M, self.coeffDict, 
                               mm, numerix.zeros(N,'d'))
                               
            del higherOrderBCs
            del mm

            b = L * lowerOrderb + b
            del lowerOrderb

            L = L * lowerOrderL
            del lowerOrderL

        elif self.order == 2:

            if not hasattr(self, 'coeffDict'):

                coeff = self._getGeomCoeff(mesh)
                minusCoeff = -coeff[0]

                coeff[0].dontCacheMe()
                minusCoeff.dontCacheMe()

                self.coeffDict = {
                    'cell 1 diag':    minusCoeff,
                    'cell 1 offdiag':  coeff[0]
                    }

                self.coeffDict['cell 2 offdiag'] = self.coeffDict['cell 1 offdiag']
                self.coeffDict['cell 2 diag'] = self.coeffDict['cell 1 diag']

                self._calcAnisotropySource(coeff, mesh, var)

                del coeff
                del minusCoeff
                
            higherOrderBCs, lowerOrderBCs = self._getBoundaryConditions(boundaryConditions)
            del lowerOrderBCs

            L, b = self._doBCs(SparseMatrix, higherOrderBCs, N, M, self.coeffDict, 
                               self._getCoefficientMatrix(SparseMatrix, mesh, self.coeffDict['cell 1 diag']), numerix.zeros(N,'d'))

            if hasattr(self, 'anisotropySource'):
                b -= self.anisotropySource
                               
            del higherOrderBCs


        else:
            
            L = SparseMatrix(mesh=mesh)
            L.addAtDiagonal(mesh.getCellVolumes())
            b = numerix.zeros((N),'d')
            
        return (L, b)
        
    def _test(self):
        r"""
        Test, 2nd order, 1 dimension, fixed flux of zero both ends.

        >>> from fipy.meshes.grid1D import Grid1D
        >>> from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix as SparseMatrix
        >>> from fipy.tools import parallel
        >>> procID = parallel.procID
        >>> mesh = Grid1D(dx = 1., nx = 2)
        >>> term = DiffusionTerm(coeff = (1,))
        >>> coeff = term._getGeomCoeff(mesh)
        >>> M = term._getCoefficientMatrix(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.getNumpyArray(), 
        ...                        (( 1., -1.), 
        ...                         (-1.,  1.))) or procID != 0
        True
        >>> from fipy.variables.cellVariable import CellVariable
        >>> L,b = term._buildMatrix(var=CellVariable(mesh=mesh), SparseMatrix=SparseMatrix)
        >>> print numerix.allclose(L.getNumpyArray(), 
        ...                        ((-1.,  1.), 
        ...                         ( 1., -1.))) or procID != 0
        True
        >>> print numerix.allclose(b, (0., 0.)) or procID != 0
        True

        The coefficient must be a `FaceVariable`, a `CellVariable` (which will
        be interpolated to a `FaceVariable`), or a scalar value 
        
        >>> from fipy.variables.faceVariable import FaceVariable
        >>> term = DiffusionTerm(coeff=FaceVariable(mesh=mesh, value=1))
        >>> coeff = term._getGeomCoeff(mesh)
        >>> M = term._getCoefficientMatrix(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.getNumpyArray(), 
        ...                        (( 1., -1.), 
        ...                         (-1.,  1.))) or procID != 0
        True
        >>> L,b = term._buildMatrix(var=CellVariable(mesh=mesh), SparseMatrix=SparseMatrix)
        >>> print numerix.allclose(L.getNumpyArray(), 
        ...                        ((-1.,  1.), 
        ...                         ( 1., -1.))) or procID != 0
        True
        >>> print numerix.allclose(b, (0., 0.)) or procID != 0
        True

        >>> term = DiffusionTerm(coeff=CellVariable(mesh=mesh, value=1))
        >>> coeff = term._getGeomCoeff(mesh)
        >>> M = term._getCoefficientMatrix(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.getNumpyArray(), 
        ...                        (( 1., -1.), 
        ...                         (-1.,  1.))) or procID != 0
        True
        >>> L,b = term._buildMatrix(var=CellVariable(mesh=mesh), SparseMatrix=SparseMatrix)
        >>> print numerix.allclose(L.getNumpyArray(), 
        ...                        ((-1.,  1.), 
        ...                         ( 1., -1.))) or procID != 0
        True
        >>> print numerix.allclose(b, (0., 0.)) or procID != 0
        True

        >>> from fipy.variables.variable import Variable
        >>> term = DiffusionTerm(coeff = Variable(value = 1))
        >>> coeff = term._getGeomCoeff(mesh)
        >>> M = term._getCoefficientMatrix(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.getNumpyArray(), 
        ...                        (( 1., -1.), 
        ...                         (-1.,  1.))) or procID != 0
        True
        >>> L,b = term._buildMatrix(var=CellVariable(mesh=mesh), SparseMatrix=SparseMatrix)
        >>> print numerix.allclose(L.getNumpyArray(), 
        ...                        ((-1.,  1.), 
        ...                         ( 1., -1.))) or procID != 0
        True
        >>> print numerix.allclose(b, (0., 0.)) or procID != 0
        True
                   
        >>> term = DiffusionTerm(coeff = ((1,2),))

        >>> term = DiffusionTerm(coeff = FaceVariable(mesh = mesh, value = (1,), rank=1))
        >>> term = DiffusionTerm(coeff = CellVariable(mesh=mesh, value=(1,), rank=1))

        Test, 2nd order, 1 dimension, fixed flux 3, fixed value of 4

        >>> from fipy.boundaryConditions.fixedFlux import FixedFlux
        >>> from fipy.boundaryConditions.fixedValue import FixedValue
        >>> bcLeft = FixedFlux(mesh.getFacesLeft(), 3.)
        >>> bcRight = FixedValue(mesh.getFacesRight(), 4.)
        >>> term = DiffusionTerm(coeff = (1.,))
        >>> coeff = term._getGeomCoeff(mesh)
        >>> M = term._getCoefficientMatrix(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.getNumpyArray(), 
        ...                        (( 1., -1.), 
        ...                         (-1.,  1.))) or procID != 0
        True
        >>> L,b = term._buildMatrix(var=CellVariable(mesh=mesh), 
        ...                         SparseMatrix=SparseMatrix ,
        ...                         boundaryConditions=(bcLeft, bcRight))
        >>> print numerix.allclose(L.getNumpyArray(), 
        ...                        ((-1.,  1.), 
        ...                         ( 1., -3.))) or procID != 0
        True
        >>> print numerix.allclose(b, (-3., -8.)) or procID != 0
        True
           
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
        >>> M = term._getCoefficientMatrix(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.getNumpyArray(), 
        ...                        (( 1., -1.), 
        ...                         (-1.,  1.))) or procID != 0
        True
        >>> L,b = term._buildMatrix(var = CellVariable(mesh = mesh), SparseMatrix=SparseMatrix,
        ...                         boundaryConditions = (bcLeft1, bcLeft2, 
        ...                                               bcRight1, bcRight2))
        >>> print numerix.allclose(L.getNumpyArray(), 
        ...                        (( 4., -6.), 
        ...                         (-4., 10.))) or procID != 0
        True
        >>> print numerix.allclose(b, (1., 21.)) or procID != 0
        True
        
        Test, 4th order, 1 dimension, x = 0; fixed flux 3, fixed curvature 2,
        x = 2, fixed value 4, fixed 3rd order -1

        >>> bcLeft1 = FixedFlux(mesh.getFacesLeft(), 3.)
        >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.getFacesLeft(), 2., 2)
        >>> bcRight1 = FixedValue(mesh.getFacesRight(), 4.)
        >>> bcRight2 =  NthOrderBoundaryCondition(mesh.getFacesRight(), -1., 3)
        >>> term = DiffusionTerm(coeff = (-1., 1.))
        >>> coeff = term._getGeomCoeff(mesh)
        >>> M = term._getCoefficientMatrix(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.getNumpyArray(), 
        ...                        ((-1.,  1.), 
        ...                         ( 1., -1.))) or procID != 0
        True
        >>> L,b = term._buildMatrix(var = CellVariable(mesh = mesh), 
        ...                         SparseMatrix=SparseMatrix,
        ...                         boundaryConditions = (bcLeft1, bcLeft2, 
        ...                                               bcRight1, bcRight2))
        >>> print numerix.allclose(L.getNumpyArray(), 
        ...                        ((-4.,  6.), 
        ...                         ( 2., -4.))) or procID != 0
        True
        >>> print numerix.allclose(b, (3., -4.)) or procID != 0
        True


        Test when dx = 0.5.

        >>> mesh = Grid1D(dx = .5, nx = 2)
        >>> bcLeft1 = FixedValue(mesh.getFacesLeft(), 0.)
        >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.getFacesLeft(), 1., 2)
        >>> bcRight1 = FixedFlux(mesh.getFacesRight(), 1.)
        >>> bcRight2 =  NthOrderBoundaryCondition(mesh.getFacesRight(), 0., 3)
        >>> term = DiffusionTerm(coeff = (1., 1.))
        >>> coeff = term._getGeomCoeff(mesh)
        >>> M = term._getCoefficientMatrix(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.getNumpyArray(), 
        ...                        (( 2., -2.), 
        ...                         (-2.,  2.))) or procID != 0
        True
        >>> L,b = term._buildMatrix(var = CellVariable(mesh = mesh), 
        ...                         SparseMatrix=SparseMatrix,
        ...                         boundaryConditions = (bcLeft1, bcLeft2, 
        ...                                               bcRight1, bcRight2))
        >>> print numerix.allclose(L.getNumpyArray(), 
        ...                        (( 80., -32.),
        ...                         (-32.,  16.))) or procID != 0
        True
        >>> print numerix.allclose(b, (-8., 4.)) or procID != 0
        True

        The following tests are to check that DiffusionTerm can take any of the four
        main Variable types.

        >>> from fipy.meshes.tri2D import Tri2D
        >>> mesh = Tri2D(nx = 1, ny = 1)
        >>> term = DiffusionTerm(CellVariable(value = 1, mesh = mesh))
        >>> print term._getGeomCoeff(mesh)[0]
        [ 6.   6.   6.   6.   1.5  1.5  1.5  1.5]
        >>> term = DiffusionTerm(FaceVariable(value = 1, mesh = mesh))
        >>> print term._getGeomCoeff(mesh)[0]
        [ 6.   6.   6.   6.   1.5  1.5  1.5  1.5]
        >>> term = DiffusionTerm(CellVariable(value=(0.5, 1), mesh=mesh, rank=1))
        >>> term = DiffusionTerm(CellVariable(value=((0.5,), (1,)), mesh=mesh, rank=1))
        >>> print term._getGeomCoeff(mesh)[0]
        [ 6.     6.     3.     3.     1.125  1.125  1.125  1.125]
        >>> term = DiffusionTerm(FaceVariable(value=(0.5, 1), mesh=mesh, rank=1))
        >>> term = DiffusionTerm(FaceVariable(value=((0.5,), (1,)), mesh=mesh, rank=1))
        >>> print term._getGeomCoeff(mesh)[0]
        [ 6.     6.     3.     3.     1.125  1.125  1.125  1.125]
        >>> mesh = Tri2D(nx = 1, ny = 1, dy = 0.1)
        >>> term = DiffusionTerm(FaceVariable(value=(0.5, 1), mesh=mesh, rank=1))
        >>> term = DiffusionTerm(FaceVariable(value=((0.5,), (1,)), mesh=mesh, rank=1))
        >>> val = (60., 60., 0.3, 0.3, 0.22277228, 0.22277228, 0.22277228, 0.22277228)
        >>> print numerix.allclose(term._getGeomCoeff(mesh)[0], val)
        1
        >>> term = DiffusionTerm(((0.5, 1),))
        >>> term = DiffusionTerm((((0.5,), (1,)),))
        >>> print numerix.allclose(term._getGeomCoeff(mesh)[0], val)
        Traceback (most recent call last):
            ...
        IndexError: diffusion coefficent tensor is not an appropriate shape for this mesh

        Anisotropy test

        >>> from fipy.meshes.tri2D import Tri2D
        >>> mesh = Tri2D(nx = 1, ny = 1)
        >>> term = DiffusionTerm((((1, 2), (3, 4)),))
        >>> print term._getGeomCoeff(mesh)
        [[ 24.          24.           6.           6.           0.           7.5
            7.5          0.        ]
         [ -3.          -3.           2.           2.          -1.41421356
            0.70710678   0.70710678  -1.41421356]]

        """
        pass

class DiffusionTermNoCorrection(DiffusionTerm):
    def _getNormals(self, mesh):
        return mesh._getFaceNormals()

    def _treatMeshAsOrthogonal(self, mesh):
        return True
        
def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
