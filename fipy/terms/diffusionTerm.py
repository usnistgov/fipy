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

from fipy.terms.baseDiffusionTerm import _BaseDiffusionTerm
from fipy.tools import numerix

class DiffusionTerm(_BaseDiffusionTerm):
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

    def _getNormals(self, mesh):
        return mesh._faceCellToCellNormals

    def _treatMeshAsOrthogonal(self, mesh):        
        return mesh._isOrthogonal()

    def _test(self):
        r"""
        Test, 2nd order, 1 dimension, fixed flux of zero both ends.

        >>> from fipy.meshes import Grid1D
        >>> from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix as SparseMatrix
        >>> from fipy.tools import parallel
        >>> from fipy.variables.cellVariable import CellVariable
        >>> procID = parallel.procID
        >>> mesh = Grid1D(dx = 1., nx = 2)
        >>> term = DiffusionTerm(coeff = (1,))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.numpyArray, 
        ...                        (( 1., -1.), 
        ...                         (-1.,  1.))) or procID != 0
        True
        >>> from fipy.variables.cellVariable import CellVariable
        >>> v,L,b = term._buildMatrix(var=CellVariable(mesh=mesh), SparseMatrix=SparseMatrix)
        >>> print numerix.allclose(L.numpyArray, 
        ...                        ((-1.,  1.), 
        ...                         ( 1., -1.))) or procID != 0
        True
        >>> print numerix.allclose(b, (0., 0.)) or procID != 0
        True

        The coefficient must be a `FaceVariable`, a `CellVariable` (which will
        be interpolated to a `FaceVariable`), or a scalar value 
        
        >>> from fipy.variables.faceVariable import FaceVariable
        >>> term = DiffusionTerm(coeff=FaceVariable(mesh=mesh, value=1))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.numpyArray, 
        ...                        (( 1., -1.), 
        ...                         (-1.,  1.))) or procID != 0
        True
        >>> v,L,b = term._buildMatrix(var=CellVariable(mesh=mesh), SparseMatrix=SparseMatrix)
        >>> print numerix.allclose(L.numpyArray, 
        ...                        ((-1.,  1.), 
        ...                         ( 1., -1.))) or procID != 0
        True
        >>> print numerix.allclose(b, (0., 0.)) or procID != 0
        True

        >>> term = DiffusionTerm(coeff=CellVariable(mesh=mesh, value=1))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.numpyArray, 
        ...                        (( 1., -1.), 
        ...                         (-1.,  1.))) or procID != 0
        True
        >>> v,L,b = term._buildMatrix(var=CellVariable(mesh=mesh), SparseMatrix=SparseMatrix)
        >>> print numerix.allclose(L.numpyArray, 
        ...                        ((-1.,  1.), 
        ...                         ( 1., -1.))) or procID != 0
        True
        >>> print numerix.allclose(b, (0., 0.)) or procID != 0
        True

        >>> from fipy.variables.variable import Variable
        >>> term = DiffusionTerm(coeff = Variable(value = 1))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.numpyArray, 
        ...                        (( 1., -1.), 
        ...                         (-1.,  1.))) or procID != 0
        True
        >>> v,L,b = term._buildMatrix(var=CellVariable(mesh=mesh), SparseMatrix=SparseMatrix)
        >>> print numerix.allclose(L.numpyArray, 
        ...                        ((-1.,  1.), 
        ...                         ( 1., -1.))) or procID != 0
        True
        >>> print numerix.allclose(b, (0., 0.)) or procID != 0
        True
                   
        >>> term = DiffusionTerm(coeff = ((1,2),))

        >>> term = DiffusionTerm(coeff = FaceVariable(mesh = mesh, value = (1,), rank=1))
        >>> term = DiffusionTerm(coeff = CellVariable(mesh=mesh, value=(1,), rank=1))

        Test, 2nd order, 1 dimension, fixed flux 3, fixed value of 4

        >>> var=CellVariable(mesh=mesh)
        >>> var.faceGrad.constrain([-3.], mesh.facesLeft)
        >>> var.constrain(4., mesh.facesRight)
        >>> term = DiffusionTerm(coeff = (1.,))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.numpyArray, 
        ...                        (( 1., -1.), 
        ...                         (-1.,  1.))) or procID != 0
        True
        >>> v,L,b = term._buildMatrix(var=var, 
        ...                         SparseMatrix=SparseMatrix)
        >>> print numerix.allclose(L.numpyArray, 
        ...                        ((-1.,  1.), 
        ...                         ( 1., -3.))) or procID != 0
        True
        >>> print numerix.allclose(b, (-3., -8.)) or procID != 0
        True
           
        Test, 4th order, 1 dimension, x = 0; fixed flux 3, fixed curvatures 0,
        x = 2, fixed value 1, fixed curvature 0

        >>> from fipy.boundaryConditions.nthOrderBoundaryCondition \
        ...     import NthOrderBoundaryCondition
        >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.facesLeft, 0., 2)
        >>> var = CellVariable(mesh=mesh)
        >>> var.faceGrad.constrain([-3.], mesh.facesLeft)
        >>> var.constrain(4., mesh.facesRight)
        >>> bcRight2 =  NthOrderBoundaryCondition(mesh.facesRight, 0., 2)
        >>> term = DiffusionTerm(coeff = (1., 1.))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.numpyArray, 
        ...                        (( 1., -1.), 
        ...                         (-1.,  1.))) or procID != 0
        True

        >>> v,L,b = term._buildMatrix(var=var, SparseMatrix=SparseMatrix,
        ...                         boundaryConditions=(bcLeft2, bcRight2))
        >>> print numerix.allclose(L.numpyArray, 
        ...                        (( 4., -6.), 
        ...                         (-4., 10.))) or procID != 0
        True
        >>> print numerix.allclose(b, (1., 21.)) or procID != 0
        True
        
        Test, 4th order, 1 dimension, x = 0; fixed flux 3, fixed curvature 2,
        x = 2, fixed value 4, fixed 3rd order -1

        >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.facesLeft, 2., 2)
        >>> bcRight2 =  NthOrderBoundaryCondition(mesh.facesRight, -1., 3)

        >>> var = CellVariable(mesh=mesh)
        >>> var.faceGrad.constrain([-3.], mesh.facesLeft)
        >>> var.constrain(4., mesh.facesRight)

        >>> term = DiffusionTerm(coeff = (-1., 1.))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.numpyArray, 
        ...                        ((-1.,  1.), 
        ...                         ( 1., -1.))) or procID != 0
        True

        >>> v,L,b = term._buildMatrix(var=var,
        ...                         SparseMatrix=SparseMatrix,
        ...                         boundaryConditions = (bcLeft2, bcRight2))
        
        >>> print numerix.allclose(L.numpyArray, 
        ...                        ((-4.,  6.), 
        ...                         ( 2., -4.))) or procID != 0
        True
        >>> print numerix.allclose(b, (3., -4.)) or procID != 0
        True


        Test when dx = 0.5.

        >>> mesh = Grid1D(dx = .5, nx = 2)

        >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.facesLeft, 1., 2)
        >>> bcRight2 =  NthOrderBoundaryCondition(mesh.facesRight, 0., 3)

        >>> var = CellVariable(mesh=mesh)
        >>> var.faceGrad.constrain([1.], mesh.facesRight)
        >>> var.constrain(0., mesh.facesLeft)
        
        >>> term = DiffusionTerm(coeff = (1., 1.))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(SparseMatrix, mesh, coeff[0])
        >>> print numerix.allclose(M.numpyArray, 
        ...                        (( 2., -2.), 
        ...                         (-2.,  2.))) or procID != 0
        True

        >>> v,L,b = term._buildMatrix(var=var, 
        ...                         SparseMatrix=SparseMatrix,
        ...                         boundaryConditions = (bcLeft2, bcRight2))

        >>> print numerix.allclose(L.numpyArray, 
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
        >>> print term._getGeomCoeff(CellVariable(mesh=mesh))[0]
        [ 6.   6.   6.   6.   1.5  1.5  1.5  1.5]
        >>> term = DiffusionTerm(FaceVariable(value = 1, mesh = mesh))
        >>> print term._getGeomCoeff(CellVariable(mesh=mesh))[0]
        [ 6.   6.   6.   6.   1.5  1.5  1.5  1.5]
        >>> term = DiffusionTerm(CellVariable(value=(0.5, 1), mesh=mesh, rank=1))
        >>> term = DiffusionTerm(CellVariable(value=((0.5,), (1,)), mesh=mesh, rank=1))
        >>> print term._getGeomCoeff(CellVariable(mesh=mesh))[0]
        [ 6.     6.     3.     3.     1.125  1.125  1.125  1.125]
        >>> term = DiffusionTerm(FaceVariable(value=(0.5, 1), mesh=mesh, rank=1))
        >>> term = DiffusionTerm(FaceVariable(value=((0.5,), (1,)), mesh=mesh, rank=1))
        >>> print term._getGeomCoeff(CellVariable(mesh=mesh))[0]
        [ 6.     6.     3.     3.     1.125  1.125  1.125  1.125]
        >>> mesh = Tri2D(nx = 1, ny = 1, dy = 0.1)
        >>> term = DiffusionTerm(FaceVariable(value=(0.5, 1), mesh=mesh, rank=1))
        >>> term = DiffusionTerm(FaceVariable(value=((0.5,), (1,)), mesh=mesh, rank=1))
        >>> val = (60., 60., 0.3, 0.3, 0.22277228, 0.22277228, 0.22277228, 0.22277228)
        >>> print numerix.allclose(term._getGeomCoeff(CellVariable(mesh=mesh))[0], val)
        1
        >>> term = DiffusionTerm(((0.5, 1),))
        >>> term = DiffusionTerm((((0.5,), (1,)),))
        >>> print numerix.allclose(term._getGeomCoeff(CellVariable(mesh=mesh))[0], val)
        Traceback (most recent call last):
            ...
        IndexError: diffusion coefficent tensor is not an appropriate shape for this mesh

        Anisotropy test

        >>> from fipy.meshes.tri2D import Tri2D
        >>> mesh = Tri2D(nx = 1, ny = 1)
        >>> term = DiffusionTerm((((1, 2), (3, 4)),))
        >>> print term._getGeomCoeff(CellVariable(mesh=mesh))
        [[ 24.          24.           6.           6.           0.           7.5
            7.5          0.        ]
         [ -3.          -3.           2.           2.          -1.41421356
            0.70710678   0.70710678  -1.41421356]]

        Negate the term.

        >>> -DiffusionTerm(coeff=[1.])
        DiffusionTerm(coeff=[-1.0])

        >>> -DiffusionTerm()
        DiffusionTerm(coeff=[-1.0])

        Testing vector diffusion terms by comparing with coupled solutions.

        >>> m = Grid1D(nx=6)
        >>> v0 = CellVariable(mesh=m)
        >>> v1 = CellVariable(mesh=m)
        >>> eq = (DiffusionTerm(coeff=1., var=v0) + DiffusionTerm(coeff=2., var=v1)) & (DiffusionTerm(coeff=3., var=v0) + DiffusionTerm(coeff=4., var=v1))
        >>> eq.cacheMatrix()
        >>> eq.solve()
        >>> print eq.matrix

        """
        pass            

from fipy.terms.diffusionTermNoCorrection import DiffusionTermNoCorrection 
        
def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
