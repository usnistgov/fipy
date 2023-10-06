from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.terms.diffusionTermNoCorrection import DiffusionTermNoCorrection
from fipy.terms.diffusionTermCorrection import DiffusionTermCorrection
from fipy.tools import numerix

class DiffusionTerm(DiffusionTermNoCorrection):
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

    def _test(self):
        r"""
        Test, 2nd order, 1 dimension, fixed flux of zero both ends.

        >>> from fipy.meshes import Grid1D
        >>> from fipy.solvers import _MeshMatrix
        >>> from fipy.variables.cellVariable import CellVariable
        >>> mesh = Grid1D(dx = 1., nx = 2)
        >>> term = DiffusionTerm(coeff = (1,))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(_MeshMatrix, CellVariable(mesh=mesh), coeff[0])
        >>> A = M.numpyArray
        >>> print(numerix.allclose(A,
        ...                        (( 1., -1.),
        ...                         (-1.,  1.)))) # doctest: +PROCESSOR_0
        True
        >>> from fipy.variables.cellVariable import CellVariable
        >>> v, L, b = term._buildMatrix(var=CellVariable(mesh=mesh), SparseMatrix=_MeshMatrix)
        >>> A = L.numpyArray
        >>> print(numerix.allclose(A,
        ...                        ((-1.,  1.),
        ...                         ( 1., -1.)))) # doctest: +PROCESSOR_0
        True
        >>> print(numerix.allclose(b, (0., 0.))) # doctest: +PROCESSOR_0
        True

        The coefficient must be a `FaceVariable`, a `CellVariable` (which will
        be interpolated to a `FaceVariable`), or a scalar value

        >>> from fipy.variables.faceVariable import FaceVariable
        >>> term = DiffusionTerm(coeff=FaceVariable(mesh=mesh, value=1))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(_MeshMatrix, CellVariable(mesh=mesh), coeff[0])
        >>> A = M.numpyArray
        >>> print(numerix.allclose(A, (( 1., -1.),
        ...                            (-1.,  1.)))) # doctest: +PROCESSOR_0
        True
        >>> v, L, b = term._buildMatrix(var=CellVariable(mesh=mesh), SparseMatrix=_MeshMatrix)
        >>> A = L.numpyArray
        >>> print(numerix.allclose(A,
        ...                        ((-1.,  1.),
        ...                         ( 1., -1.)))) # doctest: +PROCESSOR_0
        True
        >>> print(numerix.allclose(b, (0., 0.))) # doctest: +PROCESSOR_0
        True

        >>> term = DiffusionTerm(coeff=CellVariable(mesh=mesh, value=1))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(_MeshMatrix, CellVariable(mesh=mesh), coeff[0])
        >>> A = M.numpyArray
        >>> print(numerix.allclose(A,
        ...                        (( 1., -1.),
        ...                         (-1.,  1.)))) # doctest: +PROCESSOR_0
        True
        >>> v, L, b = term._buildMatrix(var=CellVariable(mesh=mesh), SparseMatrix=_MeshMatrix)
        >>> A = L.numpyArray
        >>> print(numerix.allclose(A,
        ...                        ((-1.,  1.),
        ...                         ( 1., -1.)))) # doctest: +PROCESSOR_0
        True
        >>> print(numerix.allclose(b, (0., 0.))) # doctest: +PROCESSOR_0
        True

        >>> from fipy.variables.variable import Variable
        >>> term = DiffusionTerm(coeff = Variable(value = 1))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(_MeshMatrix, CellVariable(mesh=mesh), coeff[0])
        >>> A = M.numpyArray
        >>> print(numerix.allclose(A,
        ...                        (( 1., -1.),
        ...                         (-1.,  1.)))) # doctest: +PROCESSOR_0
        True
        >>> v, L, b = term._buildMatrix(var=CellVariable(mesh=mesh), SparseMatrix=_MeshMatrix)
        >>> A = L.numpyArray
        >>> print(numerix.allclose(A,
        ...                        ((-1.,  1.),
        ...                         ( 1., -1.)))) # doctest: +PROCESSOR_0
        True
        >>> print(numerix.allclose(b, (0., 0.))) # doctest: +PROCESSOR_0
        True

        >>> term = DiffusionTerm(coeff = ((1, 2),))

        >>> term = DiffusionTerm(coeff = FaceVariable(mesh = mesh, value = (1,), rank=1))
        >>> term = DiffusionTerm(coeff = CellVariable(mesh=mesh, value=(1,), rank=1))

        Test, 2nd order, 1 dimension, fixed flux 3, fixed value of 4

        >>> var=CellVariable(mesh=mesh)
        >>> var.faceGrad.constrain([-3.], mesh.facesLeft)
        >>> var.constrain(4., mesh.facesRight)
        >>> term = DiffusionTerm(coeff = (1.,))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(_MeshMatrix, CellVariable(mesh=mesh), coeff[0])
        >>> A = M.numpyArray
        >>> print(numerix.allclose(A,
        ...                        (( 1., -1.),
        ...                         (-1.,  1.)))) # doctest: +PROCESSOR_0
        True
        >>> v, L, b = term._buildMatrix(var=var,
        ...                         SparseMatrix=_MeshMatrix)
        >>> A = L.numpyArray
        >>> print(numerix.allclose(A,
        ...                        ((-1.,  1.),
        ...                         ( 1., -3.)))) # doctest: +PROCESSOR_0
        True
        >>> print(numerix.allclose(b, (-3., -8.))) # doctest: +PROCESSOR_0
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
        >>> M = term._getCoefficientMatrixForTests(_MeshMatrix, CellVariable(mesh=mesh), coeff[0])
        >>> A = M.numpyArray
        >>> print(numerix.allclose(A,
        ...                        (( 1., -1.),
        ...                         (-1.,  1.)))) # doctest: +PROCESSOR_0
        True

        >>> v, L, b = term._buildMatrix(var=var, SparseMatrix=_MeshMatrix,
        ...                         boundaryConditions=(bcLeft2, bcRight2))
        >>> A = L.numpyArray
        >>> print(numerix.allclose(A,
        ...                        (( 4., -6.),
        ...                         (-4., 10.)))) # doctest: +PROCESSOR_0
        True
        >>> print(numerix.allclose(b, (1., 21.))) # doctest: +PROCESSOR_0
        True

        Test, 4th order, 1 dimension, x = 0; fixed flux 3, fixed curvature 2,
        x = 2, fixed value 4, fixed 3rd order `-1`

        >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.facesLeft, 2., 2)
        >>> bcRight2 =  NthOrderBoundaryCondition(mesh.facesRight, -1., 3)

        >>> var = CellVariable(mesh=mesh)
        >>> var.faceGrad.constrain([-3.], mesh.facesLeft)
        >>> var.constrain(4., mesh.facesRight)

        >>> term = DiffusionTerm(coeff = (-1., 1.))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(_MeshMatrix, CellVariable(mesh=mesh), coeff[0])
        >>> A = M.numpyArray
        >>> print(numerix.allclose(A,
        ...                        ((-1.,  1.),
        ...                         ( 1., -1.)))) # doctest: +PROCESSOR_0
        True

        >>> v, L, b = term._buildMatrix(var=var,
        ...                         SparseMatrix=_MeshMatrix,
        ...                         boundaryConditions = (bcLeft2, bcRight2))
        >>> A = L.numpyArray
        >>> print(numerix.allclose(A,
        ...                        ((-4.,  6.),
        ...                         ( 2., -4.)))) # doctest: +PROCESSOR_0
        True
        >>> print(numerix.allclose(b, (3., -4.))) # doctest: +PROCESSOR_0
        True


        Test when `dx = 0.5`.

        >>> mesh = Grid1D(dx = .5, nx = 2)

        >>> bcLeft2 =  NthOrderBoundaryCondition(mesh.facesLeft, 1., 2)
        >>> bcRight2 =  NthOrderBoundaryCondition(mesh.facesRight, 0., 3)

        >>> var = CellVariable(mesh=mesh)
        >>> var.faceGrad.constrain([1.], mesh.facesRight)
        >>> var.constrain(0., mesh.facesLeft)

        >>> term = DiffusionTerm(coeff = (1., 1.))
        >>> coeff = term._getGeomCoeff(CellVariable(mesh=mesh))
        >>> M = term._getCoefficientMatrixForTests(_MeshMatrix, CellVariable(mesh=mesh), coeff[0])
        >>> A = M.numpyArray
        >>> print(numerix.allclose(A,
        ...                        (( 2., -2.),
        ...                         (-2.,  2.)))) # doctest: +PROCESSOR_0
        True

        >>> v, L, b = term._buildMatrix(var=var,
        ...                         SparseMatrix=_MeshMatrix,
        ...                         boundaryConditions = (bcLeft2, bcRight2))
        >>> A = L.numpyArray
        >>> print(numerix.allclose(A,
        ...                        (( 80., -32.),
        ...                         (-32.,  16.)))) # doctest: +PROCESSOR_0
        True
        >>> print(numerix.allclose(b, (-8., 4.))) # doctest: +PROCESSOR_0
        True

        The following tests are to check that `DiffusionTerm` can take any of the four
        main Variable types.

        >>> from fipy.meshes.tri2D import Tri2D
        >>> mesh = Tri2D(nx = 1, ny = 1)
        >>> term = DiffusionTermCorrection(CellVariable(value = 1, mesh = mesh))
        >>> print(term._getGeomCoeff(CellVariable(mesh=mesh))[0])
        [ 6.   6.   6.   6.   1.5  1.5  1.5  1.5]
        >>> term = DiffusionTerm(FaceVariable(value = 1, mesh = mesh))
        >>> print(term._getGeomCoeff(CellVariable(mesh=mesh))[0])
        [ 6.   6.   6.   6.   1.5  1.5  1.5  1.5]
        >>> term = DiffusionTerm(CellVariable(value=((0.5,), (1,)), mesh=mesh, rank=1))
        >>> print(term._getGeomCoeff(CellVariable(mesh=mesh))[0])
        [ 6.     6.     3.     3.     1.125  1.125  1.125  1.125]
        >>> term = DiffusionTerm(FaceVariable(value=((0.5,), (1,)), mesh=mesh, rank=1))
        >>> print(term._getGeomCoeff(CellVariable(mesh=mesh))[0])
        [ 6.     6.     3.     3.     1.125  1.125  1.125  1.125]
        >>> mesh = Tri2D(nx = 1, ny = 1, dy = 0.1)
        >>> term = DiffusionTermCorrection(FaceVariable(value=((0.5,), (1,)), mesh=mesh, rank=1))
        >>> val = (60., 60., 0.3, 0.3, 0.22277228, 0.22277228, 0.22277228, 0.22277228)
        >>> print(numerix.allclose(term._getGeomCoeff(CellVariable(mesh=mesh))[0], val))
        1
        >>> term = DiffusionTerm((((0.5,), (1,)),))
        >>> print(numerix.allclose(term._getGeomCoeff(CellVariable(mesh=mesh))[0], val))
        Traceback (most recent call last):
            ...
        IndexError: diffusion coefficient tensor is not an appropriate shape for this mesh

        Anisotropy test

        >>> from fipy.meshes.tri2D import Tri2D
        >>> mesh = Tri2D(nx = 1, ny = 1)
        >>> term = DiffusionTerm((((1, 2), (3, 4)),))
        >>> print(numerix.allclose(term._getGeomCoeff(CellVariable(mesh=mesh)),
        ...                        [[ 24.,         24.,          6.,          6.,          0.,          7.5,
        ...                        7.5,         0.        ],
        ...                         [ -3.,         -3.,          2.,          2.,         -1.41421356,
        ...                        0.70710678,   0.70710678,  -1.41421356]]))
        True

        Negate the term.

        >>> -DiffusionTerm(coeff=[1.])
        DiffusionTerm(coeff=[-1.0])

        >>> -DiffusionTerm()
        DiffusionTerm(coeff=[-1.0])

        Testing vector diffusion terms by comparing with coupled solutions.

        >>> from fipy import *
        >>> m = Grid1D(nx=6)
        >>> q = CellVariable(mesh=m, elementshape=(2,))
        >>> coeff = FaceVariable(mesh=m, elementshape=(2, 2), value=[[1, 2], [3, 4]])
        >>> vectorEq = DiffusionTerm(coeff) == 0
        >>> vectorEq.cacheMatrix()
        >>> vectorEq.solve(q, solver=DummySolver())

        >>> v0 = CellVariable(mesh=m)
        >>> v1 = CellVariable(mesh=m)
        >>> coupledEq = ((DiffusionTerm(coeff=1., var=v0) + DiffusionTerm(coeff=2., var=v1)) & (DiffusionTerm(coeff=3., var=v0) + DiffusionTerm(coeff=4., var=v1)))((v0, v1))
        >>> coupledEq.cacheMatrix()
        >>> coupledEq.solve(solver=DummySolver())
        >>> print((coupledEq.matrix.numpyArray == vectorEq.matrix.numpyArray).all())
        True

        Test vector diffusion terms.

        >>> m = Grid2D(nx=2, ny=2)
        >>> q = CellVariable(mesh=m, elementshape=(3,))
        >>> coeff = FaceVariable(mesh=m, elementshape=(3, 3), value=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> vectorEq0 = DiffusionTerm(coeff) == 0
        >>> vectorEq0.cacheMatrix()
        >>> vectorEq0.solve(q, solver=DummySolver())

        Test vector diffusion terms when the coefficient is a tuple or list.

        >>> m = Grid2D(nx=2, ny=2)
        >>> q = CellVariable(mesh=m, elementshape=(3,))
        >>> vectorEq1 = DiffusionTerm(([[1, 2, 3], [4, 5, 6], [7, 8, 9]],)) == 0
        >>> vectorEq1.cacheMatrix()
        >>> vectorEq1.solve(q, solver=DummySolver())

        The following test tests for negating a tuple/list coefficient for the diffusion term.

        >>> m = Grid2D(nx=2, ny=2)
        >>> q = CellVariable(mesh=m, elementshape=(3,))
        >>> vectorEq2 = TransientTerm(1e-20) == DiffusionTerm(([[1, 2, 3], [4, 5, 6], [7, 8, 9]],))
        >>> vectorEq2.cacheMatrix()
        >>> vectorEq2.solve(q, solver=DummySolver(), dt=1.)

        Test the previous three vector examples against coupled.

        >>> v0 = CellVariable(mesh=m)
        >>> v1 = CellVariable(mesh=m)
        >>> v2 = CellVariable(mesh=m)
        >>> coupledEq = ((DiffusionTerm(coeff=1., var=v0) + DiffusionTerm(coeff=2., var=v1) + DiffusionTerm(coeff=3., var=v2)) & \
        ...              (DiffusionTerm(coeff=4., var=v0) + DiffusionTerm(coeff=5., var=v1) + DiffusionTerm(coeff=6., var=v2)) & \
        ...              (DiffusionTerm(coeff=7., var=v0) + DiffusionTerm(coeff=8., var=v1) + DiffusionTerm(coeff=9, var=v2)))((v0, v1, v2))
        >>> coupledEq.cacheMatrix()
        >>> coupledEq.solve(solver=DummySolver())
        >>> coupledMatrix = coupledEq.matrix.numpyArray
        >>> print((coupledMatrix == vectorEq0.matrix.numpyArray).all())
        True
        >>> print((coupledMatrix == vectorEq1.matrix.numpyArray).all())
        True
        >>> print((coupledMatrix == -vectorEq2.matrix.numpyArray).all())
        True

        """
        pass



def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


