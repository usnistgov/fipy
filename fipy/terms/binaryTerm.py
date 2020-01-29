from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

import os

from fipy.terms.abstractBinaryTerm import _AbstractBinaryTerm

class _BinaryTerm(_AbstractBinaryTerm):

    @property
    def _buildExplcitIfOther(self):
        return True

    def _buildAndAddMatrices(self, var, SparseMatrix,  boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None, buildExplicitIfOther=True):
        """Build matrices of constituent Terms and collect them

        Only called at top-level by `_prepareLinearSystem()`

        """

        matrix = SparseMatrix(mesh=var.mesh)
        RHSvector = 0

        for term in (self.term, self.other):

            tmpVar, tmpMatrix, tmpRHSvector = term._buildAndAddMatrices(var,
                                                                        SparseMatrix,
                                                                        boundaryConditions=boundaryConditions,
                                                                        dt=dt,
                                                                        transientGeomCoeff=transientGeomCoeff,
                                                                        diffusionGeomCoeff=diffusionGeomCoeff,
                                                                        buildExplicitIfOther=buildExplicitIfOther)

            matrix += tmpMatrix
            RHSvector += tmpRHSvector

            term._buildCache(tmpMatrix, tmpRHSvector)

        return (var, matrix, RHSvector)

    def _getDefaultSolver(self, var, solver, *args, **kwargs):
        for term in (self.term, self.other):
            defaultSolver = term._getDefaultSolver(var, solver, *args, **kwargs)
            if defaultSolver is not None:
                solver = defaultSolver

        return solver

    def __repr__(self):
        return '(' + repr(self.term) + ' + ' + repr(self.other) + ')'

    def __mul__(self, other):
        return other * self.term + other * self.other

    @property
    def _uncoupledTerms(self):
        return [self]

    def _getTransientGeomCoeff(self, var):
        return self._addNone(self.term._getTransientGeomCoeff(var), self.other._getTransientGeomCoeff(var))

    def _getDiffusionGeomCoeff(self, var):
        return self._addNone(self.term._getDiffusionGeomCoeff(var), self.other._getDiffusionGeomCoeff(var))

    __rmul__ = __mul__

    def _test(self):
        """
        >>> from fipy import *
        >>> m = Grid1D(nx=3)
        >>> v0 = CellVariable(mesh=m, value=0.)
        >>> v1 = CellVariable(mesh=m, value=1.)
        >>> eq = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0) - DiffusionTerm(coeff=2., var=v1)
        >>> var, matrix, RHSvector = eq._buildAndAddMatrices(var=v0, SparseMatrix=DefaultSolver()._matrixClass, dt=1.)
        >>> print(var)
        [ 0.  0.  0.]
        >>> print(CellVariable(mesh=m, value=RHSvector))
        [ 0.  0.  0.]
        >>> print(numerix.allequal(matrix.numpyArray, [[ 2, -1,  0],
        ...                                            [-1,  3, -1],
        ...                                            [ 0, -1,  2]]))
        True
        >>> var, matrix, RHSvector = eq._buildAndAddMatrices(var=v1, SparseMatrix=DefaultSolver()._matrixClass, dt=1.)
        >>> print(var)
        [ 1.  1.  1.]
        >>> print(CellVariable(mesh=m, value=RHSvector))
        [ 0.  0.  0.]
        >>> print(numerix.allequal(matrix.numpyArray, [[ 2, -2,  0],
        ...                                            [-2,  4, -2],
        ...                                            [ 0, -2,  2]]))
        True
        >>> print(CellVariable(mesh=m, value=eq.justResidualVector(dt=1.)))
        [ 0.  0.  0.]

        >>> m = Grid1D(nx=6)
        >>> v0 = CellVariable(mesh=m, value=1.)
        >>> v1 = CellVariable(mesh=m, value=0.)
        >>> eq = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0) - DiffusionTerm(coeff=2., var=v1)
        >>> var, matrix, RHSvector = eq._buildAndAddMatrices(var=v0, SparseMatrix=DefaultSolver()._matrixClass, dt=1.)
        >>> print(var)
        [ 1.  1.  1.  1.  1.  1.]
        >>> print(CellVariable(mesh=m, value=RHSvector))
        [ 1.  1.  1.  1.  1.  1.]
        >>> print(numerix.allequal(matrix.numpyArray, [[ 2, -1, 0, 0, 0, 0.],
        ...                                            [-1, 3, -1, 0, 0, 0.],
        ...                                            [ 0, -1, 3, -1, 0, 0.],
        ...                                            [ 0, 0, -1, 3, -1, 0.],
        ...                                            [ 0, 0, 0, -1, 3, -1.],
        ...                                            [ 0, 0, 0, 0, -1, 2.]]))
        True
        >>> var, matrix, RHSvector = eq._buildAndAddMatrices(var=v1, SparseMatrix=DefaultSolver()._matrixClass, dt=1.)
        >>> print(var)
        [ 0.  0.  0.  0.  0.  0.]
        >>> print(CellVariable(mesh=m, value=RHSvector))
        [ 0.  0.  0.  0.  0.  0.]
        >>> print(numerix.allequal(matrix.numpyArray, [[ 2, -2, 0, 0, 0, 0.],
        ...                                            [-2, 4, -2, 0, 0, 0.],
        ...                                            [ 0, -2, 4, -2, 0, 0.],
        ...                                            [ 0, 0, -2, 4, -2, 0.],
        ...                                            [ 0, 0, 0, -2, 4, -2.],
        ...                                            [ 0, 0, 0, 0, -2, 2.]]))
        True
        >>> value = eq.justResidualVector(dt=1.)
        >>> print(CellVariable(mesh=m, value=value))
        [ 0.  0.  0.  0.  0.  0.]

        >>> m = Grid1D(nx=3)
        >>> v0 = CellVariable(mesh=m, value=(0., 1., 2.))
        >>> v1 = CellVariable(mesh=m, value=(3., 4., 5.))
        >>> diffTerm = DiffusionTerm(coeff=1., var=v0)
        >>> eq00 = TransientTerm(var=v0) - diffTerm
        >>> eq0 = eq00 - DiffusionTerm(coeff=2., var=v1)
        >>> eq0.cacheMatrix()
        >>> diffTerm.cacheMatrix()
        >>> print(CellVariable(mesh=m, value=eq0.justResidualVector(dt=1.)))
        [-3.  0.  3.]
        >>> eq0.solve(var=v0, solver=DummySolver(), dt=1.)
        >>> print(numerix.allequal(eq0.matrix.numpyArray, [[ 2, -1,  0],
        ...                                                [-1,  3, -1],
        ...                                                [ 0, -1,  2]]))
        True
        >>> eq0.solve(var=v1, solver=DummySolver(), dt=1.)
        >>> print(numerix.allequal(eq0.matrix.numpyArray, [[ 2, -2,  0],
        ...                                                [-2,  4, -2],
        ...                                                [ 0, -2,  2]]))
        True
        >>> ## This correctly returns None because we lost the handle to the DiffusionTerm when it's negated.
        >>> print(diffTerm.matrix)
        None

        Testing solution for one variable in a multi-variable equation.

        >>> from fipy import *
        >>> L = 1.
        >>> nx = 3
        >>> m = Grid1D(nx=nx, dx=L / nx)
        >>> x = m.cellCenters[0]
        >>> v0 = CellVariable(mesh=m, value=0., name='v0')
        >>> v0.constrain(0., where=m.facesLeft)
        >>> v0.constrain(L, where=m.facesRight)
        >>> v1 = CellVariable(mesh=m, value=-x**2, name='v1')
        >>> v1.constrain(0.,  where=m.facesLeft)
        >>> v1.constrain(-L,  where=m.facesRight)
        >>> (DiffusionTerm(var=v0) + DiffusionTerm(var=v1)).solve(v0)
        >>> print(numerix.allclose(v0, -v1))
        True

        Testing for vector equations.

        >>> m = Grid1D(nx=6)
        >>> v = CellVariable(mesh=m, elementshape=(2,))
        >>> v[0] = 1.
        >>> v[1] = 2.
        >>> (TransientTerm(v) == ImplicitSourceTerm(-v)).solve(v, dt=1.)
        >>> print(v)
        [[ 0.5  0.5  0.5  0.5  0.5  0.5]
         [ 1.   1.   1.   1.   1.   1. ]]
        >>> (TransientTerm() == v).solve(v, dt=1.)
        >>> print(v)
        [[ 1.  1.  1.  1.  1.  1.]
         [ 2.  2.  2.  2.  2.  2.]]
        >>> (TransientTerm(v) == v).solve(v, dt=1.)
        >>> print(v)
        [[ 2.  2.  2.  2.  2.  2.]
         [ 3.  3.  3.  3.  3.  3.]]

        >>> v[0] = 1.
        >>> v[1] = 2.
        >>> (TransientTerm(v) == ImplicitSourceTerm(-v * numerix.identity(2)[..., numerix.newaxis])).solve(v, dt=1.)
        >>> print(v)
        [[ 0.5  0.5  0.5  0.5  0.5  0.5]
         [ 1.   1.   1.   1.   1.   1. ]]
        >>> (TransientTerm() == numerix.array((0.5, 1.))).solve(v, dt=1.)
        >>> print(v)
        [[ 1.  1.  1.  1.  1.  1.]
         [ 2.  2.  2.  2.  2.  2.]]
        >>> (TransientTerm(v * numerix.identity(2)[..., numerix.newaxis]) == v).solve(v, dt=1.)
        >>> print(v)
        [[ 2.  2.  2.  2.  2.  2.]
         [ 3.  3.  3.  3.  3.  3.]]

        >>> v[0] = 1.
        >>> v[1] = 2.
        >>> (TransientTerm(v) == -ImplicitSourceTerm(((1, 0), (0, 2)))).solve(v, dt=1.)
        >>> print(v)
        [[ 0.5  0.5  0.5  0.5  0.5  0.5]
         [ 1.   1.   1.   1.   1.   1. ]]

        >>> v = CellVariable(mesh=m, elementshape=(2,))
        >>> v[0] = 1
        >>> v[1] = 2
        >>> coeff = FaceVariable(mesh=m, elementshape=(1, 2, 2))
        >>> coeff[0, 0, 1] = 1
        >>> coeff[0, 1, 0] = 2
        >>> eqn = TransientTerm() + CentralDifferenceConvectionTerm(coeff=coeff)
        >>> eqn.cacheMatrix()
        >>> eqn.cacheRHSvector()
        >>> eqn.solve(v, dt=1)
        >>> print(numerix.allequal(eqn.matrix.numpyArray,
        ... [[ 1.0,   0,    0,    0,    0,    0,   0.5,  0.5,   0,    0,    0,    0,  ],
        ...  [  0,   1.0,   0,    0,    0,    0,  -0.5,   0,   0.5,   0,    0,    0,  ],
        ...  [  0,    0,   1.0,   0,    0,    0,    0,  -0.5,   0,   0.5,   0,    0,  ],
        ...  [  0,    0,    0,   1.0,   0,    0,    0,    0,  -0.5,   0,   0.5,   0,  ],
        ...  [  0,    0,    0,    0,   1.0,   0,    0,    0,    0,  -0.5,   0,   0.5, ],
        ...  [  0,    0,    0,    0,    0,   1.0,   0,    0,    0,    0,  -0.5, -0.5, ],
        ...  [ 1.0,  1.0,   0,    0,    0,    0,   1.0,   0,    0,    0,    0,    0,  ],
        ...  [-1.0,   0,   1.0,   0,    0,    0,    0,   1.0,   0,    0,    0,    0,  ],
        ...  [  0,  -1.0,   0,   1.0,   0,    0,    0,    0,   1.0,   0,    0,    0,  ],
        ...  [  0,    0,  -1.0,    0,  1.0,   0,    0,    0,    0,   1.0,   0,    0,  ],
        ...  [  0,    0,    0,  -1.0,   0,   1.0,   0,    0,    0,    0,   1.0,   0,  ],
        ...  [  0,    0,    0,    0,  -1.0, -1.0,   0,    0,    0,    0,    0,   1.0, ]]))
        True
        >>> LHS = CellVariable(mesh=m, rank=1, elementshape=(2,), value=numerix.reshape(numerix.asarray(eqn.matrix * v.value.ravel()), (2, -1))).globalValue.ravel()
        >>> RHS = CellVariable(mesh=m, rank=1, elementshape=(2,), value=numerix.reshape(eqn.RHSvector, (2, -1))).globalValue.ravel()
        >>> print(numerix.allclose(LHS, RHS))
        True

        >>> v[0] = 1
        >>> v[1] = 2
        >>> eqn = TransientTerm() + CentralDifferenceConvectionTerm(coeff=(((0, 1), (2, 0)),))
        >>> eqn.cacheMatrix()
        >>> eqn.cacheRHSvector()
        >>> eqn.solve(v, dt=1)
        >>> print(numerix.allequal(eqn.matrix.numpyArray,
        ... [[ 1.0,   0,    0,    0,    0,    0,   0.5,  0.5,   0,    0,    0,    0,  ],
        ...  [  0,   1.0,   0,    0,    0,    0,  -0.5,   0,   0.5,   0,    0,    0,  ],
        ...  [  0,    0,   1.0,   0,    0,    0,    0,  -0.5,   0,   0.5,   0,    0,  ],
        ...  [  0,    0,    0,   1.0,   0,    0,    0,    0,  -0.5,   0,   0.5,   0,  ],
        ...  [  0,    0,    0,    0,   1.0,   0,    0,    0,    0,  -0.5,   0,   0.5, ],
        ...  [  0,    0,    0,    0,    0,   1.0,   0,    0,    0,    0,  -0.5, -0.5, ],
        ...  [ 1.0,  1.0,   0,    0,    0,    0,   1.0,   0,    0,    0,    0,    0,  ],
        ...  [-1.0,   0,   1.0,   0,    0,    0,    0,   1.0,   0,    0,    0,    0,  ],
        ...  [  0,  -1.0,   0,   1.0,   0,    0,    0,    0,   1.0,   0,    0,    0,  ],
        ...  [  0,    0,  -1.0,    0,  1.0,   0,    0,    0,    0,   1.0,   0,    0,  ],
        ...  [  0,    0,    0,  -1.0,   0,   1.0,   0,    0,    0,    0,   1.0,   0,  ],
        ...  [  0,    0,    0,    0,  -1.0, -1.0,   0,    0,    0,    0,    0,   1.0, ]]))
        True
        >>> LHS = CellVariable(mesh=m, rank=1, elementshape=(2,), value=numerix.reshape(numerix.asarray(eqn.matrix * v.value.ravel()), (2, -1))).globalValue.ravel()
        >>> RHS = CellVariable(mesh=m, rank=1, elementshape=(2,), value=numerix.reshape(eqn.RHSvector, (2, -1))).globalValue.ravel()
        >>> print(numerix.allclose(LHS, RHS))
        True


        >>> X = m.faceCenters[0]
        >>> v[0] = 1
        >>> v[1] = 2
        >>> coeff[0, 0, 1] = numerix.sign(X - 3)
        >>> coeff[0, 1, 0] = -2 * numerix.sign(X - 3)
        >>> eqn = TransientTerm() + UpwindConvectionTerm(coeff=coeff)
        >>> eqn.cacheMatrix()
        >>> eqn.cacheRHSvector()
        >>> eqn.solve(v, dt=1)
        >>> print(numerix.allequal(eqn.matrix.numpyArray,
        ... [[ 1.0,   0,    0,    0,    0,    0,    0,  -1.0,   0,    0,    0,    0,  ],
        ...  [  0,   1.0,   0,    0,    0,    0,    0,   1.0, -1.0,   0,    0,    0,  ],
        ...  [  0,    0,   1.0,   0,    0,    0,    0,    0,   1.0,   0,    0,    0,  ],
        ...  [  0,    0,    0,   1.0,   0,    0,    0,    0,    0,   1.0,   0,    0,  ],
        ...  [  0,    0,    0,    0,   1.0,   0,    0,    0,    0,  -1.0,  1.0,   0, ],
        ...  [  0,    0,    0,    0,    0,   1.0,   0,    0,    0,    0,  -1.0,   0, ],
        ...  [ 2.0,   0,    0,    0,    0,    0,   1.0,   0,    0,    0,    0,    0,  ],
        ...  [-2.0,  2.0,   0,    0,    0,    0,    0,   1.0,   0,    0,    0,    0,  ],
        ...  [  0,  -2.0,   0,    0,    0,    0,    0,    0,   1.0,   0,    0,    0,  ],
        ...  [  0,    0,    0,    0,  -2.0,   0,    0,    0,    0,   1.0,   0,    0,  ],
        ...  [  0,    0,    0,    0,   2.0, -2.0,   0,    0,    0,    0,   1.0,   0,  ],
        ...  [  0,    0,    0,    0,    0,   2.0,   0,    0,    0,    0,    0,   1.0, ]]))
        True
        >>> LHS = CellVariable(mesh=m, rank=1, elementshape=(2,), value=numerix.reshape(numerix.asarray(eqn.matrix * v.value.ravel()), (2, -1))).globalValue.ravel()
        >>> RHS = CellVariable(mesh=m, rank=1, elementshape=(2,), value=numerix.reshape(eqn.RHSvector, (2, -1))).globalValue.ravel()
        >>> print(numerix.allclose(LHS, RHS))
        True

        """


def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


