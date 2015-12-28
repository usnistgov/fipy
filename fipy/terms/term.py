#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "term.py"
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

import os

from fipy.tools import numerix
from fipy.terms import AbstractBaseClassError
from fipy.terms import SolutionVariableRequiredError

__all__ = ["Term"]

class Term(object):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def __init__(self, coeff=1., var=None):
        """
        Create a `Term`.

        :Parameters:
          - `coeff`: The coefficient for the term. A `CellVariable` or number.
            `FaceVariable` objects are also acceptable for diffusion or convection terms.

        """

        if self.__class__ is Term:
            raise AbstractBaseClassError

        self.coeff = coeff
        self.geomCoeff = None
        self._cacheMatrix = False
        self._matrix = None
        self._cacheRHSvector = False
        self._RHSvector = None
        self.var = var

    def _calcVars(self):
        raise NotImplementedError

    def _checkCoeff(self, var):
        raise NotImplementedError

    def copy(self):
        return self.__class__(self.coeff, var=self.var)

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):
        raise NotImplementedError

    def _buildAndAddMatrices(self, var, SparseMatrix, boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None, buildExplicitIfOther=False):
        raise NotImplementedError

    def _checkVar(self, var):
        raise NotImplementedError

    def _buildCache(self, matrix, RHSvector):
        if self._cacheMatrix:
            self._matrix = matrix
            self._matrix.cache = True
        else:
            self._matrix = None

        if self._cacheRHSvector:
            self._RHSvector = RHSvector
        else:
            self._RHSvector = None

    def _verifyVar(self, var):
        if var is None:
            if self.var is None:
                raise SolutionVariableRequiredError
            else:
                return self.var
        else:
            return var

    @property
    def _buildExplcitIfOther(self):
        raise NotImplementedError

    def _reshapeIDs(self, var, ids):
        raise NotImplementedError

    def _vectorSize(self, var=None):
        if var is None or var.rank != 1:
            return 1
        else:
            return var.shape[0]

    def _getMatrixClass(self, solver, var):
        if self._vectorSize(var) > 1:
            from fipy.matrices.offsetSparseMatrix import OffsetSparseMatrix
            SparseMatrix =  OffsetSparseMatrix(SparseMatrix=solver._matrixClass,
                                               numberOfVariables=self._vectorSize(var),
                                               numberOfEquations=self._vectorSize(var))
        else:
            SparseMatrix = solver._matrixClass

        return SparseMatrix

    def _prepareLinearSystem(self, var, solver, boundaryConditions, dt):
        solver = self.getDefaultSolver(var, solver)

        var = self._verifyVar(var)
        self._checkVar(var)

        if type(boundaryConditions) not in (type(()), type([])):
            boundaryConditions = (boundaryConditions,)

        for bc in boundaryConditions:
            bc._resetBoundaryConditionApplied()

        if 'FIPY_DISPLAY_MATRIX' in os.environ:
            if not hasattr(self, "_viewer"):
                from fipy.viewers.matplotlibViewer.matplotlibSparseMatrixViewer import MatplotlibSparseMatrixViewer
                Term._viewer = MatplotlibSparseMatrixViewer()

        var, matrix, RHSvector = self._buildAndAddMatrices(var,
                                                           self._getMatrixClass(solver, var),
                                                           boundaryConditions=boundaryConditions,
                                                           dt=dt,
                                                           transientGeomCoeff=self._getTransientGeomCoeff(var),
                                                           diffusionGeomCoeff=self._getDiffusionGeomCoeff(var),
                                                           buildExplicitIfOther=self._buildExplcitIfOther)

        self._buildCache(matrix, RHSvector)

        solver._storeMatrix(var=var, matrix=matrix, RHSvector=RHSvector)

        if 'FIPY_DISPLAY_MATRIX' in os.environ:
            if var is None:
                name = ""
            else:
                if not hasattr(var, "name"):
                    name = ""
                else:
                    name = var.name
            self._viewer.title = r"%s %s" % (name, repr(self))
            from fipy.variables.coupledCellVariable import _CoupledCellVariable
            if isinstance(solver.RHSvector, _CoupledCellVariable):
                RHSvector = solver.RHSvector.globalValue
            else:
                RHSvector = solver.RHSvector
            self._viewer.plot(matrix=solver.matrix, RHSvector=RHSvector)
            from fipy import raw_input
            raw_input()

        return solver

    def solve(self, var=None, solver=None, boundaryConditions=(), dt=None):
        r"""
        Builds and solves the `Term`'s linear system once. This method
        does not return the residual. It should be used when the
        residual is not required.

        :Parameters:

           - `var`: The variable to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations. Defaults to `LinearPCGSolver` for Pysparse and `LinearLUSolver` for Trilinos.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.

        """

        solver = self._prepareLinearSystem(var, solver, boundaryConditions, dt)

        solver._solve()

    def sweep(self, var=None, solver=None, boundaryConditions=(), dt=None, underRelaxation=None, residualFn=None, cacheResidual=False, cacheError=False):
        r"""
        Builds and solves the `Term`'s linear system once. This method
        also recalculates and returns the residual as well as applying
        under-relaxation.

        :Parameters:

           - `var`: The variable to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations. Defaults to `LinearPCGSolver` for Pysparse and `LinearLUSolver` for Trilinos.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.
           - `underRelaxation`: Usually a value between `0` and `1` or `None` in the case of no under-relaxation
           - `residualFn`: A function that takes var, matrix, and RHSvector arguments, used to customize the residual calculation.
           - `cacheResidual`: If `True`, calculate and store the residual vector
              :math:`\vec{r}=\mathsf{L}\vec{x} - \vec{b}` in the `residualVector` member of `Term`
           - `cacheError`: If `True`, use the residual vector :math:`\vec{r}`
              to solve :math:`\mathsf{L}\vec{e}=\vec{r}` for the error vector :math:`\vec{e}`
              and store it in the `errorVector` member of `Term`

        """
        solver = self._prepareLinearSystem(var=var, solver=solver, boundaryConditions=boundaryConditions, dt=dt)
        solver._applyUnderRelaxation(underRelaxation=underRelaxation)
        residual = solver._calcResidual(residualFn=residualFn)

        if cacheResidual or cacheError:
            self.residualVector = solver._calcResidualVector(residualFn=residualFn)

        if cacheError:
            self.errorVector = solver.var.copy()
            var_tmp = solver.var
            RHS_tmp = solver.RHSvector
            solver._storeMatrix(var=self.errorVector, matrix=solver.matrix, RHSvector=self.residualVector)
            solver._solve()
            solver._storeMatrix(var=var_tmp, matrix=solver.matrix, RHSvector=RHS_tmp)

        if not cacheResidual:
            self.residualVector = None

        solver._solve()

        return residual

    def justResidualVector(self, var=None, solver=None, boundaryConditions=(), dt=None, underRelaxation=None, residualFn=None):
        r"""
        Builds the `Term`'s linear system once. This method
        also recalculates and returns the residual as well as applying
        under-relaxation.

        :Parameters:

           - `var`: The variable to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations. Defaults to `LinearPCGSolver` for Pysparse and `LinearLUSolver` for Trilinos.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.
           - `underRelaxation`: Usually a value between `0` and `1` or `None` in the case of no under-relaxation
           - `residualFn`: A function that takes var, matrix, and RHSvector arguments used to customize the residual calculation.

        `justResidualVector` returns the overlapping local value in parallel (not the non-overlapping value).

        >>> from fipy import *
        >>> m = Grid1D(nx=10)
        >>> v = CellVariable(mesh=m)
        >>> len(DiffusionTerm().justResidualVector(v)) == m.numberOfCells
        True

        """
        solver = self._prepareLinearSystem(var, solver, boundaryConditions, dt)
        solver._applyUnderRelaxation(underRelaxation)

        return solver._calcResidualVector(residualFn=residualFn)

    def residualVectorAndNorm(self, var=None, solver=None, boundaryConditions=(), dt=None, underRelaxation=None, residualFn=None):
        r"""
        Builds the `Term`'s linear system once. This method
        also recalculates and returns the residual as well as applying
        under-relaxation.

        :Parameters:

           - `var`: The variable to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations. Defaults to `LinearPCGSolver` for Pysparse and `LinearLUSolver` for Trilinos.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.
           - `underRelaxation`: Usually a value between `0` and `1` or `None` in the case of no under-relaxation
           - `residualFn`: A function that takes var, matrix, and RHSvector arguments used to customize the residual calculation.

        """
        vector = self.justResidualVector(var=var, solver=solver, boundaryConditions=boundaryConditions, dt=dt,
                                         underRelaxation=underRelaxation, residualFn=residualFn)

        L2norm = numerix.L2norm(vector)

        return vector, L2norm

    def justErrorVector(self, var=None, solver=None, boundaryConditions=(), dt=1., underRelaxation=None, residualFn=None):
        r"""
        Builds the `Term`'s linear system once. This method
        also recalculates and returns the error as well as applying
        under-relaxation.

        :Parameters:

           - `var`: The variable to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations. Defaults to `LinearPCGSolver` for Pysparse and `LinearLUSolver` for Trilinos.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.
           - `underRelaxation`: Usually a value between `0` and `1` or `None` in the case of no under-relaxation
           - `residualFn`: A function that takes var, matrix, and RHSvector arguments used to customize the residual calculation.

        `justErrorVector` returns the overlapping local value in parallel (not the non-overlapping value).

        >>> from fipy.solvers import DummySolver
        >>> from fipy import *
        >>> m = Grid1D(nx=10)
        >>> v = CellVariable(mesh=m)
        >>> len(DiffusionTerm().justErrorVector(v, solver=DummySolver())) == m.numberOfCells
        True

        """
        solver = self._prepareLinearSystem(var, solver, boundaryConditions, dt)
        solver._applyUnderRelaxation(underRelaxation)
        residualVector = solver._calcResidualVector(residualFn=residualFn)

        errorVector = solver.var.copy()
        solver._storeMatrix(var=errorVector, matrix=solver.matrix, RHSvector=residualVector)
        solver._solve()

        return errorVector

    def cacheMatrix(self):
        r"""
        Informs `solve()` and `sweep()` to cache their matrix so
        that `getMatrix()` can return the matrix.
        """
        self._cacheMatrix = True

    @property
    def matrix(self):
        r"""
        Return the matrix caculated in `solve()` or `sweep()`. The
        cacheMatrix() method should be called before `solve()` or
        `sweep()` to cache the matrix.
        """
        if not self._cacheMatrix:
            import warnings
            warnings.warn("""cacheMatrix() should be called followed by sweep() or solve()
                          to calculate the matrix. None will be returned on this occasion.""",
                          UserWarning, stacklevel=2)
            self.cacheMatrix()

        return self._matrix

    def cacheRHSvector(self):
        r"""
        Informs `solve()` and `sweep()` to cache their right hand side
        vector so that `getRHSvector()` can return it.
        """
        self._cacheRHSvector = True

    @property
    def RHSvector(self):
        r"""
        Return the RHS vector caculated in `solve()` or `sweep()`. The
        cacheRHSvector() method should be called before `solve()` or
        `sweep()` to cache the vector.
        """
        if not self._cacheRHSvector:
            import warnings
            warnings.warn("""getRHSvector should be called after cacheRHSvector() and either sweep() or solve()
                          to calculate the RHS vector. None will be returned on this occasion.""",
                          UserWarning, stacklevel=2)
            self.cacheRHSvector()

        return self._RHSvector

    def _getDefaultSolver(self, var, solver, *args, **kwargs):
        return NotImplementedError

    def getDefaultSolver(self, var=None, solver=None, *args, **kwargs):
        from fipy.solvers import DefaultSolver
        return solver or self._getDefaultSolver(var, solver, *args, **kwargs) or DefaultSolver(*args, **kwargs)

    def __add__(self, other):
        if isinstance(other, (int, float)) and other == 0:
            return self
        else:
            from fipy.terms.binaryTerm import _BinaryTerm
            return _BinaryTerm(self, other)

    __radd__ = __add__

    def __neg__(self):
        raise NotImplementedError

    def __pos__(self):
        return self

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __eq__(self, other):
        return self - other

    __hash__ = object.__hash__

    def __mul__(self, other):
        raise NotImplementedError

    __rmul__ = __mul__

    def __truediv__(self, other):
        return (1 / other) * self

    __div__ = __truediv__

    def __and__(self, other):
        if isinstance(other, Term):
            from fipy.terms.coupledBinaryTerm import _CoupledBinaryTerm
            return _CoupledBinaryTerm(self, other)
        elif other == 0:
            return self
        else:
            raise Exception, "Can only couple Term objects."

    __rand__ = __and__

    def __repr__(self):
        raise NotImplementedError

    def _calcGeomCoeff(self, var):
        raise NotImplementedError

    def _getGeomCoeff(self, var):
        if self.geomCoeff is None:
            self.geomCoeff = self._calcGeomCoeff(var)
            if self.geomCoeff is not None:
                self.geomCoeff.dontCacheMe()

        return self.geomCoeff

    def _getWeight(self, var, transientGeomCoeff=None, diffusionGeomCoeff=None):
        raise NotImplementedError

    def _getDiagonalSign(self, transientGeomCoeff=None, diffusionGeomCoeff=None):
        raise NotImplementedError

    def _getDiffusionGeomCoeff(self, var):
        return None

    def _getTransientGeomCoeff(self, var):
        return None

    def _getNormals(self, mesh):
        raise NotImplementedError

    def _getOldAdjacentValues(self, oldArray, id1, id2, dt):
        raise NotImplementedError

    def _getDifferences(self, adjacentValues, cellValues, oldArray, cellToCellIDs, mesh):
        raise NotImplementedError

    def _alpha(self, P):
        raise NotImplementedError

    def _checkDt(self, dt):
        raise NotImplementedError

    def _test(self):
        """
        Test stuff.

        >>> from fipy import *
        >>> L = 1.
        >>> nx = 100
        >>> m = Grid1D(nx=nx, dx=L / nx)
        >>> v = CellVariable(mesh=m, value=1.)
        >>> eqn = DiffusionTerm() - v

        >>> v.constrain(0.,  m.facesLeft)
        >>> v.constrain(1.,  m.facesRight)

        >>> res = 1.
        >>> sweep = 0
        >>> while res > 1e-8 and sweep < 100:
        ...     res = eqn.sweep(v)
        ...     sweep += 1
        >>> x = m.cellCenters[0]
        >>> answer = (numerix.exp(x) - numerix.exp(-x)) / (numerix.exp(L) - numerix.exp(-L))
        >>> print numerix.allclose(v, answer, rtol=2e-5)
        True

        >>> v.setValue(0.)
        >>> eqn = DiffusionTerm(0.2) * 5. - 5. * ImplicitSourceTerm(0.2)
        >>> eqn.solve(v)
        >>> print numerix.allclose(v, answer, rtol=2e-5)
        True

        >>> v.setValue(0.)
        >>> eqn = 2. * (DiffusionTerm(1.) - ImplicitSourceTerm(.5)) - DiffusionTerm(1.)
        >>> eqn.solve(v)
        >>> print numerix.allclose(v, answer, rtol=2e-5)
        True

        >>> from fipy import Grid1D, CellVariable, DiffusionTerm, TransientTerm
        >>> mesh = Grid1D(nx=3)
        >>> A = CellVariable(mesh=mesh, name="A", value=[1., 2., 3.])
        >>> B = CellVariable(mesh=mesh, name="B", value=[3., 4., 5.])
        >>> C = CellVariable(mesh=mesh, name="C")

        >>> eq = DiffusionTerm(coeff=1.)
        >>> solver = eq._prepareLinearSystem(var=None, solver=None, boundaryConditions=(), dt=1.) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        ...
        SolutionVariableRequiredError: The solution variable needs to be specified.
        >>> solver = eq._prepareLinearSystem(var=A, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = solver.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[-1, 1, 0],
        ...                                      [ 1,-2, 1],
        ...                                      [ 0, 1,-1]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(solver.RHSvector, [0, 0, 0]) # doctest: +PROCESSOR_0
        True

        >>> eq = DiffusionTerm(coeff=1., var=A)
        >>> solver = eq._prepareLinearSystem(var=None, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = solver.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[-1, 1, 0],
        ...                                      [ 1,-2, 1],
        ...                                      [ 0, 1,-1]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(solver.RHSvector, [0, 0, 0]) # doctest: +PROCESSOR_0
        True
        >>> solver = eq._prepareLinearSystem(var=B, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = solver.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[ 0, 0, 0],
        ...                                      [ 0, 0, 0],
        ...                                      [ 0, 0, 0]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(solver.RHSvector, [0, 0, 0]) # doctest: +PROCESSOR_0
        True

        >>> eq = TransientTerm(coeff=1.) == DiffusionTerm(coeff=1.)
        >>> solver = eq._prepareLinearSystem(var=None, solver=None, boundaryConditions=(), dt=1.) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        ...
        SolutionVariableRequiredError: The solution variable needs to be specified.
        >>> solver = eq._prepareLinearSystem(var=A, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = solver.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[ 2,-1, 0],
        ...                                      [-1, 3,-1],
        ...                                      [ 0,-1, 2]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(solver.RHSvector, [1, 2, 3]) # doctest: +PROCESSOR_0
        True

        >>> eq = TransientTerm(coeff=1., var=A) == DiffusionTerm(coeff=1., var=A)
        >>> solver = eq._prepareLinearSystem(var=None, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = solver.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[ 2,-1, 0],
        ...                                      [-1, 3,-1],
        ...                                      [ 0,-1, 2]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(solver.RHSvector, [1, 2, 3]) # doctest: +PROCESSOR_0
        True
        >>> solver = eq._prepareLinearSystem(var=A, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = solver.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[ 2,-1, 0],
        ...                                      [-1, 3,-1],
        ...                                      [ 0,-1, 2]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(solver.RHSvector, [1, 2, 3]) # doctest: +PROCESSOR_0
        True
        >>> solver = eq._prepareLinearSystem(var=B, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = solver.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[ 0, 0, 0],
        ...                                      [ 0, 0, 0],
        ...                                      [ 0, 0, 0]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(solver.RHSvector, [1, 0, -1]) # doctest: +PROCESSOR_0
        True

        >>> eq = TransientTerm(coeff=1., var=A) == DiffusionTerm(coeff=1., var=B)
        >>> print eq
        (TransientTerm(coeff=1.0, var=A) + DiffusionTerm(coeff=[-1.0], var=B))
        >>> A in set(eq._vars) and B in set(eq._vars) ## _getVars() is unordered for _BinaryTerm's.
        True
        >>> print (eq.term, eq.other)
        (TransientTerm(coeff=1.0, var=A), DiffusionTerm(coeff=[-1.0], var=B))
        >>> res = eq.justResidualVector(boundaryConditions=(), dt=1.)
        >>> print numerix.allequal(res, [-1, 0, 1])  # doctest: +PROCESSOR_0
        True
        >>> solver = eq._prepareLinearSystem(var=A, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = solver.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[1, 0, 0],
        ...                                      [0, 1, 0],
        ...                                      [0, 0, 1]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(solver.RHSvector, [2, 2, 2]) # doctest: +PROCESSOR_0
        True
        >>> solver = eq._prepareLinearSystem(var=B, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = solver.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[ 1,-1, 0],
        ...                                      [-1, 2,-1],
        ...                                      [ 0,-1, 1]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(solver.RHSvector, [0, 0, 0]) # doctest: +PROCESSOR_0
        True
        >>> solver = eq._prepareLinearSystem(var=C, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = solver.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[ 0, 0, 0],
        ...                                      [ 0, 0, 0],
        ...                                      [ 0, 0, 0]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(solver.RHSvector, [1, 0, -1]) # doctest: +PROCESSOR_0
        True

        >>> eq = TransientTerm(coeff=1.) == DiffusionTerm(coeff=1., var=B) + 10.  # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        ExplicitVariableError: Terms with explicit Variables cannot mix with Terms with implicit Variables.
        >>> eq = DiffusionTerm(coeff=1., var=B) + 10. == 0
        >>> print eq
        (DiffusionTerm(coeff=[1.0], var=B) + 10.0)
        >>> print eq._vars
        [B]
        >>> print (eq.term, eq.other)
        (DiffusionTerm(coeff=[1.0], var=B), 10.0)
        >>> solver = eq._prepareLinearSystem(var=B, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = solver.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[-1, 1, 0],
        ...                                      [1, -2, 1],
        ...                                      [0, 1, -1]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(solver.RHSvector, [-10, -10, -10])  # doctest: +PROCESSOR_0
        True

        >>> from fipy.solvers import DummySolver
        >>> eq.solve(var=B, solver=DummySolver())

        >>> m = Grid1D(nx=2)
        >>> A = CellVariable(mesh=m, name='A')
        >>> B = CellVariable(mesh=m, name='B')
        >>> C = CellVariable(mesh=m, name='C')
        >>> DiffusionTerm().solve() # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableRequiredError: The solution variable needs to be specified.
        >>> DiffusionTerm().solve(A, solver=DummySolver())
        >>> DiffusionTerm(var=A).solve(A, solver=DummySolver())
        >>> (DiffusionTerm(var=A) + DiffusionTerm()) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        ExplicitVariableError: Terms with explicit Variables cannot mix with Terms with implicit Variables.
        >>> (DiffusionTerm(var=A) + DiffusionTerm(var=B)).solve(solver=DummySolver())
        >>> (DiffusionTerm(var=A) + DiffusionTerm(var=B)).solve(A, solver=DummySolver())
        >>> DiffusionTerm() & DiffusionTerm() # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableNumberError: Different number of solution variables and equations.
        >>> DiffusionTerm(var=A) & DiffusionTerm() # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        ExplicitVariableError: Terms with explicit Variables cannot mix with Terms with implicit Variables.

        >>> A = CellVariable(mesh=m, name='A', value=1)
        >>> B = CellVariable(mesh=m, name='B')
        >>> C = CellVariable(mesh=m, name='C')

        >>> eq = (DiffusionTerm(coeff=1., var=A)) & (DiffusionTerm(coeff=2., var=B))
        >>> eq.cacheMatrix()
        >>> eq.cacheRHSvector()
        >>> eq.solve(solver=DummySolver())
        >>> numpyMatrix = eq.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[-1, 1, 0, 0],
        ...                                      [1, -1, 0, 0],
        ...                                      [0, 0, -2, 2],
        ...                                      [0, 0, 2, -2]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print eq.RHSvector.globalValue
        [ 0.  0.  0.  0.]
        >>> print eq._vars
        [A, B]
        >>> solver = eq._prepareLinearSystem(var=None, solver=None, boundaryConditions=(), dt=1.)
        >>> numpyMatrix = eq.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[-1, 1, 0, 0],
        ...                                      [1, -1, 0, 0],
        ...                                      [0, 0, -2, 2],
        ...                                      [0, 0, 2, -2]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print eq.RHSvector.globalValue
        [ 0.  0.  0.  0.]
        >>> solver = eq._prepareLinearSystem(var=A, solver=None, boundaryConditions=(), dt=1.) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableNumberError: The solution variable should not be specified.
        >>> solver = eq._prepareLinearSystem(var=C, solver=None, boundaryConditions=(), dt=1.) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableNumberError: The solution variable should not be specified.

        >>> DiffusionTerm(var=A) & DiffusionTerm(var=A) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableNumberError: Different number of solution variables and equations.
        >>> DiffusionTerm() & DiffusionTerm() # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableNumberError: Different number of solution variables and equations.
        >>> (DiffusionTerm(var=A) & DiffusionTerm(var=B)).solve(A) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableNumberError: The solution variable should not be specified.
        >>> DiffusionTerm(var=A) & DiffusionTerm(var=B) & DiffusionTerm(var=B) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableNumberError: Different number of solution variables and equations.
        >>> (DiffusionTerm(var=A) & DiffusionTerm(var=B) & DiffusionTerm(var=C)).solve(solver=DummySolver())
        >>> (DiffusionTerm(var=A) & DiffusionTerm(var=B) & DiffusionTerm(var=C)).solve(A) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableNumberError: The solution variable should not be specified.
        >>> (DiffusionTerm(var=A) & (DiffusionTerm(var=B) + DiffusionTerm(var=C))).solve(A) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableNumberError: The solution variable should not be specified.
        >>> eq = (DiffusionTerm(coeff=1., var=A) + DiffusionTerm(coeff=2., var=B)) & (DiffusionTerm(coeff=2., var=B) + DiffusionTerm(coeff=3., var=C)) & (DiffusionTerm(coeff=3., var=C) + DiffusionTerm(coeff=1., var=A))
        >>> eq.cacheMatrix()
        >>> eq.cacheRHSvector()
        >>> eq.solve(solver=DummySolver())
        >>> numpyMatrix = eq.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[-1, 1, -2, 2, 0, 0],
        ...                                      [1, -1, 2, -2, 0, 0],
        ...                                      [0, 0, -2, 2, -3, 3],
        ...                                      [0, 0, 2, -2, 3, -3],
        ...                                      [-1, 1, 0, 0, -3, 3],
        ...                                      [1, -1, 0, 0, 3, -3]])
        ... # doctest: +PROCESSOR_0
        True
        >>> print eq.RHSvector.globalValue
        [ 0.  0.  0.  0.  0.  0.]
        >>> print eq._vars
        [A, B, C]
        >>> eq = DiffusionTerm(var=A)
        >>> print (0 & eq) is eq
        True
        >>> eq = 0
        >>> eq &= DiffusionTerm(coeff=1., var=A)
        >>> eq &= DiffusionTerm(coeff=2., var=B)
        >>> eq.cacheMatrix()
        >>> eq.cacheRHSvector()
        >>> eq.solve(solver=DummySolver())
        >>> numpyMatrix = eq.matrix.numpyArray
        >>> print numerix.allequal(numpyMatrix, [[-1, 1, 0, 0],
        ...                                      [1, -1, 0, 0],
        ...                                      [0, 0, -2, 2],
        ...                                      [0, 0, 2, -2]])
        ... # doctest: +PROCESSOR_0
        True

        Test for ticket:658, the test should run without an error.

        >>> m = Grid1D(nx=3)
        >>> v = CellVariable(mesh=m, elementshape=(2,))
        >>> v.constrain([[0], [1]], m.facesLeft)
        >>> v.constrain([[1], [0]], m.facesRight)
        >>> eqn = TransientTerm() == DiffusionTerm([[[0.01, -1],[1, 0.01]]])
        >>> res = eqn.sweep(var=v, dt=1.)

        """

class __Term(Term):
    """
    Dummy subclass for tests
    """
    pass

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
