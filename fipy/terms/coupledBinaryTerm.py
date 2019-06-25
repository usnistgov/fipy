from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.terms.abstractBinaryTerm import _AbstractBinaryTerm
from fipy.variables.coupledCellVariable import _CoupledCellVariable
from fipy.variables.cellVariable import CellVariable
from fipy.tools import numerix
from fipy.terms import SolutionVariableNumberError
from fipy.matrices.offsetSparseMatrix import OffsetSparseMatrix

class _CoupledBinaryTerm(_AbstractBinaryTerm):
    """
    Test to ensure that `_getTransientGeomCoeff` and `_getDiffusionGeomCoeff` return sensible results for coupled equations.

    >>> from fipy import *
    >>> m = Grid1D(nx=1)
    >>> v0 = CellVariable(mesh=m)
    >>> v1 = CellVariable(mesh=m)
    >>> eq0 = TransientTerm(1, var=v0) == DiffusionTerm(2, var=v1)
    >>> eq1 = TransientTerm(3, var=v1) == DiffusionTerm(4, var=v0)
    >>> eq = eq0 & eq1
    >>> print(eq._getTransientGeomCoeff(v0))
    None
    >>> print(eq._getDiffusionGeomCoeff(v1))
    None
    >>> tranCoeff = eq._uncoupledTerms[0]._getTransientGeomCoeff(v0)
    >>> print(numerix.allequal(tranCoeff, [1])) # doctest: +PROCESSOR_0
    True
    >>> diffCoeff = eq._uncoupledTerms[1]._getDiffusionGeomCoeff(v0)
    >>> print(numerix.allequal(diffCoeff, [[-8, -8]]))
    True

    """
    def __init__(self, term, other):
        _AbstractBinaryTerm.__init__(self, term, other)
        if len(self._vars) < len(self._uncoupledTerms):
            raise SolutionVariableNumberError

    @property
    def _uncoupledTerms(self):
        return self.term._uncoupledTerms + self.other._uncoupledTerms

    def _verifyVar(self, var):
        if var is not None:
            raise SolutionVariableNumberError('The solution variable should not be specified.')

        if len(self._vars) != len(self._uncoupledTerms):
            raise SolutionVariableNumberError

        return _AbstractBinaryTerm._verifyVar(self, _CoupledCellVariable(self._vars))

    @property
    def _buildExplcitIfOther(self):
        return False

    def _buildAndAddMatrices(self, var, SparseMatrix,  boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None, buildExplicitIfOther=False):
        """Build matrices of constituent Terms and collect them

        Only called at top-level by `_prepareLinearSystem()`

        """

        from fipy.matrices.offsetSparseMatrix import OffsetSparseMatrix
        SparseMatrix =  OffsetSparseMatrix(SparseMatrix=SparseMatrix,
                                           numberOfVariables=len(self._vars),
                                           numberOfEquations=len(self._uncoupledTerms))
        matrix = SparseMatrix(mesh=var.mesh)
        RHSvectors = []

        for equationIndex, uncoupledTerm in enumerate(self._uncoupledTerms):

            SparseMatrix.equationIndex = equationIndex
            termRHSvector = 0
            termMatrix = SparseMatrix(mesh=var.mesh)

            for varIndex, tmpVar in enumerate(var.vars):

                SparseMatrix.varIndex = varIndex

                tmpVar, tmpMatrix, tmpRHSvector = uncoupledTerm._buildAndAddMatrices(tmpVar,
                                                                                     SparseMatrix,
                                                                                     boundaryConditions=(),
                                                                                     dt=dt,
                                                                                     transientGeomCoeff=uncoupledTerm._getTransientGeomCoeff(tmpVar),
                                                                                     diffusionGeomCoeff=uncoupledTerm._getDiffusionGeomCoeff(tmpVar),
                                                                                     buildExplicitIfOther=buildExplicitIfOther)

                termMatrix += tmpMatrix
                termRHSvector += tmpRHSvector

            uncoupledTerm._buildCache(termMatrix, termRHSvector)
            RHSvectors += [CellVariable(value=termRHSvector, mesh=var.mesh)]
            matrix += termMatrix

        return (var, matrix, _CoupledCellVariable(RHSvectors))

    def __repr__(self):
        return '(' + repr(self.term) + ' & ' + repr(self.other) + ')'

    def _getDefaultSolver(self, var, solver, *args, **kwargs):
        if solver and not solver._canSolveAsymmetric():
            import warnings
            warnings.warn("%s cannot solve asymmetric matrices" % solver)
        from fipy.solvers import DefaultAsymmetricSolver
        return solver or DefaultAsymmetricSolver(*args, **kwargs)

    def _calcVars(self):
        """
        This method returns the equations variables ordered by, transient terms,
        diffusion terms and other terms. Currently, this won't cure all of
        coupled equations ills, but it fixes the majority of issues, mainly with
        preconditioning. Some tests to make sure it works correctly.

        >>> from fipy import *
        >>> m = Grid1D(nx=2)
        >>> v0 = CellVariable(mesh=m, name='v0')
        >>> v1 = CellVariable(mesh=m, name='v1')
        >>> v2 = CellVariable(mesh=m, name='v2')
        >>> ConvectionTerm = CentralDifferenceConvectionTerm
        >>> eq0 = TransientTerm(var=v0) + VanLeerConvectionTerm(var=v0) + DiffusionTerm(var=v1)
        >>> eq1 = TransientTerm(var=v2) + ConvectionTerm(var=v2) + DiffusionTerm(var=v2) + ConvectionTerm(var=v1) + ImplicitSourceTerm(var=v1)
        >>> eq2 = ImplicitSourceTerm(var=v1) + 1 + ImplicitSourceTerm(var=v0) + 1 + DiffusionTerm(var=v1)
        >>> print((eq0 & eq1 & eq2)._vars)
        [v0, v2, v1]
        >>> print((eq0 & eq2 & eq1)._vars)
        [v0, v1, v2]
        >>> eq0 =  DiffusionTerm(var=v1) + TransientTerm(var=v0) + VanLeerConvectionTerm(var=v0)
        >>> print((eq0 & eq2 & eq1)._vars)
        [v0, v1, v2]
        >>> print((eq2 & eq0 & eq1)._vars)
        [v1, v0, v2]
        >>> print((eq2 & eq0 & eq1)([v1, v2, v0])._vars)
        [v1, v2, v0]
        >>> print((eq2 & eq0 & eq1)([v1, v2, v0, v2])._vars) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableNumberError: Different number of solution variables and equations.
        >>> print((eq2 & eq0 & eq1)([v1, v2, 1])._vars) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableNumberError: Variable not in previously defined variables for this coupled equation.
        >>> print((eq2 & eq0 & eq1)([v1, v2, v1])._vars) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableNumberError: Different number of solution variables and equations.
        >>> print((eq2 & eq0 & eq1)([v1, v2])._vars) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        SolutionVariableNumberError: Different number of solution variables and equations.

        """

        ## set() is used to force comparison by reference rather than value
        unorderedVars = _AbstractBinaryTerm._calcVars(self)
        uncoupledTerms = self._uncoupledTerms

        if len(unorderedVars) == len(uncoupledTerms):
            unorderedVars = set(unorderedVars)
            orderedVars = [None] * len(uncoupledTerms)

            for fnc in (lambda index, term: term._transientVars,
                        lambda index, term: term._diffusionVars,
                        lambda index, term: list(unorderedVars)):
                for index, term in enumerate(uncoupledTerms):
                    if orderedVars[index] is None:
                        _vars = fnc(index, term)
                        if  _vars != [] and _vars[0] in unorderedVars:
                            orderedVars[index] = _vars[0]
                            unorderedVars.remove(_vars[0])

            return orderedVars
        else:
            ## Constituent _CoupledBinaryTerms don't necessarily have the same
            ## number of equations and variables so ordering is unnecessary.
            return unorderedVars

    def __call__(self, _vars):
        _vars = list(_vars)

        if len(_vars) != len(self._vars) or len(set(_vars)) != len(self._vars):
            raise SolutionVariableNumberError

        for var in _vars:
            if var not in set(self._vars):
                raise SolutionVariableNumberError('Variable not in previously defined variables for this coupled equation.')

        self._internalVars = _vars
        return self

    def _test(self):
        """
        Offset tests

        >>> from fipy import *
        >>> m = Grid1D(nx=3)
        >>> v0 = CellVariable(mesh=m, value=0.)
        >>> v1 = CellVariable(mesh=m, value=1.)
        >>> eq0 = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0) - DiffusionTerm(coeff=2., var=v1)
        >>> eq1 = TransientTerm(var=v1) - DiffusionTerm(coeff=3., var=v0) - DiffusionTerm(coeff=4., var=v1)
        >>> eq = eq0 & eq1
        >>> var, matrix, RHSvector = eq._buildAndAddMatrices(var=eq._verifyVar(None), SparseMatrix=DefaultSolver()._matrixClass, dt=1.)
        >>> print(var.globalValue)
        [ 0.  0.  0.  1.  1.  1.]
        >>> print(RHSvector.globalValue)
        [ 0.  0.  0.  1.  1.  1.]
        >>> print(numerix.allequal(matrix.numpyArray,
        ...                        [[ 2, -1,  0,  2, -2,  0],
        ...                         [-1,  3, -1, -2,  4, -2],
        ...                         [ 0, -1,  2,  0, -2,  2],
        ...                         [ 3, -3,  0,  5, -4,  0],
        ...                         [-3,  6, -3, -4,  9, -4],
        ...                         [ 0, -3,  3,  0, -4,  5]]))
        True

        >>> m = Grid1D(nx=6)
        >>> v0 = CellVariable(mesh=m, value=0.)
        >>> v1 = CellVariable(mesh=m, value=1.)
        >>> eq0 = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0) - DiffusionTerm(coeff=2., var=v1)
        >>> eq1 = TransientTerm(var=v1) - DiffusionTerm(coeff=3., var=v0) - DiffusionTerm(coeff=4., var=v1)
        >>> eq = eq0 & eq1
        >>> var, matrix, RHSvector = eq._buildAndAddMatrices(var=eq._verifyVar(None), SparseMatrix=DefaultSolver()._matrixClass, dt=1.)
        >>> print(var.globalValue)
        [ 0.  0.  0.  0.  0.  0.  1.  1.  1.  1.  1.  1.]
        >>> print(RHSvector.globalValue)
        [ 0.  0.  0.  0.  0.  0.  1.  1.  1.  1.  1.  1.]
        >>> print(numerix.allequal(matrix.numpyArray,
        ...                        [[ 2, -1,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0],
        ...                         [-1,  3, -1,  0,  0,  0, -2,  4, -2,  0,  0,  0],
        ...                         [ 0, -1,  3, -1,  0,  0,  0, -2,  4, -2,  0,  0],
        ...                         [ 0,  0, -1,  3, -1,  0,  0,  0, -2,  4, -2,  0],
        ...                         [ 0,  0,  0, -1,  3, -1,  0,  0,  0, -2,  4, -2],
        ...                         [ 0,  0,  0,  0, -1,  2,  0,  0,  0,  0, -2,  2],
        ...                         [ 3, -3,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0],
        ...                         [-3,  6, -3,  0,  0,  0, -4,  9, -4,  0,  0,  0],
        ...                         [ 0, -3,  6, -3,  0,  0,  0, -4,  9, -4,  0,  0],
        ...                         [ 0,  0, -3,  6, -3,  0,  0,  0, -4,  9, -4,  0],
        ...                         [ 0,  0,  0, -3,  6, -3,  0,  0,  0, -4,  9, -4],
        ...                         [ 0,  0,  0,  0, -3,  3,  0,  0,  0,  0, -4,  5]]))
        True

        >>> m = Grid1D(nx=3)
        >>> v0 = CellVariable(mesh=m, value=0.)
        >>> v1 = CellVariable(mesh=m, value=1.)
        >>> diffTerm = DiffusionTerm(coeff=1., var=v0)
        >>> eq00 = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0)
        >>> eq0 = eq00 - DiffusionTerm(coeff=2., var=v1)
        >>> eq1 = TransientTerm(var=v1) - DiffusionTerm(coeff=3., var=v0) - DiffusionTerm(coeff=4., var=v1)
        >>> eq0.cacheMatrix()
        >>> diffTerm.cacheMatrix()
        >>> (eq0 & eq1).solve(dt=1.)
        >>> print(numerix.allequal(eq0.matrix.numpyArray,
        ...                        [[ 2, -1,  0,  2, -2,  0],
        ...                         [-1,  3, -1, -2,  4, -2],
        ...                         [ 0, -1,  2,  0, -2,  2],
        ...                         [ 0,  0,  0,  0,  0,  0],
        ...                         [ 0,  0,  0,  0,  0,  0],
        ...                         [ 0,  0,  0,  0,  0,  0]]))
        True
        >>> ## This correctly returns None because we lost the handle to the DiffusionTerm when it's negated.
        >>> print(diffTerm.matrix)
        None

        Check `diffusionGeomCoeff` is determined correctly (ticket:329)

        >>> m = Grid1D(nx=3)
        >>> v0 = CellVariable(mesh=m, value=0.)
        >>> v1 = CellVariable(mesh=m, value=1.)
        >>> eq0 = PowerLawConvectionTerm(coeff=100., var=v0) - DiffusionTerm(coeff=1., var=v0)
        >>> eq1 = PowerLawConvectionTerm(coeff=-0.001, var=v1) - DiffusionTerm(coeff=3., var=v1)
        >>> eq = eq0 & eq1
        >>> var, matrix, RHSvector = eq._buildAndAddMatrices(var=eq._verifyVar(None), SparseMatrix=DefaultSolver()._matrixClass)
        >>> print(var.globalValue)
        [ 0.  0.  0.  1.  1.  1.]
        >>> print(RHSvector.globalValue)
        [ 0.  0.  0.  0.  0.  0.]
        >>> print(numerix.allclose(matrix.numpyArray,
        ...                        [[ 100, 1e-15,      0,       0,       0,       0],
        ...                         [-100,   100,  1e-15,       0,       0,       0],
        ...                         [   0,  -100, -1e-15,       0,       0,       0],
        ...                         [   0,     0,      0,  2.9995, -3.0005,       0],
        ...                         [   0,     0,      0, -2.9995,  6.0000, -3.0005],
        ...                         [   0,     0,      0,       0, -2.9995,  3.0005]]))
        True
        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()

