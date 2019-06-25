from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.solvers.solver import Solver
from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix

class _PysparseMatrixSolver(Solver):

    """
    A class consolidating methods for solver packages which use
    `_PysparseMeshMatrix` for their matrix class.

    Subclasses have a `_solve_` method, which is called by `_solve`. Typically,
    `_solve_` returns the new value of `self.var` to `_solve` and solve sets the
    var accordingly.

    A solution function `solveFnc`, usually of the form `solve(A, x, b)`, is
    implemented in most leaf-node child classes.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    solveFnc = None

    @property
    def _matrixClass(self):
        return _PysparseMeshMatrix

    def _solve(self):
        """
        Call `_solve_` for the new value of `self.var`.

        In certain cases, `_solve_` won't return anything, e.g.
        `fipy.solvers.pysparse.linearLUSolver`. In these cases, we preserve the
        value of `self.var.numericValue`.
        """

        if self.var.mesh.communicator.Nproc > 1:
            raise Exception("%ss cannot be used with multiple processors" \
                            % self.__class__)

        array = self.var.numericValue
        newArr = self._solve_(self.matrix, array, self.RHSvector)

        if newArr is not None:
            array = newArr

        factor = self.var.unit.factor

        if factor != 1:
            array /= self.var.unit.factor

        self.var[:] = array
