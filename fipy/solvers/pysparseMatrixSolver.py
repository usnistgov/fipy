from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.solvers.solver import Solver
from fipy.matrices.pysparseMatrix import _PysparseMeshMatrix

class PysparseMatrixSolver(Solver):

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

    @property
    def _matrixClass(self):
        return _PysparseMeshMatrix

    @staticmethod
    def solveFnc(*args, **kwargs):
        raise NotImplementedError()
