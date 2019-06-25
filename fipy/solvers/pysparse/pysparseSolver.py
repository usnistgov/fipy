from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

import os
from fipy.solvers.pysparseMatrixSolver import _PysparseMatrixSolver

__all__ = ["PysparseSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class PysparseSolver(_PysparseMatrixSolver):
    """
    The base `pysparseSolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, *args, **kwargs):
        if self.__class__ is PysparseSolver:
            raise NotImplementedError("can't instantiate abstract base class")

        super(PysparseSolver, self).__init__(*args, **kwargs)

    def _solve_(self, L, x, b):
        """
        `_solve_` is only for use by solvers which may use
        preconditioning. If you are writing a solver which
        doesn't use preconditioning, this must be overridden.

        Parameters
        ----------
        L : ~fipy.matrices.pysparseMatrix._PysparseMeshMatrix
            Matrix
        x : ndarray
            Solution vector
        b : ndarray
            Right hand side vector
        """

        A = L.matrix

        if self.preconditioner is None:
            P = None
        else:
            P, A = self.preconditioner._applyToMatrix(A)

        info, iter, relres = self.solveFnc(A, b, x, self.tolerance,
                                           self.iterations, P)

        self._raiseWarning(info, iter, relres)

        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            from fipy.tools.debug import PRINT
            PRINT('iterations: %d / %d' % (iter, self.iterations))

            if info < 0:
                PRINT('failure', self._warningList[info].__class__.__name__)
            PRINT('relres:', relres)

    def _solve(self):

        if self.var.mesh.communicator.Nproc > 1:
            raise Exception("Pysparse solvers cannot be used with multiple processors")

        array = self.var.numericValue.ravel()

        from fipy.terms import SolutionVariableNumberError

        if ((self.matrix == 0)
            or (self.matrix.matrix.shape[0] != self.matrix.matrix.shape[1])
            or (self.matrix.matrix.shape[0] != len(array))):

            raise SolutionVariableNumberError

        self._solve_(self.matrix, array, self.RHSvector)
        factor = self.var.unit.factor
        if factor != 1:
            array /= self.var.unit.factor

        self.var[:] = array.reshape(self.var.shape)
