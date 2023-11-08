from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from ..pysparseMatrixSolver import _PysparseMatrixSolver
from fipy.tools import numerix

__all__ = ["PysparseSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class PysparseSolver(_PysparseMatrixSolver):
    """
    The base `pysparseSolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def __init__(self, tolerance=1e-10, criterion="default",
                 iterations=1000, precon=None):
        """
        Create a `Solver` object.

        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'unscaled', 'RHS', 'matrix', 'initial', 'legacy'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.pysparse.preconditioners.preconditioner.Preconditioner
            Preconditioner to use.
        """
        if self.__class__ is PysparseSolver:
            raise NotImplementedError("can't instantiate abstract base class")

        super(PysparseSolver, self).__init__(tolerance=tolerance, criterion=criterion,
                                             iterations=iterations, precon=precon)

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

        tolerance_scale, _ = self._adaptTolerance(L, x, b)

        # Pysparse returns the relative residual,
        # which changes depending on which solver is used
        legacy_norm = self._legacyNorm(L, x, b)

        self._log.debug("BEGIN solve")

        info, iter, relres = self.solveFnc(A, b, x,
                                           self.tolerance * tolerance_scale,
                                           self.iterations, P)

        self._log.debug("END solve")

        self._setConvergence(suite="pysparse",
                             code=info,
                             iterations=iter,
                             tolerance_scale=tolerance_scale,
                             residual=relres * legacy_norm)

        self.convergence.warn()

    def _rhsNorm(self, L, x, b):
        return numerix.L2norm(b)

    def _matrixNorm(self, L, x, b):
        return L.matrix.norm('inf')

    def _residualVectorAndNorm(self, L, x, b):
        residualVector = L * x - b

        return residualVector, numerix.L2norm(residualVector)

    @property
    def _Lxb(self):
        return (self.matrix, self.var.ravel(), numerix.array(self.RHSvector))

    def _adaptUnscaledTolerance(self, L, x, b):
        factor = 1. / self._legacyNorm(L, x, b)
        return (factor, None)

    def _adaptRHSTolerance(self, L, x, b):
        factor = self._rhsNorm(L, x, b) / self._legacyNorm(L, x, b)
        return (factor, None)

    def _adaptMatrixTolerance(self, L, x, b):
        factor = self._matrixNorm(L, x, b) / self._legacyNorm(L, x, b)
        return (factor, None)

    def _adaptInitialTolerance(self, L, x, b):
        factor = self._residualNorm(L, x, b) / self._legacyNorm(L, x, b)
        return (factor, None)

    def _adaptLegacyTolerance(self, L, x, b):
        return (1., None)

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
