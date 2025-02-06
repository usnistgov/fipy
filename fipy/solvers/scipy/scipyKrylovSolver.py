from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = ["ScipyKrylovSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

import os
import warnings
import scipy

from .scipySolver import ScipySolver
from fipy.tools import numerix
from fipy.tools.version import Version, parse_version

scipy_has_tol = (parse_version(scipy.__version__) < Version("1.12"))

class ScipyKrylovSolver(ScipySolver):
    """
    The base `ScipyKrylovSolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def _countIterations(self, xk):
        self.actualIterations += 1

    def _adaptLegacyTolerance(self, L, x, b):
        return self._adaptRHSTolerance(L, x, b)

    def _adaptUnscaledTolerance(self, L, x, b):
        factor = 1. / self._rhsNorm(L, x, b)
        return (factor, None)

    def _adaptRHSTolerance(self, L, x, b):
        return (1., None)

    def _adaptMatrixTolerance(self, L, x, b):
        factor = self._matrixNorm(L, x, b) / self._rhsNorm(L, x, b)
        return (factor, None)

    def _adaptInitialTolerance(self, L, x, b):
        factor = self._residualNorm(L, x, b) / self._rhsNorm(L, x, b)
        return (factor, None)

    def _solve_(self, L, x, b):
        tolerance_scale, _ = self._adaptTolerance(L, x, b)

        self.actualIterations = 0

        rtol = self.scale_tolerance(self.tolerance, tolerance_scale)

        A = L.matrix

        self._log.debug("BEGIN precondition")

        if self.preconditioner is None:
            M = None
        else:
            M, _ = self.preconditioner._applyToMatrix(A)

        self._log.debug("END precondition")
        self._log.debug("BEGIN solve")

        if scipy_has_tol:
            tolerance = dict(tol=rtol)
        else:
            tolerance = dict(rtol=rtol)

        x, info = self.solveFnc(L, b, x,
                                atol=self.absolute_tolerance,
                                maxiter=self.iterations,
                                M=M,
                                callback=self._countIterations,
                                **tolerance)

        self._log.debug("END solve")

        self._setConvergence(suite="scipy",
                             code=numerix.sign(info),
                             actual_code=info,
                             iterations=self.actualIterations,
                             tolerance_scale=tolerance_scale,
                             residual=self._residualNorm(L, x, b))

        self.convergence.warn()

        return x
