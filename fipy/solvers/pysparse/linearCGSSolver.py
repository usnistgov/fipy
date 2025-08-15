from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

import sys
import warnings

from pysparse.itsolvers import krylov

from .linearRHSSolver import LinearRHSSolver

__all__ = ["LinearCGSSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearCGSSolver(LinearRHSSolver):

    """Interface to the conjugate gradient squared (:term:`CGS`) solver in
    :ref:`Pysparse`.
    """

    solveFnc = staticmethod(krylov.cgs)

    def __init__(self, tolerance="default", criterion="default",
                 iterations="default", precon="default"):
        """
        Create a `LinearCGSSolver` object.

        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'unscaled', 'RHS', 'matrix', 'initial', 'legacy'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.pysparse.preconditioners.pysparsePreconditioner.PysparsePreconditioner
            Preconditioner to use.
        """
        warnings.warn("The Pysparse CGS solver may return incorrect results for some matrices", UserWarning)
        super(LinearCGSSolver, self).__init__(tolerance=tolerance, criterion=criterion,
                                              iterations=iterations, precon=precon)
