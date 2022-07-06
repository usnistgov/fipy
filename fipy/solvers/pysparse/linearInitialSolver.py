from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from .pysparseSolver import PysparseSolver

class LinearInitialSolver(PysparseSolver):
    """Wrapper for solvers that normalize the residual by the initial value.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    criteria = {
        "default": None,
        "initial": None
    }

    def __init__(self, tolerance=1e-10, criterion="default",
                 iterations=1000, precon=None):
        """
        Create a `Solver` object.

        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'initial'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.pysparse.preconditioners.preconditioner.Preconditioner
            Preconditioner to use.
        """
        super(LinearInitialSolver, self).__init__(tolerance=tolerance, criterion=criterion,
                                                  iterations=iterations, precon=precon)
