from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from pysparse.itsolvers import krylov

from .linearInitialSolver import LinearInitialSolver
from .preconditioners import JacobiPreconditioner

__all__ = ["LinearGMRESSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearGMRESSolver(LinearInitialSolver):
    """

    The `LinearGMRESSolver` solves a linear system of equations using the
    generalized minimal residual method (GMRES) with Jacobi
    preconditioning. GMRES solves systems with a general non-symmetric
    coefficient matrix.

    The `LinearGMRESSolver` is a wrapper class for the the Pysparse_
    `itsolvers.gmres()` and `precon.jacobi()` methods.

    .. _Pysparse: http://pysparse.sourceforge.net

    """

    solveFnc = staticmethod(krylov.gmres)

    def __init__(self, tolerance=1e-10, criterion="default",
                 iterations=1000, precon=JacobiPreconditioner()):
        """
        Create a `LinearGMRESSolver` object.

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
            (default :class:`fipy.solvers.pysparse.preconditioners.JacobiPreconditioner`).
        """
        super(LinearGMRESSolver, self).__init__(tolerance=tolerance, criterion=criterion,
                                                iterations=iterations, precon=precon)
