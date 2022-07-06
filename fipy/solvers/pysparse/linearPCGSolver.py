from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from pysparse.itsolvers import krylov

from .linearRHSSolver import LinearRHSSolver
from .preconditioners import SsorPreconditioner

__all__ = ["LinearPCGSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearPCGSolver(LinearRHSSolver):
    """

    The `LinearPCGSolver` solves a linear system of equations using the
    preconditioned conjugate gradient method (PCG) with symmetric successive
    over-relaxation (SSOR) preconditioning by default. Alternatively,
    Jacobi preconditioning can be specified through `precon`.
    The PCG method solves systems with
    a symmetric positive definite coefficient matrix.

    The `LinearPCGSolver` is a wrapper class for the the Pysparse_
    `itsolvers.pcg()` and `precon.ssor()` methods.

    .. _Pysparse: http://pysparse.sourceforge.net

    """

    solveFnc = staticmethod(krylov.pcg)

    def __init__(self, precon=SsorPreconditioner(), *args, **kwargs):
        """
        Parameters
        ----------
        precon : ~fipy.solvers.pysparse.preconditioners.preconditioner.Preconditioner, optional
        """
        super(LinearPCGSolver, self).__init__(precon=precon, *args, **kwargs)

    def __init__(self, tolerance=1e-10, criterion="default",
                 iterations=1000, precon=SsorPreconditioner()):
        """
        Create a `LinearPCGSolver` object.

        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'RHS'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.pysparse.preconditioners.preconditioner.Preconditioner
            Preconditioner to use
            (default :class:`fipy.solvers.pysparse.preconditioners.SsorPreconditioner`).
        """
        super(LinearPCGSolver, self).__init__(tolerance=tolerance, criterion=criterion,
                                              iterations=iterations, precon=precon)

    def _canSolveAsymmetric(self):
        return False
