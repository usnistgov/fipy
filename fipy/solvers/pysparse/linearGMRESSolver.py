from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from pysparse.itsolvers import krylov

from fipy.solvers.pysparse.preconditioners import JacobiPreconditioner
from fipy.solvers.pysparse.pysparseSolver import PysparseSolver

__all__ = ["LinearGMRESSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearGMRESSolver(PysparseSolver):
    """

    The `LinearGMRESSolver` solves a linear system of equations using the
    generalized minimal residual method (GMRES) with Jacobi
    preconditioning. GMRES solves systems with a general non-symmetric
    coefficient matrix.

    The `LinearGMRESSolver` is a wrapper class for the the Pysparse_
    `itsolvers.gmres()` and `precon.jacobi()` methods.

    .. _Pysparse: http://pysparse.sourceforge.net

    """

    def __init__(self, precon=JacobiPreconditioner(), *args, **kwargs):
        """
        Parameters
        ----------
        precon : ~fipy.solvers.pysparse.preconditioners.preconditioner.Preconditioner, optional
        """
        super(LinearGMRESSolver, self).__init__(precon=precon, *args, **kwargs)
        self.solveFnc = krylov.gmres
