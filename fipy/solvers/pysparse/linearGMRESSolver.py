


__docformat__ = 'restructuredtext'

from fipy.solvers.pysparse.preconditioners import JacobiPreconditioner
from pysparse import itsolvers

from fipy.solvers.pysparse.pysparseSolver import PysparseSolver

__all__ = ["LinearGMRESSolver"]

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
        :Parameters:
          - `precon`: Preconditioner to use
        """
        super(LinearGMRESSolver, self).__init__(precon=precon, *args, **kwargs)
        self.solveFnc = itsolvers.gmres
