#!/usr/bin/env python


__docformat__ = 'restructuredtext'

from pysparse import itsolvers

from fipy.solvers.pysparse.preconditioners import SsorPreconditioner
from fipy.solvers.pysparse.pysparseSolver import PysparseSolver

__all__ = ["LinearPCGSolver"]

class LinearPCGSolver(PysparseSolver):
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

    def __init__(self, precon=SsorPreconditioner(), *args, **kwargs):
        """
        :Parameters:
          - `precon`: Preconditioner to use
        """
        super(LinearPCGSolver, self).__init__(precon=precon, *args, **kwargs)
        self.solveFnc = itsolvers.pcg

    def _canSolveAsymmetric(self):
        return False
