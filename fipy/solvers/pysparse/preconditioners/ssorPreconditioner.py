from pysparse import precon

from fipy.solvers.pysparse.preconditioners.preconditioner import Preconditioner

__all__ = ["SsorPreconditioner"]

class SsorPreconditioner(Preconditioner):
    """
    SSOR preconditioner for PySparse.
    Really just a wrapper class for pysparse.precon.jacobi.
    """
    def _applyToMatrix(self, A):
        """
        Returns (preconditioning matrix, resulting matrix)
        """
        A = A.to_sss()
        return precon.ssor(A), A
