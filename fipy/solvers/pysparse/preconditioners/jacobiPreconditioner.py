


from pysparse import precon

from fipy.solvers.pysparse.preconditioners.preconditioner import Preconditioner

__all__ = ["JacobiPreconditioner"]

class JacobiPreconditioner(Preconditioner):
    """
    Jacobi preconditioner for PySparse.
    Really just a wrapper class for pysparse.precon.jacobi.
    """
    def _applyToMatrix(self, A):
        """
        Returns (preconditioning matrix, resulting matrix)
        """
        return precon.jacobi(A), A.to_csr()
