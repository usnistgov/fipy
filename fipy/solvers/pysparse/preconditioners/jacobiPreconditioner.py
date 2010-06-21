from fipy.solvers.pysparse.preconditioners.preconditioner import Preconditioner
from pysparse import precon

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

