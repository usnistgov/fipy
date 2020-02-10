from __future__ import unicode_literals
from pysparse.precon import precon

from fipy.solvers.pysparse.preconditioners.preconditioner import Preconditioner

__all__ = ["SsorPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class SsorPreconditioner(Preconditioner):
    """
    SSOR preconditioner for Pysparse.
    Really just a wrapper class for `pysparse.precon.jacobi`.
    """
    def _applyToMatrix(self, A):
        """
        Returns (preconditioning matrix, resulting matrix)
        """
        A = A.to_sss()
        return precon.ssor(A), A
