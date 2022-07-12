# Adapted from https://stackoverflow.com/q/46876951/2019542

from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext'

__all__ = ["JacobiPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

from scipy.sparse import diags
from scipy.sparse.linalg import LinearOperator, spsolve

from .preconditioner import Preconditioner

class JacobiPreconditioner(Preconditioner):
    
    def _applyToMatrix(self, A):
        P = diags(A.diagonal()).tocsc()
        Mx = lambda x: spsolve(P, x)

        return LinearOperator(A.shape, Mx)
