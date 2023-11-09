from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from .preconditioner import Preconditioner

__all__ = ["DefaultPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class DefaultPreconditioner(Preconditioner):
    """
    Apply PETSc's default preconditioning.

    "The default preconditioner for sparse matrices is PCILU or PCICC with
    0 fill on one process and block Jacobi (PCBJACOBI) with PCILU or PCICC
    in parallel." [#PETSc_Default_Preconditioner]_

    .. [#PETSc_Default_Preconditioner] https://petsc.org/main/manualpages/PC/PCCreate/#note
    """

    pctype = None

    def _applyToSolver(self, solver, matrix):
        """Leave solver alone.
        """
        pass
