__docformat__ = 'restructuredtext'

from .petscPreconditioner import PETScPreconditioner

__all__ = ["SSORPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class SSORPreconditioner(PETScPreconditioner):
    """
    SSOR preconditioner for :class:`~fipy.solvers.petsc.petscSolver.PETScSolver`.

    """

    pctype = "sor"
