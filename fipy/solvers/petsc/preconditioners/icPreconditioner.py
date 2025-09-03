__docformat__ = 'restructuredtext'

from .petscPreconditioner import PETScPreconditioner

__all__ = ["ICPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class ICPreconditioner(PETScPreconditioner):
    """Incomplete Choleski preconditioner for :class:`~fipy.solvers.petsc.petscSolver.PETScSolver`.
    """

    pctype = "icc"
