__docformat__ = 'restructuredtext'

from .petscPreconditioner import PETScPreconditioner

__all__ = ["ILUPreconditioner"]

class ILUPreconditioner(PETScPreconditioner):
    """Incomplete LU preconditioner for :class:`~fipy.solvers.petsc.petscSolver.PETScSolver`.
    """

    pctype = "ilu"
