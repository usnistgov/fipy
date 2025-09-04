__docformat__ = 'restructuredtext'

from .petscPreconditioner import PETScPreconditioner

__all__ = ["LUPreconditioner"]

class LUPreconditioner(PETScPreconditioner):
    """LU preconditioner for :class:`~fipy.solvers.petsc.petscSolver.PETScSolver`.

    """

    pctype = "lu"
