__docformat__ = 'restructuredtext'

from .petscPreconditioner import PETScPreconditioner

__all__ = ["HYPREPreconditioner"]

class HYPREPreconditioner(PETScPreconditioner):
    """
    HYPRE preconditioner for :class:`~fipy.solvers.petsc.petscSolver.PETScSolver`.

    """

    pctype = "hypre"
