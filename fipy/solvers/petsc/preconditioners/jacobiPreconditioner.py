__docformat__ = 'restructuredtext'

from .petscPreconditioner import PETScPreconditioner

__all__ = ["JacobiPreconditioner"]

class JacobiPreconditioner(PETScPreconditioner):
    """Jacobi preconditioner for :class:`~fipy.solvers.petsc.petscSolver.PETScSolver`.
    """

    pctype = "jacobi"
