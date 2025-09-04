__docformat__ = 'restructuredtext'

from .petscPreconditioner import PETScPreconditioner

__all__ = ["DefaultPreconditioner"]

class DefaultPreconditioner(PETScPreconditioner):
    """Apply PETSc's default preconditioning to :class:`~fipy.solvers.petsc.petscSolver.PETScSolver`.

    "The default preconditioner for sparse matrices is `PCILU` or `PCICC` with
    0 fill on one process and block Jacobi (`PCBJACOBI`) with `PCILU` or `PCICC`
    in parallel." [#PETSc_Default_Preconditioner]_

    .. [#PETSc_Default_Preconditioner] https://petsc.org/main/manualpages/PC/PCCreate/#note
    """

    pctype = None

    def _applyToSolver(self, solver, matrix):
        """Leave solver alone.
        """
        pass
