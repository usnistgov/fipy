__docformat__ = 'restructuredtext'

from fipy.solvers.preconditioner import SolverModifyingPreconditioner

__all__ = ["PETScPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class PETScPreconditioner(SolverModifyingPreconditioner):
    """Base class preconditioners of  for :class:`~fipy.solvers.petsc.petscSolver.PETScSolver`.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def _applyToSolver(self, solver, matrix):
        solver.getPC().setType(self.pctype)
