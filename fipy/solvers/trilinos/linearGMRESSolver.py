


__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.trilinosAztecOOSolver import TrilinosAztecOOSolver
from fipy.solvers.trilinos.preconditioners.multilevelDDPreconditioner import MultilevelDDPreconditioner

__all__ = ["LinearGMRESSolver"]

class LinearGMRESSolver(TrilinosAztecOOSolver):

    """
    The `LinearGMRESSolver` is an interface to the GMRES solver in Trilinos,
    using a the `MultilevelDDPreconditioner` by default.

    """

    def __init__(self, tolerance=1e-10, iterations=1000, precon=MultilevelDDPreconditioner()):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use.

        """
        TrilinosAztecOOSolver.__init__(self, tolerance=tolerance,
                                       iterations=iterations, precon=precon)
        self.solver = AztecOO.AZ_gmres
