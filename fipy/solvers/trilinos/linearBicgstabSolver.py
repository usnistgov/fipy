#!/usr/bin/env python


__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.trilinosAztecOOSolver import TrilinosAztecOOSolver
from fipy.solvers.trilinos.preconditioners.jacobiPreconditioner import JacobiPreconditioner

__all__ = ["LinearBicgstabSolver"]

class LinearBicgstabSolver(TrilinosAztecOOSolver):

    """
    The `LinearBicgstabSolver` is an interface to the biconjugate gradient
    stabilized solver in Trilinos, using the `JacobiPreconditioner` by default.

    """

    def __init__(self, tolerance=1e-10, iterations=1000, precon=JacobiPreconditioner()):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use.

        """
        TrilinosAztecOOSolver.__init__(self, tolerance=tolerance,
                                       iterations=iterations, precon=precon)
        self.solver = AztecOO.AZ_bicgstab
