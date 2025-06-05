from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.trilinosAztecOOSolver import TrilinosAztecOOSolver
from fipy.solvers.trilinos.preconditioners.jacobiPreconditioner import JacobiPreconditioner

__all__ = ["LinearBicgstabSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearBicgstabSolver(TrilinosAztecOOSolver):

    """
    The `LinearBicgstabSolver` is an interface to the biconjugate gradient
    stabilized solver in Trilinos, using the `JacobiPreconditioner` by default.

    """

    solver = AztecOO.AZ_bicgstab

    DEFAULT_PRECONDITIONER = JacobiPreconditioner
