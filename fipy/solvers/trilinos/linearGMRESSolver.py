from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.trilinosAztecOOSolver import TrilinosAztecOOSolver
from fipy.solvers.trilinos.preconditioners.multilevelDDPreconditioner import MultilevelDDPreconditioner

__all__ = ["LinearGMRESSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearGMRESSolver(TrilinosAztecOOSolver):

    """
    The `LinearGMRESSolver` is an interface to the GMRES solver in Trilinos,
    using the `MultilevelDDPreconditioner` by default.

    """

    solver = AztecOO.AZ_gmres

    DEFAULT_PRECONDITIONER = MultilevelDDPreconditioner
