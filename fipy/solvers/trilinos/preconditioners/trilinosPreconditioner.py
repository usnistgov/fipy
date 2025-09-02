from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.solvers.preconditioner import SolverModifyingPreconditioner

__all__ = ["TrilinosPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class TrilinosPreconditioner(SolverModifyingPreconditioner):
    """Base class of preconditioners for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    pass
