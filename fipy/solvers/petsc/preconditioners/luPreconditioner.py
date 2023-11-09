from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from .preconditioner import Preconditioner

__all__ = ["LUPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LUPreconditioner(Preconditioner):
    """
    LU Preconditioner for PETSc solvers.

    """

    pctype = "lu"
