from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.solvers.scipy.linearLUSolver import LinearLUSolver as SciPyLinearLUSolver

__all__ = ["LinearLUSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearLUSolver(SciPyLinearLUSolver):
    pass
