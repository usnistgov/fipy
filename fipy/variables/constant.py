from __future__ import unicode_literals
from builtins import str
__all__ = []

from fipy.variables.variable import Variable

class _Constant(Variable):
    def __repr__(self):
        return str(self)
