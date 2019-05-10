from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext'

__all__ = ["Constraint"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Constraint(object):
    def __init__(self, value, where=None):
        """Object to hold a `Variable` to `value` at `where`

        see :meth:`~fipy.variables.variable.Variable.constrain`
        """
        self.value = value
        self.where = where

    def __repr__(self):
        return "Constraint(value=%s, where=%s)" % (repr(self.value), repr(self.where))
