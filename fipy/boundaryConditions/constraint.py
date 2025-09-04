"""Restriction on value of a :class:`~fipy.variables.variable.Variable`
"""
from builtins import object
__docformat__ = 'restructuredtext'

__all__ = ["Constraint"]

class Constraint(object):
    """Holds a :class:`~fipy.variables.variable.Variable` to `value` at `where`

    see :meth:`~fipy.variables.variable.Variable.constrain`
    """

    def __init__(self, value, where=None):
        self.value = value
        self.where = where

    def __repr__(self):
        return "Constraint(value=%s, where=%s)" % (repr(self.value), repr(self.where))
