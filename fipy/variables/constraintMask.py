from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix

def _ConstraintMask(var):
    class _ConstraintMaskVariable(var._variableClass):
        def __init__(self, var):
            super(_ConstraintMaskVariable, self).__init__(mesh=var.mesh, rank=0, value=False)
            for constraint in var.constraints:
                self._requires(constraint.where)
            self.var = var

        def _calcValue(self):
            returnMask = numerix.zeros(self.shape, dtype=bool)
            for constraint in self.var.constraints:
                returnMask = returnMask | numerix.asarray(constraint.where)
            return returnMask

    return _ConstraintMaskVariable(var)
