from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix

def _BinaryOperatorVariable(operatorClass=None):
    """
    Test `binOp` pickling

        >>> from fipy import Grid1D, FaceVariable, CellVariable, dump, Variable
        >>> import os, sys
        >>> m = Grid1D()
        >>> vs = (CellVariable(mesh=m, value=2.), FaceVariable(mesh=m, value=3.), Variable(4))
        >>> tmp = []
        >>> for v in vs:
        ...     (f, n) = dump.write(v * v)
        ...     tmp += [dump.read(n)]
        ...     if f is not None:
        ...         os.close(f)
        ...         os.remove(n)
        >>> for v in tmp:
        ...     print(v.__class__)
        <class 'fipy.variables.cellVariable.CellVariable'>
        <class 'fipy.variables.faceVariable.FaceVariable'>
        <class 'fipy.variables.variable.Variable'>
        >>> print(tmp[0].allclose(4.))
        True
        >>> print(tmp[1].allclose(9.))
        True
        >>> print(tmp[2].allclose(16))
        True

    """
    # declare a binary operator class with the desired base class
    class binOp(operatorClass):

        def _calcValue_(self):
            from fipy.variables.variable import Variable
            if isinstance(self.var[1], Variable):
                val1 = self.var[1].value
            else:
                if isinstance(self.var[1], type('')):
                    self.var[1] = physicalField.PhysicalField(value=self.var[1])
                val1 = self.var[1]

            value = self.op(self.var[0].value, val1)
            if self.return_scalar:
                value = value[()]

            return value

        @property
        def unit(self):
            if self._unit is None:
                try:
                    var = self._varProxy
                    return self._extractUnit(self.op(var[0], var[1]))
                except:
                    return self._extractUnit(self._calcValue_())
            else:
                return self._unit

        def _getRepresentation(self, style="__repr__", argDict={}, id=id, freshen=False):
            self.id = id
            if (style == "__repr__") and hasattr(self, '_name') and len(self._name) > 0:
                return self._name
            else:
                return "(" + operatorClass._getRepresentation(self, style=style, argDict=argDict, id=id, freshen=freshen) + ")"

    return binOp

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()

