from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

def _UnaryOperatorVariable(operatorClass=None):
    """
    Test `binOp` pickling

    >>> from fipy import Grid1D, FaceVariable, CellVariable, dump, Variable
    >>> import os, sys
    >>> m = Grid1D()
    >>> vs = (CellVariable(mesh=m, value=2.), FaceVariable(mesh=m, value=3.), Variable(4))
    >>> tmp = []
    >>> for v in vs:
    ...     (f, n) = dump.write(-v)
    ...     tmp += [dump.read(n)]
    ...     if f is not None:
    ...         os.close(f)
    ...         os.remove(n)
    >>> for v in tmp:
    ...     print(v.__class__)
    <class 'fipy.variables.cellVariable.CellVariable'>
    <class 'fipy.variables.faceVariable.FaceVariable'>
    <class 'fipy.variables.variable.Variable'>
    >>> print(tmp[0].allclose(-2.))
    True
    >>> print(tmp[1].allclose(-3.))
    True
    >>> print(tmp[2].allclose(-4.))
    True
    """

    class unOp(operatorClass):
        def _calcValue_(self):
            value = self.op(self.var[0].value)
            if self.return_scalar:
                value = value[()]

            return value

        @property
        def unit(self):
            assert(hasattr(self, "_unit") == True)
            if self._unit is None:
                try:
                    var = self._varProxy
                    return self._extractUnit(self.op(var[0]))
                except:
                    return self._extractUnit(self._calcValue())
            else:
                return self._unit

    return unOp

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()

