from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix

class _CoupledCellVariable(object):
    def __init__(self, vars):
        self.vars = vars

    @property
    def shape(self):
        return (len(self.vars) * self.mesh.numberOfCells,)

    @property
    def rank(self):
        return self.vars[0].rank

    def __repr__(self):
        return "(" + ", ".join([repr(var) for var in self.vars]) + ")"

    @property
    def mesh(self):
        meshes = list(set([var.mesh for var in self.vars]))

        if len(meshes) == 0:
            raise Exception("There are no Meshes defined")
        elif len(meshes) > 1:
            raise Exception("All Variables must be defined on the same Mesh")
        else:
            return meshes[0]

    def __getitem__(self, index):
        return numerix.concatenate([numerix.array(var[index]) for var in self.vars])

    def __setitem__(self, index, value):
        N = self.mesh.numberOfCells
        for i, var in enumerate(self.vars):
            if numerix.shape(value) == ():
                var[index] = value
            else:
                var[index] = value[i * N:(i + 1) * N]

    def _getValue(self):
        return numerix.concatenate([numerix.array(var.value) for var in self.vars])

    def _setValue(self, value):
        self[:] = value

    value = property(_getValue, _setValue)

    @property
    def globalValue(self):
        return numerix.concatenate([numerix.array(var.globalValue) for var in self.vars])

    @property
    def numericValue(self):
        return numerix.concatenate([var.numericValue for var in self.vars])

    @property
    def unit(self):
        from fipy.tools.dimensions import physicalField
        return physicalField._unity

    def __array__(self, dtype=None, copy=None):
        """
        Attempt to convert the `_CoupledCellVariable` to a numerix `array` object

        >>> from fipy import *
        >>> mesh = Grid1D(nx=2)
        >>> v1 = CellVariable(mesh=mesh, value=[2, 3])
        >>> v2 = CellVariable(mesh=mesh, value=[4, 5])
        >>> v = _CoupledCellVariable(vars=(v1, v2))
        >>> print(numerix.allequal([2, 3, 4, 5], numerix.array(v))) # doctest: +PROCESSOR_0
        True
        >>> v[:] = (6, 7, 8, 9)
        >>> print(v1)
        [6 7]
        >>> print(v2)
        [8 9]
        >>> issubclass(v.dtype.type, numerix.integer)
        True

        """
        if not copy:
            copy = numerix.copy_if_needed

        return numerix.array(self.value, dtype=dtype, copy=copy)

    def __neg__(self):
        return _CoupledCellVariable([-var for var in self.vars])

    def __abs__(self):
        return _CoupledCellVariable([abs(var) for var in self.vars])

    def __iter__(self):
        return iter(self.value)

    @property
    def dtype(self):
        return self.numericValue.dtype

    def ravel(self):
        return self.value.ravel()

    def copy(self):
        return self.__class__(vars=[var.copy() for var in self.vars])

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


