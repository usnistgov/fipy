from __future__ import unicode_literals
__all__ = []

from fipy.variables.faceVariable import FaceVariable
from fipy.tools import numerix

class _FaceGradContributions(FaceVariable):
    """
    Test case

    >>> from fipy import *
    >>> m = Grid2D(nx=3, ny=3)
    >>> x, y = m.cellCenters
    >>> v = CellVariable(mesh=m, elementshape=(3,))
    >>> v[0] = x
    >>> v[1] = y
    >>> v[2] = x**2
    >>> out = _FaceGradContributions(v)
    >>> v0 = CellVariable(mesh=m, value=x)
    >>> v1 = CellVariable(mesh=m, value=y)
    >>> v2 = CellVariable(mesh=m, value=x**2)
    >>> print(numerix.allequal(_FaceGradContributions(v).globalValue.shape, (2, 3, 24)))
    True
    >>> print(numerix.allequal(_FaceGradContributions(v0).globalValue.shape, (2, 24)))
    True
    >>> print(_FaceGradContributions(v0).allclose([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -0.5, 1.,  2.,
    ...                                              2.5, -0.5, 1.,  2.,  2.5, -0.5, 1.,  2.,  2.5],
    ...                                            [-0.5, -1.5, -2.5, 0.5, 1.5, 2.5, 0.5, 1.5, 2.5, 0.5, 1.5, 2.5, 0.,  0.,  0.,
    ...                                              0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0. ]]))
    True
    >>> print((_FaceGradContributions(v0).globalValue == out.globalValue[:, 0]).all())
    True
    >>> print((_FaceGradContributions(v1).globalValue == out.globalValue[:, 1]).all())
    True
    >>> print((_FaceGradContributions(v2).globalValue == out.globalValue[:, 2]).all())
    True

    """

    def __init__(self, var):
        FaceVariable.__init__(self, mesh=var.mesh, elementshape=(var.mesh.dim,) + var.shape[:-1])
        self.var = self._requires(var)

    def _calcValue(self):
        faceValue = self.var.arithmeticFaceValue.numericValue
        return self.mesh._areaProjections[(slice(0, None, None),) + (numerix.newaxis,) * (len(faceValue.shape) - 1) + (slice(0, None, None),)] * faceValue[numerix.newaxis]

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


