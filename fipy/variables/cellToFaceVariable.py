from __future__ import unicode_literals
__all__ = []

from fipy.variables.faceVariable import FaceVariable
from fipy.tools import numerix

class _CellToFaceVariable(FaceVariable):
    def __init__(self, var):
        FaceVariable.__init__(self, mesh=var.mesh, elementshape=var.shape[:-1])
        self.var = self._requires(var)

    def _calcValue(self):
        alpha = self.mesh._faceToCellDistanceRatio
        id1, id2 = self.mesh._adjacentCellIDs

        return self._calcValue_(alpha=alpha, id1=id1, id2=id2)

    def release(self, constraint):
        """Remove `constraint` from `self`

        >>> from fipy import *
        >>> m = Grid1D(nx=3)
        >>> v = CellVariable(mesh=m, value=m.cellCenters[0])
        >>> c0 = Constraint(0., where=m.facesLeft)
        >>> v.constrain(c0)
        >>> c1 = Constraint(3., where=m.facesRight)
        >>> v.faceValue.constrain(c1)
        >>> print(v.faceValue)
        [ 0.  1.  2.  3.]
        >>> v.faceValue.release(constraint=c0)
        >>> print(v.faceValue)
        [ 0.5  1.   2.   3. ]
        >>> v.faceValue.release(constraint=c1)
        >>> print(v.faceValue)
        [ 0.5  1.   2.   2.5]
        """
        try:
            super(FaceVariable, self).constraints.remove(constraint)
        except ValueError:
            self.var.release(constraint=constraint)

    @property
    def constraints(self):
        if hasattr(self.var, "faceConstraints"):
            faceConstraints = self.var.faceConstraints
        else:
            faceConstraints = []

        return super(_CellToFaceVariable, self).constraints + faceConstraints

    @property
    def constraintMask(self):
        if hasattr(self, '_constraintMask'):
            if hasattr(self.var, 'faceConstraints'):
                for constraint in self.var.faceConstraints:
                    self._constraintMask._requires(constraint.where)

        return super(_CellToFaceVariable, self).constraintMask

    def __getstate__(self):
        return dict(var=self.var)

    def __setstate__(self, dict):
        self.__init__(**dict)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()

