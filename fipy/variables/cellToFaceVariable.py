#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "cellToFaceVariable.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

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
        >>> print v.faceValue
        [ 0.  1.  2.  3.]
        >>> v.faceValue.release(constraint=c0)
        >>> print v.faceValue
        [ 0.5  1.   2.   3. ]
        >>> v.faceValue.release(constraint=c1)
        >>> print v.faceValue
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
