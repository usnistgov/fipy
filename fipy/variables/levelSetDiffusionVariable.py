#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "levelSetDiffusionVariable.py"
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

__docformat__ = 'restructuredtext'

__all__ = []

from fipy.variables.cellToFaceVariable import _CellToFaceVariable
from fipy.tools import inline
from fipy.tools import numerix

class _LevelSetDiffusionVariable(_CellToFaceVariable):
    r"""
    This variable sets it's face value to zero if either of the
    surrounding cell values are zero else it uses the value of the
    diffusion coefficient. The diffusion coefficient is given by,

    .. math::

        D = \begin{cases}
            D_c & \text{when $\phi > 0$} \\
            0  & \text{when $\phi \le 0$}

    Here is a simple 1D test case:

    >>> from fipy.meshes import Grid2D
    >>> mesh = Grid2D(dx = 1., nx = 3)
    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(mesh = mesh, value = (-1, 1, 1))
    >>> from fipy.variables.faceVariable import FaceVariable
    >>> answer = FaceVariable(mesh=mesh, value=(0,1,1,0,1,1,0,0,1,1))
    >>> print _LevelSetDiffusionVariable(var, 1).allclose(answer)
    True
    """
    def __init__(self, distanceVariable = None, diffusionCoeff = None):
        """
        Creates a `_LevelSetDiffusionVariable`.

        :Parameters:
          - `distanceVariable` : A `DistanceVariable` object
          - `diffusionCoeff` : The `coeff` value.

        """
        _CellToFaceVariable.__init__(self, distanceVariable)
        self.diffusionCoeff = diffusionCoeff

    if inline.doInline:
        def _calcValue_(self, alpha, id1, id2):
            val = self._array.copy()

            inline._runInline("""
                int ID1 = id1[i];
                int ID2 = id2[i];
                double	cell1 = var[ID1];
                double	cell2 = var[ID2];

                if (cell1 < 0 || cell2 < 0) {
                    val[i] = 0;
                } else {
                    val[i] = diffusionCoeff;
                }
            """,
            var = numerix.array(self.var),
            val = val,
            id1 = id1, id2 = id2,
            diffusionCoeff = self.diffusionCoeff,
            ni = self.mesh.numberOfFaces
            )

            return self._makeValue(value = val)
    else:
        def _calcValue_(self, alpha, id1, id2):
            distance = numerix.array(self.var)
            cell1 = numerix.take(distance, id1)
            cell2 = numerix.take(distance, id2)

            return numerix.where(numerix.logical_or(cell1 < 0, cell2 < 0),
                                 0,
                                 self.diffusionCoeff)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
