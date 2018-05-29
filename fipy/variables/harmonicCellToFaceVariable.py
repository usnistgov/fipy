#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "harmonicCellToFaceVariable.py"
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
from fipy.tools import numerix
from fipy.tools import inline

class _HarmonicCellToFaceVariable(_CellToFaceVariable):
    if inline.doInline:
        def _calcValue_(self, alpha, id1, id2):
            val = self._array.copy()

            inline._runIterateElementInline("""
                int ID1 = ITEM(id1, i, NULL);
                int ID2 = ITEM(id2, i, NULL);
                double cell1 = ITEM(var, ID1, vec);
                double cell2 = ITEM(var, ID2, vec);
                double cell1Xcell2 = cell1 * cell2;
                double tmp = ((cell2 - cell1) * ITEM(alpha, i, NULL) + cell1);
                if (tmp != 0 && cell1Xcell2 > 0.) {
                    ITEM(val, i, vec) = cell1Xcell2 / tmp;
                } else {
                    ITEM(val, i, vec) = 0.;
                }
            """,
            var = self.var.numericValue,
            val = val,
            alpha = alpha,
            id1 = id1, id2 = id2,
            shape=numerix.array(numerix.shape(val)),
            ni = self.mesh.numberOfFaces)

            return self._makeValue(value = val)
    else:
        def _calcValue_(self, alpha, id1, id2):
            cell1 = numerix.take(self.var,id1, axis=-1)
            cell2 = numerix.take(self.var,id2, axis=-1)
            value = ((cell2 - cell1) * alpha + cell1)
            eps = 1e-20
            value = (value == 0.) * eps + (value != 0.) * value
            cell1Xcell2 = cell1 * cell2
            value = ((value > eps) | (value < -eps)) * cell1Xcell2 / value
            value = (cell1Xcell2 >= 0.) * value

            return value
