from __future__ import division
from __future__ import unicode_literals
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
            cell1 = numerix.take(self.var, id1, axis=-1)
            cell2 = numerix.take(self.var, id2, axis=-1)
            value = ((cell2 - cell1) * alpha + cell1)
            eps = 1e-20
            value = (value == 0.) * eps + (value != 0.) * value
            cell1Xcell2 = cell1 * cell2
            value = ((value > eps) | (value < -eps)) * cell1Xcell2 / value
            value = (cell1Xcell2 >= 0.) * value

            return value
