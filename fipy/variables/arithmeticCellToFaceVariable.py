from __future__ import unicode_literals
__all__ = []

from fipy.variables.cellToFaceVariable import _CellToFaceVariable
from fipy.tools import numerix
from fipy.tools import inline

class _ArithmeticCellToFaceVariable(_CellToFaceVariable):
    if inline.doInline:
        def _calcValue_(self, alpha, id1, id2):
            val = self._array.copy()

            inline._runIterateElementInline("""
                int ID1 = ITEM(id1, i, NULL);
                int ID2 = ITEM(id2, i, NULL);
                double cell1 = ITEM(var, ID1, vec);
                double cell2 = ITEM(var, ID2, vec);
                ITEM(val, i, vec) = (cell2 - cell1) * ITEM(alpha, i, NULL) + cell1;
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
            return (cell2 - cell1) * alpha + cell1
