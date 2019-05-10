from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix
from fipy.tools import inline
from fipy.variables.arithmeticCellToFaceVariable import _ArithmeticCellToFaceVariable

class _ModCellToFaceVariable(_ArithmeticCellToFaceVariable):
    def __init__(self, var, modIn):
        _ArithmeticCellToFaceVariable.__init__(self, var)
        self.modIn = modIn

    if inline.doInline:
        def  _calcValue_(self, alpha, id1, id2):
            val = self._array.copy()

            inline._runInline(self.modIn + """
            int ID1 = id1[i];
            int ID2 = id2[i];
            double cell2 = var[ID2];
            val[i] = mod(cell2 - var[ID1]) * alpha[i] + var[ID1];
            """, var = self.var.numericValue,
                val = val,
                alpha = alpha,
                id1 = id1, id2 = id2,
                ni = self.mesh.numberOfFaces)

            return self._makeValue(value = val)
