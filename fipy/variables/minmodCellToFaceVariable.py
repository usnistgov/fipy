from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.variables.cellToFaceVariable import _CellToFaceVariable
from fipy.tools import numerix

class _MinmodCellToFaceVariable(_CellToFaceVariable):
    def _calcValue_(self, alpha, id1, id2):
        cell1 = numerix.take(self.var, id1, axis=-1)
        cell2 = numerix.take(self.var, id2, axis=-1)
        return numerix.where((cell1 > 0) & (cell2 > 0),
                             numerix.minimum(cell1, cell2),
                             numerix.where((cell1 < 0) & (cell2 < 0),
                                           numerix.maximum(cell1, cell2),
                                           0))
