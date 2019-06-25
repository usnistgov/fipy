from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable
from fipy.tools.numerix import MA
from fipy.tools import numerix

__all__ = []

class _InterfaceAreaVariable(CellVariable):
    def __init__(self, distanceVar):
        """
        Creates an `_InterfaceAreaVariable` object.

        Parameters
        ----------
        distanceVar : ~fipy.variables.distanceVariable.DistanceVariable
        """
        CellVariable.__init__(self, distanceVar.mesh, hasOld=False)
        self.distanceVar = self._requires(distanceVar)

    def _calcValue(self):
        normals = numerix.array(MA.filled(self.distanceVar._cellInterfaceNormals, 0))
        areas = numerix.array(MA.filled(self.mesh._cellAreaProjections, 0))
        return numerix.sum(abs(numerix.dot(normals, areas)), axis=0)
