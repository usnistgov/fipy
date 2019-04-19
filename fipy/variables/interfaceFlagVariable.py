__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable
from fipy.tools.numerix import MA
from fipy.tools import numerix

class _InterfaceFlagVariable(CellVariable):
    def __init__(self, distanceVar):
        """
        Creates an `_InterfaceFlagVariable` object.

        :Parameters:
          - `distanceVar` : A `DistanceVariable` object.

        """
        CellVariable.__init__(self, distanceVar.mesh, hasOld=False)
        self.distanceVar = self._requires(distanceVar)

    def _calcValue(self):
        flag = MA.filled(numerix.take(self.distanceVar._interfaceFlag, self.mesh.cellFaceIDs), 0)
        flag = numerix.sum(flag, axis=0)
        return numerix.where(numerix.logical_and(self.distanceVar.value > 0, flag > 0), 1, 0)
