from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable
from fipy.variables.distanceVariable import DistanceVariable

class GapFillDistanceVariable(DistanceVariable):

    def extendVariable(self, extensionVariable, order=2):
        if not hasattr(self, 'fineDistanceVariable'):
            self.fineDistanceVariable = DistanceVariable(mesh=self.mesh.fineMesh)
        if not hasattr(self, 'fineExtensionVariable'):
            self.fineExtensionVariable = CellVariable(mesh=self.mesh.fineMesh)
        self.fineDistanceVariable[:] = self(self.mesh.fineMesh.cellCenters)
        self.fineExtensionVariable[:] = extensionVariable(self.mesh.fineMesh.cellCenters)
        self.fineDistanceVariable.extendVariable(self.fineExtensionVariable, order=order)
        extensionVariable[:] = self.fineExtensionVariable(self.mesh.cellCenters)

    def calcDistanceFunction(self, order=2):
        if not hasattr(self, 'fineDistanceVariable'):
            self.fineDistanceVariable = DistanceVariable(mesh=self.mesh.fineMesh)
        self.fineDistanceVariable[:] = self(self.mesh.fineMesh.cellCenters)
        self.fineDistanceVariable.calcDistanceFunction(order=order)
        self[:] = self.fineDistanceVariable(self.mesh.cellCenters)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
