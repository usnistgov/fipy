
from fipy.tools import numerix
from meshGeometry1D import MeshGeometry1D

class PeriodicGridGeometry1D(MeshGeometry1D):

    @property
    def cellCenters(self):
        return super(PeriodicGridGeometry1D, self).cellCenters \
                % numerix.sum(self.mesh.globalNumberOfCells * self.mesh.args['dx']) 

