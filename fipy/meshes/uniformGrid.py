from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.meshes.abstractMesh import AbstractMesh

__all__ = ["UniformGrid"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class UniformGrid(AbstractMesh):
    """Wrapped scaled geometry properties"""
    @property
    def _scaledFaceAreas(self):
        return self._faceAreas

    @property
    def _scaledCellVolumes(self):
        return self._cellVolumes

    @property
    def _scaledCellCenters(self):
        return self._cellCenters

    @property
    def _scaledCellDistances(self):
        return self._cellDistances

    @property
    def _scaledCellToCellDistances(self):
        return self._cellToCellDistances

    @property
    def _scaledFaceToCellDistances(self):
        return self._faceToCellDistances

    """Geometry properties common to 1D, 2D, 3D"""
    @property
    def _orientedFaceNormals(self):
        return self.faceNormals

    @property
    def _faceCellToCellNormals(self):
        return self.faceNormals

    def _getFaceToCellDistances(self):
        return self._internalFaceToCellDistances

    def _setFaceToCellDistances(self, v):
        self._internalFaceToCellDistances = v
        self._setScaledValues()

    _faceToCellDistances = property(_getFaceToCellDistances,
                                    _setFaceToCellDistances)
