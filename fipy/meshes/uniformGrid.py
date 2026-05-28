__docformat__ = 'restructuredtext'

from fipy.meshes.abstractMesh import AbstractMesh
from fipy.tools import numerix
from fipy.tools.numerix import MA

__all__ = ["UniformGrid"]

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
    def _cellToFaceDistanceVectors(self):
        """Vectors from each adjacent cell center to each face center.

        Shape ``(dim, 2, numberOfFaces)``: the leading axis indexes the
        spatial component, the second axis selects which of the two
        cells adjacent to the face is the source (the second row is
        masked along boundary faces, matching the non-uniform grids).

        Needed for upwinding and for Robin boundary conditions; see #706.
        """
        tmp = MA.repeat(self._faceCenters[..., numerix.NewAxis, :], 2, 1)
        tmp = tmp - numerix.take(self._cellCenters, self.faceCellIDs, axis=1)
        return tmp

    @property
    def _cellDistanceVectors(self):
        """Vector between adjacent cell centers across each face.

        Shape ``(dim, numberOfFaces)``: ``c1 - c0`` where ``c0, c1`` are
        the two cells adjacent to a face.  At boundary faces (one cell
        missing) the result falls back to the cell-to-face vector of
        the present cell, matching the non-uniform grid convention.
        See #706.
        """
        tmp = numerix.take(self._cellCenters, self.faceCellIDs, axis=1)
        tmp = tmp[..., 1, :] - tmp[..., 0, :]
        tmp = MA.filled(MA.where(MA.getmaskarray(tmp),
                                 self._cellToFaceDistanceVectors[:, 0],
                                 tmp))
        return tmp

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
