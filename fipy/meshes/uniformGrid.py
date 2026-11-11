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

    """Distance vectors, calculated as in :class:`~fipy.meshes.mesh.Mesh`"""
    @property
    def _cellToFaceDistanceVectors(self):
        faceCenters = MA.repeat(self._faceCenters[..., numerix.NewAxis, :], 2, 1)
        return faceCenters - numerix.take(self._cellCenters, self.faceCellIDs, axis=1)

    @property
    def _cellDistanceVectors(self):
        """Vector from the first to the second cell of each face, or from
        the cell to the face on the boundary

        Uniform grids give the same vectors as the equivalent non-uniform
        grids, which are calculated by :class:`~fipy.meshes.mesh.Mesh`
        (#706).

            >>> from fipy import Grid1D, Grid2D, Grid3D
            >>> for uniform, nonuniform in (
            ...         (Grid1D(nx=3, dx=0.5),
            ...          Grid1D(dx=[0.5] * 3)),
            ...         (Grid2D(nx=3, ny=2, dx=0.5, dy=2.),
            ...          Grid2D(dx=[0.5] * 3, dy=[2.] * 2)),
            ...         (Grid3D(nx=3, ny=2, nz=2, dx=0.5, dy=2., dz=1.5),
            ...          Grid3D(dx=[0.5] * 3, dy=[2.] * 2, dz=[1.5] * 2))):
            ...     u = uniform.cellToFaceDistanceVectors
            ...     n = nonuniform.cellToFaceDistanceVectors
            ...     print(type(uniform).__name__,
            ...           numerix.allclose(uniform.cellDistanceVectors,
            ...                            nonuniform.cellDistanceVectors),
            ...           numerix.allequal(MA.getmaskarray(u), MA.getmaskarray(n)),
            ...           numerix.allclose(MA.filled(u, 0.), MA.filled(n, 0.))) # doctest: +SERIAL
            UniformGrid1D True True True
            UniformGrid2D True True True
            UniformGrid3D True True True
        """
        vectors = numerix.take(self._cellCenters, self.faceCellIDs, axis=1)
        vectors = vectors[..., 1, :] - vectors[..., 0, :]
        return MA.filled(MA.where(MA.getmaskarray(vectors),
                                  self._cellToFaceDistanceVectors[:, 0],
                                  vectors))
