from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix

from fipy.meshes.topologies.abstractTopology import _AbstractTopology

class _MeshTopology(_AbstractTopology):

    @property
    def _concatenatedClass(self):
        from fipy.meshes.mesh import Mesh
        return Mesh

    @property
    def _isOrthogonal(self):
        return False

    @property
    def _cellTopology(self):
        """return a map of the topology of each cell"""
        facesPerCell = self.mesh._facesPerCell
        nodesPerFace = self.mesh._nodesPerFace

        def faceCountsMatch(targetCounts):
            if len(targetCounts) > nodesPerFace.shape[0]:
                # pad nodesPerFace with zeros
                paddedNodesPerFace = numerix.zeros((len(targetCounts), nodesPerFace.shape[1]), dtype=numerix.INT_DTYPE)
                paddedNodesPerFace[:nodesPerFace.shape[0],:] = nodesPerFace

                paddedTargetCounts = numerix.array(targetCounts)[..., numerix.newaxis]
            else:
                # pad target face node count with zeros
                paddedTargetCounts = numerix.concatenate((targetCounts,
                                                          [0] * (self.mesh._maxFacesPerCell - len(targetCounts))))
                paddedTargetCounts = paddedTargetCounts[..., numerix.newaxis]

                paddedNodesPerFace = nodesPerFace

            return ((facesPerCell == len(targetCounts))
                    & (paddedNodesPerFace == paddedTargetCounts).all(axis=0))

        cellTopology = numerix.empty((self.mesh.numberOfCells,), dtype=numerix.ubyte)

        t = self._elementTopology

        if self.mesh.dim == 1:
            cellTopology[:] = t["line"]
        elif self.mesh.dim == 2:
            cellTopology[:] = t["polygon"]
            cellTopology[faceCountsMatch([2, 2, 2])] = t["triangle"]
            cellTopology[faceCountsMatch([2, 2, 2, 2])] = t["quadrangle"]
        else:
            cellTopology[:] = t["unknown"]
            cellTopology[faceCountsMatch([3, 3, 3, 3])] = t["tetrahedron"]
            cellTopology[faceCountsMatch([4, 4, 4, 4, 4, 4])] = t["hexahedron"]
            cellTopology[faceCountsMatch([4, 4, 4, 3, 3])] = t["prism"]
            cellTopology[faceCountsMatch([4, 3, 3, 3, 3])] = t["pyramid"]

        return cellTopology

class _Mesh1DTopology(_MeshTopology):

    @property
    def _concatenatedClass(self):
        from fipy.meshes.mesh1D import Mesh1D
        return Mesh1D

    @property
    def _isOrthogonal(self):
        return True

    @property
    def _cellTopology(self):
        """return a map of the topology of each cell of grid"""
        cellTopology = numerix.empty((self.mesh.numberOfCells,), dtype=numerix.ubyte)
        cellTopology[:] = self._elementTopology["line"]

        return cellTopology

class _Mesh2DTopology(_MeshTopology):

    @property
    def _concatenatedClass(self):
        from fipy.meshes.mesh2D import Mesh2D
        return Mesh2D

    @property
    def _cellTopology(self):
        """return a map of the topology of each cell"""
        cellTopology = numerix.empty((self.mesh.numberOfCells,), dtype=numerix.ubyte)

        t = self._elementTopology
        cellTopology[:] = t["polygon"]

        facesPerCell = self.mesh._facesPerCell
        cellTopology[facesPerCell == 3] = t["triangle"]
        cellTopology[facesPerCell == 4] = t["quadrangle"]

        return cellTopology
