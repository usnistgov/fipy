"""Generic mesh class using numerix to do the calculations

    Meshes contain cells, faces, and vertices.

    This is built for a non-mixed element mesh.
"""
from __future__ import unicode_literals

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools import serialComm

from fipy.meshes.mesh import Mesh
from fipy.meshes.representations.meshRepresentation import _MeshRepresentation
from fipy.meshes.topologies.meshTopology import _Mesh1DTopology

__all__ = ["Mesh1D"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Mesh1D(Mesh):
    def __init__(self, vertexCoords, faceVertexIDs, cellFaceIDs, communicator=serialComm, _RepresentationClass=_MeshRepresentation, _TopologyClass=_Mesh1DTopology):
        super(Mesh1D, self).__init__(vertexCoords=vertexCoords, faceVertexIDs=faceVertexIDs, cellFaceIDs=cellFaceIDs, communicator=communicator,
                                     _RepresentationClass=_RepresentationClass, _TopologyClass=_TopologyClass)

    def _calcScaleArea(self):
        return 1.

    def _calcScaleVolume(self):
        return self.scale['length']

    def _calcFaceAreas(self):
        return numerix.ones(self.numberOfFaces, 'd')

    def _calcFaceNormals(self):
        faceNormals = numerix.array((numerix.ones(self.numberOfFaces, 'd'),))
        # The left-most face has neighboring cells None and the left-most cell.
        # We must reverse the normal to make fluxes work correctly.
        if self.numberOfFaces > 0:
            faceNormals[..., 0] = -faceNormals[..., 0]
        return faceNormals

    def _calcFaceTangents(self):
        faceTangents1 = numerix.zeros(self.numberOfFaces, 'd')[numerix.NewAxis, ...]
        faceTangents2 = numerix.zeros(self.numberOfFaces, 'd')[numerix.NewAxis, ...]
        return faceTangents1, faceTangents2

    def _translate(self, vector):
        newCoords = self.vertexCoords + vector
        newmesh = Mesh1D(newCoords, numerix.array(self.faceVertexIDs), numerix.array(self.cellFaceIDs))
        return newmesh

    def __mul__(self, factor):
        newCoords = self.vertexCoords * factor
        newmesh = Mesh1D(newCoords, numerix.array(self.faceVertexIDs), numerix.array(self.cellFaceIDs))
        return newmesh

    @property
    def _VTKCellType(self):
        try:
            from tvtk.api import tvtk
        except ImportError as e:
            from enthought.tvtk.api import tvtk
        return tvtk.Line().cell_type
