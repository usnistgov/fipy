from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.meshes.representations.abstractRepresentation import _AbstractRepresentation

class _MeshRepresentation(_AbstractRepresentation):

    def getstate(self):
        """Collect the necessary information to ``pickle`` the `Mesh` to persistent storage.
        """
        return dict(vertexCoords=self.mesh.vertexCoords *  self.mesh.scale['length'],
                    faceVertexIDs=self.mesh.faceVertexIDs,
                    cellFaceIDs=self.mesh.cellFaceIDs,
                    _RepresentationClass=self.__class__)

    @staticmethod
    def setstate(mesh, state):
        """Populate a new `Mesh` from ``pickled`` persistent storage.
        """
        from fipy.meshes.mesh import Mesh
        Mesh.__init__(mesh, **state)

    def repr(self):
        return "%s()" % self.mesh.__class__.__name__
