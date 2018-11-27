


__docformat__ = 'restructuredtext'

from fipy.terms.abstractDiffusionTerm import _AbstractDiffusionTerm

__all__ = ["DiffusionTermNoCorrection"]

class DiffusionTermNoCorrection(_AbstractDiffusionTerm):
    def _getNormals(self, mesh):
        return mesh.faceNormals

    def _treatMeshAsOrthogonal(self, mesh):
        return True
