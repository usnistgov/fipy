from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.terms.abstractDiffusionTerm import _AbstractDiffusionTerm

__all__ = ["DiffusionTermNoCorrection"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class DiffusionTermNoCorrection(_AbstractDiffusionTerm):
    def _getNormals(self, mesh):
        return mesh.faceNormals

    def _treatMeshAsOrthogonal(self, mesh):
        return True
