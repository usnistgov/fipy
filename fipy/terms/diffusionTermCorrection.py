from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.terms.abstractDiffusionTerm import _AbstractDiffusionTerm
from fipy.tools import numerix

__all__ = ["DiffusionTermCorrection"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class DiffusionTermCorrection(_AbstractDiffusionTerm):

    def _getNormals(self, mesh):
        return mesh._faceCellToCellNormals

    def _treatMeshAsOrthogonal(self, mesh):
        return mesh._isOrthogonal
