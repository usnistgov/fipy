#!/usr/bin/env python


__docformat__ = 'restructuredtext'

from fipy.terms.abstractDiffusionTerm import _AbstractDiffusionTerm
from fipy.tools import numerix

__all__ = ["DiffusionTermCorrection"]

class DiffusionTermCorrection(_AbstractDiffusionTerm):

    def _getNormals(self, mesh):
        return mesh._faceCellToCellNormals

    def _treatMeshAsOrthogonal(self, mesh):
        return mesh._isOrthogonal
