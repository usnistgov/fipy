from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.meshes.builders.grid1DBuilder import _NonuniformGrid1DBuilder

class _PeriodicGrid1DBuilder(_NonuniformGrid1DBuilder):

    def buildGridData(self, *args, **kwargs):
        kwargs["cacheOccupiedNodes"] = True
        return super(_PeriodicGrid1DBuilder, self).buildGridData(*args,
                                                                **kwargs)

    def _buildOverlap(self, overlap, procID, occupiedNodes):
        if occupiedNodes == 1:
            return super(_PeriodicGrid1DBuilder, self)._buildOverlap(overlap,
                     procID, occupiedNodes)
        else:
            return (overlap, overlap, {'left': overlap, 'right': overlap})
