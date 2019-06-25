"""
Periodic 1D Mesh
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.meshes.nonUniformGrid1D import NonUniformGrid1D
from fipy.meshes.builders import _PeriodicGrid1DBuilder
from fipy.meshes.topologies.gridTopology import _PeriodicGrid1DTopology

__all__ = ["PeriodicGrid1D"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class PeriodicGrid1D(NonUniformGrid1D):
    """

    Creates a Periodic grid mesh.

        >>> mesh = PeriodicGrid1D(dx = (1, 2, 3))

        >>> print(numerix.allclose(numerix.nonzero(mesh.exteriorFaces)[0],
        ...                        [3])) # doctest: +PROCESSOR_0
        True

        >>> print(numerix.allclose(mesh.faceCellIDs.filled(-999),
        ...                        [[2, 0, 1, 2],
        ...                         [0, 1, 2, -999]])) # doctest: +PROCESSOR_0
        True

        >>> print(numerix.allclose(mesh._cellDistances,
        ...                        [ 2., 1.5, 2.5, 1.5])) # doctest: +PROCESSOR_0
        True

        >>> print(numerix.allclose(mesh._cellToCellDistances,
        ...                        [[ 2.,   1.5,  2.5],
        ...                         [ 1.5,  2.5,  2. ]])) # doctest: +PROCESSOR_0
        True

        >>> print(numerix.allclose(mesh.faceNormals,
        ...                        [[ 1.,  1.,  1.,  1.]])) # doctest: +PROCESSOR_0
        True

        >>> print(numerix.allclose(mesh._cellVertexIDs,
        ...                        [[1, 2, 2],
        ...                        [0, 1, 0]])) # doctest: +PROCESSOR_0
        True
    """
    def __init__(self, dx = 1., nx = None, overlap=2, *args, **kwargs):

        super(PeriodicGrid1D, self).__init__(dx = dx, nx = nx, overlap=overlap,
                                             _BuilderClass=_PeriodicGrid1DBuilder,
                                             _TopologyClass=_PeriodicGrid1DTopology, *args, **kwargs)
        self._nonPeriodicCellFaceIDs = numerix.array(super(PeriodicGrid1D, self).cellFaceIDs)
        self._makePeriodic()

    def _makePeriodic(self):
        if self.occupiedNodes == 1:
            self._connectFaces(numerix.nonzero(self.facesLeft),
                               numerix.nonzero(self.facesRight))

    @property
    def cellCenters(self):
        """Defined outside of a geometry class since we need the `CellVariable`
        version of `cellCenters`; that is, the `cellCenters` defined in
        :mod:`fipy.meshes.mesh` and not in any geometry (since a `CellVariable` requires
        a reference to a mesh)."""
        return super(PeriodicGrid1D, self).cellCenters \
                % numerix.sum(self.globalNumberOfCells * self.args['dx'])

    def _translate(self, vector):
        """
        Test for ticket:298.

        >>> from fipy import *
        >>> m = PeriodicGrid1D(nx=2) + [[-1]]
        >>> print(CellVariable(mesh=m, value=m.cellCenters[0]))
        [-0.5  0.5]

        """
        newCoords = self.vertexCoords + vector
        newmesh = self.__class__(**self.args)
        from fipy.meshes.mesh1D import Mesh1D
        Mesh1D.__init__(newmesh, newCoords, numerix.array(self.faceVertexIDs), self._nonPeriodicCellFaceIDs, communicator=self.communicator)
        newmesh._makePeriodic()
        return newmesh

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()

