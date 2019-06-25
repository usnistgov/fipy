"""
2D periodic rectangular Mesh
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools import parallelComm
from fipy.meshes.nonUniformGrid2D import NonUniformGrid2D

__all__ = ["PeriodicGrid2D", "PeriodicGrid2DLeftRight", "PeriodicGrid2DTopBottom"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class _BasePeriodicGrid2D(NonUniformGrid2D):
    def __init__(self, dx = 1., dy = 1., nx = None, ny = None, overlap=2, communicator=parallelComm, *args, **kwargs):
        super(_BasePeriodicGrid2D, self).__init__(dx = dx, dy = dy, nx = nx, ny = ny, overlap=overlap, communicator=communicator, *args, **kwargs)
        self._nonPeriodicCellVertexIDs = super(_BasePeriodicGrid2D, self)._cellVertexIDs
        self._orderedCellVertexIDs_data = super(_BasePeriodicGrid2D, self)._orderedCellVertexIDs
        self._nonPeriodicCellFaceIDs = numerix.array(super(_BasePeriodicGrid2D, self).cellFaceIDs)
        self._makePeriodic()

    @property
    def _cellVertexIDs(self):
        return self._nonPeriodicCellVertexIDs

    def _translate(self, vector):
        """
        Test for ticket:298.

        >>> from fipy import *
        >>> m = PeriodicGrid2DLeftRight(nx=2, ny=2) + [[-1], [0]]
        >>> orderedCellVertexIDs = [[1, 2, 4, 5],
        ...                         [4, 5, 7, 8],
        ...                         [3, 4, 6, 7],
        ...                         [0, 1, 3, 4]]
        >>> print(numerix.allclose(m._orderedCellVertexIDs, orderedCellVertexIDs))  # doctest: +PROCESSOR_0
        True
        >>> print(CellVariable(mesh=m, value=m.cellCenters[0]))
        [-0.5  0.5 -0.5  0.5]
        """
        newCoords = self.vertexCoords + vector
        newmesh = self.__class__(**self.args)
        from fipy.meshes.mesh2D import Mesh2D
        Mesh2D.__init__(newmesh, newCoords, self.faceVertexIDs, self._nonPeriodicCellFaceIDs, communicator=self.communicator)
        newmesh._makePeriodic()
        return newmesh

class PeriodicGrid2D(_BasePeriodicGrid2D):
    """
    Creates a periodic 2D grid mesh with horizontal faces numbered
    first and then vertical faces. Vertices and cells are numbered
    in the usual way.

        >>> from fipy import numerix

        >>> mesh = PeriodicGrid2D(dx = 1., dy = 0.5, nx = 2, ny = 2)

        >>> print(numerix.allclose(numerix.nonzero(mesh.exteriorFaces)[0],
        ...                        [ 4,  5,  8, 11]))  # doctest: +PROCESSOR_0
        True

        >>> print(numerix.allclose(mesh.faceCellIDs.filled(-1),
        ...                        [[2, 3, 0, 1, 2, 3, 1, 0, 1, 3, 2, 3],
        ...                         [0, 1, 2, 3, -1, -1, 0, 1, -1, 2, 3, -1]])) # doctest: +PROCESSOR_0
        True

        >>> print(numerix.allclose(mesh._cellDistances,
        ...                        [ 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 1., 1., 0.5, 1., 1., 0.5])) # doctest: +PROCESSOR_0
        True

        >>> print(numerix.allclose(mesh.cellFaceIDs,
        ...                        [[0, 1, 2, 3],
        ...                         [7, 6, 10, 9],
        ...                         [2, 3, 0, 1],
        ...                         [6, 7, 9, 10]])) # doctest: +PROCESSOR_0
        True

        >>> print(numerix.allclose(mesh._cellToCellDistances,
        ...                        [[ 0.5, 0.5, 0.5, 0.5],
        ...                         [ 1., 1., 1., 1. ],
        ...                         [ 0.5, 0.5, 0.5, 0.5],
        ...                         [ 1., 1., 1., 1. ]])) # doctest: +PROCESSOR_0
        True

        >>> normals = [[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
        ...            [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]]

        >>> print(numerix.allclose(mesh.faceNormals, normals)) # doctest: +PROCESSOR_0
        True

        >>> print(numerix.allclose(mesh._cellVertexIDs,
        ...                        [[4, 5, 7, 8],
        ...                         [3, 4, 6, 7],
        ...                         [1, 2, 4, 5],
        ...                         [0, 1, 3, 4]])) # doctest: +PROCESSOR_0
        True
    """

    def _makePeriodic(self):
        self._connectFaces(numerix.nonzero(self.facesLeft),
                           numerix.nonzero(self.facesRight))
        self._connectFaces(numerix.nonzero(self.facesBottom),
                           numerix.nonzero(self.facesTop))

class PeriodicGrid2DLeftRight(_BasePeriodicGrid2D):
    def _makePeriodic(self):
        self._connectFaces(numerix.nonzero(self.facesLeft),
                           numerix.nonzero(self.facesRight))

class PeriodicGrid2DTopBottom(_BasePeriodicGrid2D):
    def _makePeriodic(self):
        self._connectFaces(numerix.nonzero(self.facesBottom),
                           numerix.nonzero(self.facesTop))

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()

