from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

__all__ = []

class _AbstractTopology(object):
    _concatenatedClass = None

    def __init__(self, mesh):
        self.mesh = mesh

    @property
    def _isOrthogonal(self):
        raise NotImplementedError

    @property
    def _globalNonOverlappingCellIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Does not include the IDs of boundary cells.

        E.g., would return [0, 1, 4, 5] for mesh A

        ```
            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfCells)

    @property
    def _globalOverlappingCellIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Includes the IDs of boundary cells.

        E.g., would return [0, 1, 2, 4, 5, 6] for mesh A

        ```
            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfCells)

    @property
    def _localNonOverlappingCellIDs(self):
        """Return the IDs of the local mesh in isolation.

        Does not include the IDs of boundary cells.

        E.g., would return [0, 1, 2, 3] for mesh A

        ```
            A        B
        ------------------
        | 3 | 4 || 4 | 5 |
        ------------------
        | 0 | 1 || 1 | 2 |
        ------------------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfCells)

    @property
    def _localOverlappingCellIDs(self):
        """Return the IDs of the local mesh in isolation.

        Includes the IDs of boundary cells.

        E.g., would return [0, 1, 2, 3, 4, 5] for mesh A

        ```
            A        B
        ------------------
        | 3 | 4 || 5 |   |
        ------------------
        | 0 | 1 || 2 |   |
        ------------------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfCells)

    @property
    def _globalNonOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Does not include the IDs of faces of boundary cells.

        E.g., would return [0, 1, 4, 5, 8, 9, 12, 13, 14, 17, 18, 19]
        for mesh A

        ```
            A   ||   B
        --8---9---10--11--
       17   18  19  20   21
        --4---5----6---7--
       12   13  14  15   16
        --0---1----2---3--
                ||
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfFaces)

    @property
    def _globalOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Includes the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13,
        14, 15, 17, 18, 19, 20] for mesh A

        ```
            A   ||   B
        --8---9---10--11--
       17   18  19  20   21
        --4---5----6---7--
       12   13  14  15   16
        --0---1----2---3--
                ||
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfFaces)

    @property
    def _localNonOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in isolation.

        Does not include the IDs of faces of boundary cells.

        E.g., would return [0, 1, 3, 4, 6, 7, 9, 10, 11, 13, 14, 15]
        for mesh A

        ```
            A   ||   B
        --6---7-----7---8--
       13   14 15/14 15   16
        --3---4-----4---5--
        9   10 11/10 11   12
        --0---1-----1---2--
                ||
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfFaces)

    @property
    def _localOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in isolation.

        Includes the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        12, 13, 14, 15, 16] for mesh A

        ```
            A   ||   B
        --6---7----8------
       13   14  15  16   |
        --3---4----5------
        9   10  11  12   |
        --0---1----2------
                ||
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfFaces)

    # abstract element types mutually understood by FiPy and other meshing systems
    # (VTK, Gmsh, etc.)
    _elementTopology = dict([(k, v) for (v, k) in enumerate(("vertex",
                                                             "line",
                                                             "triangle",
                                                             "quadrangle",
                                                             "pixel",
                                                             "polygon",
                                                             "tetrahedron",
                                                             "hexahedron",
                                                             "voxel",
                                                             "prism",
                                                             "pyramid",
                                                             "unknown"))])

    @property
    def _cellTopology(self):
        """return a map of the topology of each cell of grid"""
        raise NotImplementedError
