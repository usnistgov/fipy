from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext'

from ...tools import numerix
from ...variables.cellVariable import CellVariable

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
    def _cellProcID(self):
        """
        Return the processor that "owns" each cell.
        Ownership goes to the lowest `procID` of a neighboring cell.
        Boundary faces try not to be duplicated

        E.g., would return [0, 0, 1, 0, 0, 1] for mesh A
                       and [0, 1, 1, 0, 1, 1] for mesh B

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
        procID = CellVariable(mesh=self.mesh,
                              value=self.mesh.communicator.procID)
        procID._updateGhosts()

        return procID

    @property
    def _ownedFaceIDs(self):
        """
        Return the local face IDs of the mesh "owned" by the current
        processor.  Ownership goes to the lowest `procID` of a neighboring
        cell.  Boundary faces try not to be duplicated

        E.g., would return [0, 1, 3, 4, 6, 7, 9, 10, 11, 13, 14, 15] for mesh A
                       and [1, 2, 4, 5, 7, 8, 11, 12, 15, 16] for mesh B

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
        minproc = numerix.take(self._cellProcID,
                               self.mesh.faceCellIDs).min(axis=0)

        return numerix.where(minproc == self.mesh.communicator.procID)[0]

    @property
    def _nonOverlappingFaces(self):
        """Return mask of faces on local mesh.

        False for faces only belonging to ghost cells.

        E.g., would return [True, True, False, True, True, False, True,
        True, False, True, True, True, False, True, True, True, False]
        for mesh A

        ```
            A   ||   B
        --6---7----8------
       13   14  15  16   |
        --3---4----5------
        9   10  11  12   |
        --0---1----2------
                ||
        ```
        """
        return (numerix.isin(self.mesh.faceCellIDs[0],
                             self._localNonOverlappingCellIDs).filled(False)
                | numerix.isin(self.mesh.faceCellIDs[1],
                               self._localNonOverlappingCellIDs).filled(False))

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

    @property
    def _vertexCellIDs(self):
        """Return cell IDs bounded by each vertex

        E.g., would return

        [[ 0  0  1  0  0  1  2  2  3]
         [--  1 --  2  1  3 --  3 --]
         [-- -- -- --  2 -- -- -- --]
         [-- -- -- --  3 -- -- -- --]]

        given

        ```
        6-------7-------8
        |       |       |
        |   2   |   3   |
        |       |       |
        3-------4-------5
        |       |       |
        |   0   |   1   |
        |       |       |
        0-------1-------2
        ```
        """
        return numerix._invert_indices(self.mesh._cellVertexIDs)

    @property
    def _vertexFaceIDs(self):
        """Return face IDs bounded by each vertex

        E.g., would return

        [[  0   0   1   2   2   3   4   4   5]
         [  6   1   8   6   3   8   9   5  11]
         [ --   7  --   9   7  11  --  10  --]
         [ --  --  --  --  10  --  --  --  --]]

        given

        ```
        6---4---7---5---8
        |       |       |
        9       10      11
        |       |       |
        3---2---4---3---5
        |       |       |
        6       7       8
        |       |       |
        0---0---1---1---2
        ```
        """
        return numerix._invert_indices(self.mesh.faceVertexIDs)

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
