from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix

from fipy.meshes.mesh1D import Mesh1D
from fipy.meshes.mesh2D import Mesh2D
from fipy.meshes.mesh import Mesh

from fipy.meshes.topologies.abstractTopology import _AbstractTopology

class _GridTopology(_AbstractTopology):

    @property
    def _isOrthogonal(self):
        return True

class _Grid1DTopology(_GridTopology):

    _concatenatedClass = Mesh1D

    @property
    def _globalNonOverlappingCellIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Does not include the IDs of boundary cells.

        E.g., would return [0, 1] for mesh A

        ```
            A        B
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        ```

        .. note:: Trivial except for parallel meshes
        """

        return numerix.arange(self.mesh.offset + self.mesh.overlap['left'],
                              self.mesh.offset + self.mesh.nx - self.mesh.overlap['right'])

    @property
    def _globalOverlappingCellIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Includes the IDs of boundary cells.

        E.g., would return [0, 1, 2] for mesh A

        ```
            A        B
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.offset, self.mesh.offset + self.mesh.nx)

    @property
    def _localNonOverlappingCellIDs(self):
        """Return the IDs of the local mesh in isolation.

        Does not include the IDs of boundary cells.

        E.g., would return [0, 1] for mesh A

        ```
            A        B
        ------------------
        | 0 | 1 || 1 | 2 |
        ------------------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.overlap['left'],
                              self.mesh.nx - self.mesh.overlap['right'])

    @property
    def _localOverlappingCellIDs(self):
        """Return the IDs of the local mesh in isolation.

        Includes the IDs of boundary cells.

        E.g., would return [0, 1, 2] for mesh A

        ```
            A        B
        ------------------
        | 0 | 1 || 2 |   |
        ------------------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(0, self.mesh.nx)

    @property
    def _globalNonOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Does not include the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2] for mesh A

        ```
            A    ||   B
        ------------------
        0   1    2   3   4
        ------------------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.offset + self.mesh.overlap['left'],
                              self.mesh.offset + self.mesh.numberOfFaces - self.mesh.overlap['right'])

    @property
    def _globalOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Includes the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2, 3] for mesh A

        ```
            A    ||   B
        ------------------
        0   1    2   3   4
        ------------------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.offset, self.mesh.offset + self.mesh.numberOfFaces)

    @property
    def _localNonOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in isolation.

        Does not include the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2] for mesh A

        ```
            A    ||   B
        ------------------
        0   1   2/1  2   3
        ------------------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.overlap['left'],
                              self.mesh.numberOfFaces - self.mesh.overlap['right'])

    @property
    def _localOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in isolation.

        Includes the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2, 3] for mesh A

        ```
            A   ||   B
        ------------------
        0   1   2   3    |
        ------------------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(0, self.mesh.numberOfFaces)

    @property
    def _cellTopology(self):
        """return a map of the topology of each cell of grid"""
        cellTopology = numerix.empty((self.mesh.numberOfCells,), dtype=numerix.ubyte)
        cellTopology[:] = self._elementTopology["line"]

        return cellTopology

class _Grid2DTopology(_GridTopology):

    _concatenatedClass = Mesh2D

    @property
    def _globalNonOverlappingCellIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Does not include the IDs of boundary cells.

        E.g., would return [0, 1] for mesh A

        ```
        ---------
        | 4 | 5 |
        ---------  B
        | 2 | 3 |
        =========
        | 0 | 1 |  A
        ---------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange((self.mesh.offset[1] + self.mesh.overlap['bottom']) * self.mesh.nx,
                              (self.mesh.offset[1] + self.mesh.ny - self.mesh.overlap['top']) * self.mesh.nx)

    @property
    def _globalOverlappingCellIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Includes the IDs of boundary cells.

        E.g., would return [0, 1, 2, 3] for mesh A

        ```
        ---------
        | 4 | 5 |
        ---------  B
        | 2 | 3 |
        =========
        | 0 | 1 |  A
        ---------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.offset[1] * self.mesh.nx, (self.mesh.offset[1] + self.mesh.ny) * self.mesh.nx)

    @property
    def _localNonOverlappingCellIDs(self):
        """Return the IDs of the local mesh in isolation.

        Does not include the IDs of boundary cells.

        E.g., would return [0, 1] for mesh A

        ```
        ---------
        | 4 | 5 |
        ---------  B
        | 2 | 3 |
        =========
        | 0 | 1 |  A
        ---------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.overlap['bottom'] * self.mesh.nx,
                              (self.mesh.ny - self.mesh.overlap['top']) * self.mesh.nx)

    @property
    def _localOverlappingCellIDs(self):
        """Return the IDs of the local mesh in isolation.

        Includes the IDs of boundary cells.

        E.g., would return [0, 1, 2, 3] for mesh A

        ```
        ---------
        |   |   |
        ---------  B
        | 2 | 3 |
        =========
        | 0 | 1 |  A
        ---------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(0, self.mesh.ny * self.mesh.nx)

    def _calcFaceIDs(self, y0, ny, global_horz_faces):
        prev_horz_faces = y0 * self.mesh.nx
        horz = numerix.arange(prev_horz_faces,
                              prev_horz_faces + (ny + 1) * self.mesh.nx)
        prev_vert_faces = y0 * (self.mesh.nx + 1)
        vert = numerix.arange(global_horz_faces + prev_vert_faces,
                              global_horz_faces + prev_vert_faces
                              + ny * (self.mesh.nx + 1))

        return numerix.concatenate((horz, vert))

    @property
    def _globalNonOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Does not include the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2, 3, 8, 9, 10] for mesh A

        ```
        --6---7--
        14  15  16
        --4---5--  B
        11  12  13
        ==2===3==
        8   9   10 A
        --0---1--
        ```

        .. note:: Trivial except for parallel meshes
        """
        return self._calcFaceIDs(y0=(self.mesh.offset[1]
                                     + self.mesh.overlap['bottom']),
                                 ny=(self.mesh.ny
                                     - self.mesh.overlap['bottom']
                                     - self.mesh.overlap['top']),
                                 global_horz_faces=(self.mesh.nx
                                                    * (self.mesh.args['ny'] + 1)))

    @property
    def _globalOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Includes the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13] for mesh A

        ```
        --6---7--
        14  15  16
        --4---5--  B
        11  12  13
        ==2===3==
        8   9   10 A
        --0---1--
        ```

        .. note:: Trivial except for parallel meshes
        """
        return self._calcFaceIDs(y0=self.mesh.offset[1],
                                 ny=self.mesh.ny,
                                 global_horz_faces=(self.mesh.nx
                                                    * (self.mesh.args['ny'] + 1)))

    @property
    def _localNonOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in isolation.

        Does not include the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2, 3, 6, 7, 8] for mesh A
        and [2, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16] for mesh B

        ```
        --6---7--
        14  15  16
        --4---5--  B
        11  12  13
        ==2===3==
        6   7   8  A
        --0---1--
        ```

        .. note:: Trivial except for parallel meshes
        """
        return self._calcFaceIDs(y0=self.mesh.overlap['bottom'],
                                 ny=(self.mesh.ny
                                     - self.mesh.overlap['bottom']
                                     - self.mesh.overlap['top']),
                                 global_horz_faces=0)

    @property
    def _localOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in isolation.

        Includes the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] for mesh A
        and [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16] for mesh B

        ```
        --6---7--
        14  15  16
        --4---5--  B
        11  12  13
        ==2===3==
        6   7   8  A
        --0---1--
        ```

        .. note:: Trivial except for parallel meshes
        """
        return self._calcFaceIDs(y0=0, ny=self.mesh.ny, global_horz_faces=0)

    @property
    def _cellTopology(self):
        """return a map of the topology of each cell of grid"""
        cellTopology = numerix.empty((self.mesh.numberOfCells,), dtype=numerix.ubyte)
        cellTopology[:] = self._elementTopology["pixel"]

        return cellTopology


class _Grid3DTopology(_GridTopology):

    _concatenatedClass = Mesh

    @property
    def _globalNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Does not include the IDs of boundary cells.

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange((self.mesh.offset[2] + self.mesh.overlap['front']) * self.mesh.nx * self.mesh.ny,
                              (self.mesh.offset[2] + self.mesh.nz - self.mesh.overlap['back']) * self.mesh.nx * self.mesh.ny)

    @property
    def _globalOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Includes the IDs of boundary cells.

        .. note:: Trivial except for parallel meshes
        """

        return numerix.arange(self.mesh.offset[2] * self.mesh.nx * self.mesh.ny,
                              (self.mesh.offset[2] + self.mesh.nz) * self.mesh.nx * self.mesh.ny)

    @property
    def _localNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation.
        Does not include the IDs of boundary cells.

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.overlap['front'] * self.mesh.nx * self.mesh.ny,
                              (self.mesh.nz - self.mesh.overlap['back']) * self.mesh.nx * self.mesh.ny)

    @property
    def _localOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation.
        Includes the IDs of boundary cells.

        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(0, self.mesh.ny * self.mesh.nx * self.mesh.nz)

    def _calcFaceIDs(self, z0, nz, global_xy_faces, global_xz_faces):
        prev_xy_faces = z0 * self.mesh.nx * self.mesh.ny
        xy = numerix.arange(prev_xy_faces,
                            prev_xy_faces +
                            (nz + 1) * self.mesh.nx * self.mesh.ny)

        prev_xz_faces = z0 * self.mesh.nx * (self.mesh.ny + 1)
        xz = numerix.arange(global_xy_faces + prev_xz_faces,
                            global_xy_faces + prev_xz_faces +
                            self.mesh.nx * (self.mesh.ny + 1) * nz)

        prev_yz_faces = z0 * (self.mesh.nx + 1) * self.mesh.ny
        yz = numerix.arange(global_xy_faces + global_xz_faces + prev_yz_faces,
                            global_xy_faces + global_xz_faces + prev_yz_faces +
                            (self.mesh.nx + 1) * self.mesh.ny * nz)

        return numerix.concatenate((xy, xz, yz))

    @property
    def _globalNonOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Does not include the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15
                            20, 21, 22, 23, 24, 25] for mesh A
        and [4, 5, 6, 7, 16, 17, 18, 19, 26, 27, 28] for mesh B

        ```
              ---------
             /   /   / |
            ---------| /
           /   /   / |/
          =========: / A
         /   /   / :/
        ---------  /         y
        |   |   | /  B       |
        ---------            --x
                            /
                           z
        XY    ---------
              | 0 | 1 |
            --------- -
            | 2 | 3 |
          ========= -  A
          | 4 | 5 |
        --------- =
        | 6 | 7 |    B
        ---------

        XZ    ---------
             / 10/ 11/
            --------- -
           / 14/ 15/ /
          ========= -  A
         / 18/ 19/ /
        --------- =
         / 16/ 17/   B
        ---------

        YZ   /|  /|  /|
            /2| /2| /2|
           /|0//|1//|2/
          /2|//2|//2|/ A
         /:3//:4//:5/
        /2://2://2:/
        |6/ |7/ |8/  B
        |/  |/  |/

              ---------
             / 10/ 11/2|
            ---------|2/
           / 14/ 15/2|/
          =========|5/ A
         / 18/ 19/:|/
        ---------2:/
        | 6 | 7 |8/  B
        ---------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return self._calcFaceIDs(z0=(self.mesh.offset[2]
                                     + self.mesh.overlap['front']),
                                 nz=(self.mesh.nz
                                     - self.mesh.overlap['front']
                                     - self.mesh.overlap['back']),
                                 global_xy_faces=(self.mesh.nx
                                                  * self.mesh.ny
                                                  * (self.mesh.args['nz'] + 1)),
                                 global_xz_faces=(self.mesh.nx
                                                  * (self.mesh.ny + 1)
                                                  * self.mesh.args['nz']))

    @property
    def _globalOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh.

        Includes the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2, 3, 4, 5, 6, 7,
                            8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                            20, 21, 22, 23, 24, 25, 26, 27, 28] for mesh A
        and [2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 16, 17, 18, 19,
             23, 24, 25, 26, 27, 28] for mesh B

        ```
              ---------
             /   /   / |
            ---------| /
           /   /   / |/
          =========: / A
         /   /   / :/
        ---------  /         y
        |   |   | /  B       |
        ---------            --x
                            /
                           z
        XY    ---------
              | 0 | 1 |
            --------- -
            | 2 | 3 |
          ========= -  A
          | 4 | 5 |
        --------- =
        | 6 | 7 |    B
        ---------

        XZ    ---------
             / 10/ 11/
            --------- -
           / 14/ 15/ /
          ========= -  A
         / 18/ 19/ /
        --------- =
         / 16/ 17/   B
        ---------

        YZ   /|  /|  /|
            /2| /2| /2|
           /|0//|1//|2/
          /2|//2|//2|/ A
         /:3//:4//:5/
        /2://2://2:/
        |6/ |7/ |8/  B
        |/  |/  |/

              ---------
             / 10/ 11/2|
            ---------|2/
           / 14/ 15/2|/
          =========|5/ A
         / 18/ 19/:|/
        ---------2:/
        | 6 | 7 |8/  B
        ---------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return self._calcFaceIDs(z0=self.mesh.offset[2],
                                 nz=self.mesh.nz,
                                 global_xy_faces=(self.mesh.nx
                                                  * self.mesh.ny
                                                  * (self.mesh.args['nz'] + 1)),
                                 global_xz_faces=(self.mesh.nx
                                                  * (self.mesh.ny + 1)
                                                  * self.mesh.args['nz']))

    @property
    def _localNonOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in isolation.

        Does not include the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15
                            20, 21, 22, 23, 24, 25] for mesh A
        and [2, 3, 4, 5, 10, 11, 12, 13, 17, 18, 19] for mesh B

        ```
              ---------
             /   /   / |
            ---------| /
           /   /   / |/
          =========: / A
         /   /   / :/
        ---------  /         y
        |   |   | /  B       |
        ---------            --x
                            /
                           z
        XY    ---------
              | 0 | 1 |
            --------- -
            | 2 | 3 |
          ========= -  A
          |4/2|5/3|
        --------- =
        | 4 | 5 |    B
        ---------

        XZ    ---------
             / 10/ 11/
            --------- -
           / 14/ 15/ /
          ========= -  A
         / 12/ 13/ /
        --------- =
         / 10/ 11/   B
        ---------

        YZ   /|  /|  /|
            /2| /2| /2|
           /|0//|1//|2/
          /2|//2|//2|/ A
         /:3//:4//:5/
        /1://1://1:/
        |7/ |8/ |9/  B
        |/  |/  |/

              ---------
             / 10/ 11/2|
            ---------|2/
           / 14/ 15/2|/
          =========|5/ A
         / 12/ 13/:|/
        ---------1:/
        | 4 | 5 |9/  B
        ---------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return self._calcFaceIDs(z0=self.mesh.overlap['front'],
                                 nz=(self.mesh.nz
                                     - self.mesh.overlap['front']
                                     - self.mesh.overlap['back']),
                                 global_xy_faces=0,
                                 global_xz_faces=0)

    @property
    def _localOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in isolation.

        Includes the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2, 3, 4, 5, 6, 7,
                            8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                            20, 21, 22, 23, 24, 25, 26, 27, 28] for mesh A
        and [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19] for mesh B

        ```
              ---------
             /   /   / |
            ---------| /
           /   /   / |/
          =========: / A
         /   /   / :/
        ---------  /         y
        |   |   | /  B       |
        ---------            --x
                            /
                           z
        XY    ---------
              | 0 | 1 |
            --------- -
            | 2 | 3 |
          ========= -  A
          |4/2|5/3|
        --------- =
        | 4 | 5 |    B
        ---------

        XZ    ---------
             / 10/ 11/
            --------- -
           / 14/ 15/ /
          ========= -  A
         / 12/ 13/ /
        --------- =
         / 10/ 11/   B
        ---------

        YZ   /|  /|  /|
            /2| /2| /2|
           /|0//|1//|2/
          /2|//2|//2|/ A
         /:3//:4//:5/
        /1://1://1:/
        |7/ |8/ |9/  B
        |/  |/  |/

              ---------
             / 10/ 11/2|
            ---------|2/
           / 14/ 15/2|/
          =========|5/ A
         / 12/ 13/:|/
        ---------1:/
        | 4 | 5 |9/  B
        ---------
        ```

        .. note:: Trivial except for parallel meshes
        """
        return self._calcFaceIDs(z0=0,
                                 nz=self.mesh.nz,
                                 global_xy_faces=0,
                                 global_xz_faces=0)

    @property
    def _cellTopology(self):
        """return a map of the topology of each cell of grid"""
        cellTopology = numerix.empty((self.mesh.numberOfCells,), dtype=numerix.ubyte)
        cellTopology[:] = self._elementTopology["voxel"]

        return cellTopology

class _PeriodicGrid1DTopology(_Grid1DTopology):

    @property
    def _globalOverlappingCellIDs(self):
        return super(_PeriodicGrid1DTopology, self)._globalOverlappingCellIDs % self.mesh.args['nx']
