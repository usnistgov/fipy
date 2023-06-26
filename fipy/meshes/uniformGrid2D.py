"""
2D rectangular Mesh with constant spacing in x and constant spacing in y
"""
from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools import inline
from fipy.tools import parallelComm

from fipy.meshes.uniformGrid import UniformGrid
from fipy.meshes.builders import _UniformGrid2DBuilder
from fipy.meshes.builders import _Grid2DBuilder
from fipy.meshes.representations.gridRepresentation import _Grid2DRepresentation
from fipy.meshes.topologies.gridTopology import _Grid2DTopology

__all__ = ["UniformGrid2D"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class UniformGrid2D(UniformGrid):
    """
    Creates a 2D grid mesh with horizontal faces numbered
    first and then vertical faces.
    """
    def __init__(self, dx=1., dy=1., nx=1, ny=1, origin=((0,), (0,)),
                       overlap=2, communicator=parallelComm,
                       _RepresentationClass=_Grid2DRepresentation,
                       _TopologyClass=_Grid2DTopology):


        super(UniformGrid2D, self).__init__(communicator=communicator,
                                            _RepresentationClass=_RepresentationClass,
                                            _TopologyClass=_TopologyClass)

        builder = _UniformGrid2DBuilder()

        self.args = {
            'dx': dx,
            'dy': dy,
            'nx': nx,
            'ny': ny,
            'origin': origin,
            'overlap': overlap
        }

        builder.buildGridData([dx, dy], [nx, ny], overlap, communicator,
                              origin)

        ([self.dx, self.dy],
         [self.nx, self.ny],
         self.dim,
         scale,
         self.globalNumberOfCells,
         self.globalNumberOfFaces,
         self.overlap,
         self.offset,
         self.numberOfVertices,
         self.numberOfFaces,
         self.numberOfCells,
         self.shape,
         self.physicalShape,
         self._meshSpacing,
         self.numberOfHorizontalRows,
         self.numberOfVerticalColumns,
         self.numberOfHorizontalFaces,
         self.numberOfVerticalFaces,
         self.origin) = builder.gridData

    """
    Topology set and calculate
    """

    @property
    def _exteriorFaces(self):
        """
        Return only the faces that have one neighboring cell.
        """
        exteriorIDs = numerix.concatenate((numerix.arange(0, self.nx),
                                           numerix.arange(0, self.nx) + self.nx * self.ny,
                                           numerix.arange(0, self.ny) * self.numberOfVerticalColumns + self.numberOfHorizontalFaces,
                                           numerix.arange(0, self.ny) * self.numberOfVerticalColumns + self.numberOfHorizontalFaces + self.nx))

        from fipy.variables.faceVariable import FaceVariable
        exteriorFaces = FaceVariable(mesh=self, value=False)
        exteriorFaces[exteriorIDs] = True
        return exteriorFaces

    @property
    def _interiorFaces(self):
        """
        Return only the faces that have two neighboring cells.
        """
        Hids = numerix.arange(0, self.numberOfHorizontalFaces)
        Hids = numerix.reshape(Hids, (self.numberOfHorizontalRows, self.nx))
        Hids = Hids[1:-1, ...]

        Vids = numerix.arange(self.numberOfHorizontalFaces, self.numberOfFaces)
        Vids = numerix.reshape(Vids, (self.ny, self.numberOfVerticalColumns))
        Vids = Vids[..., 1:-1]

        interiorIDs = numerix.concatenate((numerix.reshape(Hids, (self.nx * (self.ny - 1),)),
                                           numerix.reshape(Vids, ((self.nx - 1) * self.ny,))))

        from fipy.variables.faceVariable import FaceVariable
        interiorFaces = FaceVariable(mesh=self, value=False)
        interiorFaces[interiorIDs] = True
        return interiorFaces

    @property
    def _cellToFaceOrientations(self):
        cellFaceOrientations = numerix.ones((4, self.numberOfCells), 'l')
        if self.numberOfCells > 0:
            cellFaceOrientations[0, self.nx:] = -1
            cellFaceOrientations[3,:] = -1
            cellFaceOrientations[3, ::self.nx] = 1
        return cellFaceOrientations

    if inline.doInline:
        @property
        def _adjacentCellIDs(self):
            faceCellIDs0 =  numerix.zeros(self.numberOfFaces, 'l')
            faceCellIDs1 =  numerix.zeros(self.numberOfFaces, 'l')

            inline._runInline("""
                int ID = j * ni + i;

                faceCellIDs0[ID] = ID - ni;
                faceCellIDs1[ID] = ID;

                faceCellIDs0[ID + Nhor + j] = ID - 1;
                faceCellIDs1[ID + Nhor + j] = ID;

                if (j == 0) {
                    faceCellIDs0[ID] = ID;
                }

                if (j == nj - 1) {
                    faceCellIDs0[ID + ni] = ID;
                    faceCellIDs1[ID + ni] = ID;
                }

                if (i == 0) {
                    faceCellIDs0[ID + Nhor + j] = ID;
                }

                if ( i == ni - 1 ) {
                    faceCellIDs0[ID + Nhor + j + 1] = ID;
                    faceCellIDs1[ID + Nhor + j + 1] = ID;
                }

            """,
            Nhor=self.numberOfHorizontalFaces,
            faceCellIDs0=faceCellIDs0,
            faceCellIDs1=faceCellIDs1,
            ni=self.nx,
            nj=self.ny)

            return (faceCellIDs0, faceCellIDs1)

    else:
        @property
        def _adjacentCellIDs(self):
            Hids = numerix.zeros((self.numberOfHorizontalRows, self.nx, 2), 'l')
            indices = numerix.indices((self.numberOfHorizontalRows, self.nx))

            Hids[..., 1] = indices[1] + indices[0] * self.nx
            Hids[..., 0] = Hids[..., 1] - self.nx

            if self.numberOfHorizontalRows > 0:
                Hids[0, ..., 0] = Hids[0, ..., 1]
                Hids[0, ..., 1] = Hids[0, ..., 0]
                Hids[-1, ..., 1] = Hids[-1, ..., 0]

            Vids = numerix.zeros((self.ny, self.numberOfVerticalColumns, 2), 'l')
            indices = numerix.indices((self.ny, self.numberOfVerticalColumns))
            Vids[..., 1] = indices[1] + indices[0] * self.nx
            Vids[..., 0] = Vids[..., 1] - 1

            if self.numberOfVerticalColumns > 0:
                Vids[..., 0, 0] = Vids[..., 0, 1]
                Vids[..., 0, 1] = Vids[..., 0, 0]
                Vids[..., -1, 1] = Vids[..., -1, 0]

            faceCellIDs =  numerix.concatenate((numerix.reshape(Hids, (self.numberOfHorizontalFaces, 2)),
                                                numerix.reshape(Vids, (self.numberOfFaces - self.numberOfHorizontalFaces, 2))))

            return (faceCellIDs[:, 0], faceCellIDs[:, 1])

    @property
    def _cellToCellIDs(self):
        ids = MA.zeros((4, self.nx, self.ny), 'l')
        indices = numerix.indices((self.nx, self.ny))
        ids[0] = indices[0] + (indices[1] - 1) * self.nx
        ids[1] = (indices[0] + 1) + indices[1] * self.nx
        ids[2] = indices[0] + (indices[1] + 1) * self.nx
        ids[3] = (indices[0] - 1) + indices[1] * self.nx

        if self.ny > 0:
            ids[0, ..., 0] = MA.masked
            ids[2, ..., -1] = MA.masked
        if self.nx > 0:
            ids[1, -1, ...] = MA.masked
            ids[3, 0, ...] = MA.masked

        return MA.reshape(ids.swapaxes(1, 2), (4, self.numberOfCells))

    @property
    def _cellToCellIDsFilled(self):
        N = self.numberOfCells
        M = self._maxFacesPerCell
        cellIDs = numerix.repeat(numerix.arange(N)[numerix.newaxis, ...], M, axis=0)
        cellToCellIDs = self._cellToCellIDs
        return MA.where(MA.getmaskarray(cellToCellIDs), cellIDs, cellToCellIDs)

    """
    Geometry set and calculate
    """

    @property
    def _orientedAreaProjections(self):
        return self._areaProjections

    if inline.doInline:
        @property
        def _areaProjections(self):
            areaProjections = numerix.zeros((2, self.numberOfFaces), 'd')

            inline._runInline("""
                              if (i < nx) {
                                  areaProjections[i + 1 * ni] = faceAreas[i] * faceNormals[i + 1 * ni];
                              } else if (i < Nhor) {
                                  areaProjections[i + 1 * ni] = faceAreas[i] * faceNormals[i + 1 * ni];
                              } else if ( (i - Nhor) % (nx + 1) == 0 ) {
                                  areaProjections[i + 0 * ni] = faceAreas[i] * faceNormals[i + 0 * ni];
                              } else {
                                  areaProjections[i + 0 * ni] = faceAreas[i] * faceNormals[i + 0 * ni];
                              }
                              """,
                              dx = float(self.dx), # horrible hack to get around
                              dy = float(self.dy), # http://www.scipy.org/scipy/scipy/ticket/496
                              nx = self.nx,
                              Nhor = self.numberOfHorizontalFaces,
                              areaProjections = areaProjections,
                              ni = self.numberOfFaces,
                              faceNormals = self.faceNormals,
                              faceAreas = self._faceAreas)

            return areaProjections

    else:
        @property
        def _areaProjections(self):
            return self.faceNormals * self._faceAreas

    @property
    def _faceAspectRatios(self):
        return self._faceAreas / self._cellDistances

    @property
    def _faceAreas(self):
        faceAreas = numerix.zeros(self.numberOfFaces, 'd')
        faceAreas[:self.numberOfHorizontalFaces] = self.dx
        faceAreas[self.numberOfHorizontalFaces:] = self.dy
        return faceAreas

    @property
    def faceNormals(self):
        normals = numerix.zeros((2, self.numberOfFaces), 'd')

        normals[1, :self.numberOfHorizontalFaces] = 1
        normals[1, :self.nx] = -1

        normals[0, self.numberOfHorizontalFaces:] = 1
        if self.numberOfVerticalColumns > 0:
            normals[0, self.numberOfHorizontalFaces::self.numberOfVerticalColumns] = -1

        return normals

    @property
    def _cellVolumes(self):
        return numerix.ones(self.numberOfCells, 'd') * self.dx * self.dy

    @property
    def _cellCenters(self):
        centers = numerix.zeros((2, self.nx, self.ny), 'd')
        indices = numerix.indices((self.nx, self.ny))
        centers[0] = (indices[0] + 0.5) * self.dx
        centers[1] = (indices[1] + 0.5) * self.dy
        ccs = centers.reshape((2, self.numberOfCells),
                               order='F') + self.origin
        return ccs

    @property
    def _cellDistances(self):
        Hdis = numerix.repeat((self.dy,), self.numberOfHorizontalFaces)
        Hdis = numerix.reshape(Hdis, (self.nx, self.numberOfHorizontalRows))
        if self.numberOfHorizontalRows > 0:
            Hdis[..., 0] = self.dy / 2.
            Hdis[..., -1] = self.dy / 2.

        Vdis = numerix.repeat((self.dx,), self.numberOfFaces - self.numberOfHorizontalFaces)
        Vdis = numerix.reshape(Vdis, (self.numberOfVerticalColumns, self.ny))
        if self.numberOfVerticalColumns > 0:
            Vdis[0, ...] = self.dx / 2.
            Vdis[-1, ...] = self.dx / 2.

        return numerix.concatenate((numerix.reshape(numerix.swapaxes(Hdis, 0, 1), (self.numberOfHorizontalFaces,)),
                                    numerix.reshape(numerix.swapaxes(Vdis, 0, 1), (self.numberOfFaces - self.numberOfHorizontalFaces,))))

    @property
    def _faceToCellDistanceRatio(self):
        """how far face is from first to second cell
        
        distance from center of face to center of first cell divided by distance
        between cell centers
        """
        faceToCellDistanceRatios = numerix.zeros(self.numberOfFaces, 'd')
        faceToCellDistanceRatios[:] = 0.5
        faceToCellDistanceRatios[:self.nx] = 1.
        faceToCellDistanceRatios[self.numberOfHorizontalFaces - self.nx:self.numberOfHorizontalFaces] = 1.
        if self.numberOfVerticalColumns > 0:
            faceToCellDistanceRatios[self.numberOfHorizontalFaces::self.numberOfVerticalColumns] = 1.
            faceToCellDistanceRatios[(self.numberOfHorizontalFaces + self.nx)::self.numberOfVerticalColumns] = 1.
        return faceToCellDistanceRatios

    def _getFaceToCellDistances(self):
        if hasattr(self, "_internalFaceToCellDistances"):
            """faces have been connected."""
            return self._internalFaceToCellDistances
        else:
            faceToCellDistances = numerix.zeros((2, self.numberOfFaces), 'd')
            distances = self._cellDistances
            ratios = self._faceToCellDistanceRatio
            faceToCellDistances[0] = distances * ratios
            faceToCellDistances[1] = distances * (1 - ratios)
            return faceToCellDistances

    def _setFaceToCellDistances(self, v):
        """Exists only to allow `_connectFaces`."""
        self._internalFaceToCellDistances = v

    _faceToCellDistances = property(_getFaceToCellDistances,
                                    _setFaceToCellDistances)

    @property
    def _faceTangents1(self):
        tangents = numerix.zeros((2, self.numberOfFaces), 'd')

        if self.numberOfFaces > 0:
            tangents[0, :self.numberOfHorizontalFaces] = -1
            tangents[0, :self.nx] = 1
            tangents[1, self.numberOfHorizontalFaces:] = 1
            tangents[1, self.numberOfHorizontalFaces::self.numberOfVerticalColumns] = -1

        return tangents

    @property
    def _faceTangents2(self):
        return numerix.zeros((2, self.numberOfFaces), 'd')

    @property
    def _cellToCellDistances(self):
        distances = numerix.zeros((4, self.nx, self.ny), 'd')
        distances[0] = self.dy
        distances[1] = self.dx
        distances[2] = self.dy
        distances[3] = self.dx

        if self.ny > 0:
            distances[0, ..., 0] = self.dy / 2.
            distances[2, ..., -1] = self.dy / 2.
        if self.nx > 0:
            distances[3, 0, ...] = self.dx / 2.
            distances[1, -1, ...] = self.dx / 2.

        return distances.reshape((4, self.numberOfCells), order='F')


    @property
    def _cellNormals(self):
        normals = numerix.zeros((2, 4, self.numberOfCells), 'd')
        normals[:, 0] = [[ 0], [-1]]
        normals[:, 1] = [[ 1], [ 0]]
        normals[:, 2] = [[ 0], [ 1]]
        normals[:, 3] = [[-1], [ 0]]

        return normals

    @property
    def _cellAreas(self):
        areas = numerix.ones((4, self.numberOfCells), 'd')
        areas[0] = self.dx
        areas[1] = self.dy
        areas[2] = self.dx
        areas[3] = self.dy
        return areas

    @property
    def _cellAreaProjections(self):
        return self._cellAreas * self._cellNormals

    @property
    def _faceCenters(self):
        Hcen = numerix.zeros((2, self.nx, self.numberOfHorizontalRows), 'd')
        indices = numerix.indices((self.nx, self.numberOfHorizontalRows))
        Hcen[0, ...] = (indices[0] + 0.5) * self.dx
        Hcen[1, ...] = indices[1] * self.dy

        Vcen = numerix.zeros((2, self.numberOfVerticalColumns, self.ny), 'd')
        indices = numerix.indices((self.numberOfVerticalColumns, self.ny))
        Vcen[0, ...] = indices[0] * self.dx
        Vcen[1, ...] = (indices[1] + 0.5) * self.dy

        return numerix.concatenate((Hcen.reshape((2, self.numberOfHorizontalFaces), order='F'),
                                    Vcen.reshape((2,
                                        self.numberOfVerticalFaces),
                                        order='F')), axis=1) + self.origin

    def _translate(self, vector):
        return self.__class__(dx = self.args['dx'], nx = self.args['nx'],
                              dy = self.args['dy'], ny = self.args['ny'],
                             origin = numerix.array(self.args['origin']) + vector, overlap=self.args['overlap'])

    def __mul__(self, factor):
        dxdy = numerix.array([[self.args['dx']],
                              [self.args['dy']]])
        dxdy *= factor

        return UniformGrid2D(dx=dxdy[0,0], nx=self.args['nx'],
                             dy=dxdy[1,0], ny=self.args['ny'],
                             origin=numerix.array(self.args['origin']) * factor, overlap=self.args['overlap'])

    @property
    def _concatenableMesh(self):
        from fipy.meshes.nonUniformGrid2D import NonUniformGrid2D
        args = self.args.copy()
        origin = args['origin']
        from fipy.tools import serialComm
        args['communicator'] = serialComm
        del args['origin']
        return NonUniformGrid2D(**args) + origin

    @property
    def _cellFaceIDs(self):
        return _Grid2DBuilder.createCells(self.nx, self.ny,
                                          self.numberOfFaces,
                                          self.numberOfHorizontalFaces,
                                          self.numberOfVerticalColumns)

    @property
    def _maxFacesPerCell(self):
        return 4

    @property
    def vertexCoords(self):
        return _Grid2DBuilder.createVertices(self.nx, self.ny,
                                             self.dx, self.dy,
                                             self.numberOfVertices,
                                             self.numberOfVerticalColumns) \
                 + self.origin

    if inline.doInline:
        @property
        def faceCellIDs(self):
            faceCellIDs = numerix.zeros((2, self.numberOfFaces), 'l')
            mask = numerix.zeros((2, self.numberOfFaces), 'l')

            inline._runInline("""
                int ID = j * ni + i;
                int rowlength = ni * nj + Nhor + nj;

                faceCellIDs[ID + 0 * rowlength] = ID - ni;
                faceCellIDs[ID + 1 * rowlength] = ID;

                faceCellIDs[ID + Nhor + j + 0 * rowlength] = ID - 1;
                faceCellIDs[ID + Nhor + j + 1 * rowlength] = ID;

                if (j == 0) {
                    faceCellIDs[ID + 0 * rowlength] = ID;
                    mask[ID + 1 * rowlength] = 1;
                }

                if (j == nj - 1) {
                    faceCellIDs[ID + ni + 0 * rowlength] = ID;
                    mask[ID + ni + 1 * rowlength] = 1;
                }

                if (i == 0) {
                    faceCellIDs[ID + Nhor + j + 0 * rowlength] = ID;
                    mask[ID + Nhor + j + 1 * rowlength] = 1;
                }

                if ( i == ni - 1 ) {
                    faceCellIDs[ID + Nhor + j + 1 + 0 * rowlength] = ID;
                    mask[ID + Nhor + j + 1 + 1 * rowlength] = 1;
                }
            """,
            Nhor=self.numberOfHorizontalFaces,
            mask=mask,
            faceCellIDs=faceCellIDs,
            ni=self.nx,
            nj=self.ny)

            return MA.masked_where(mask, faceCellIDs)
    else:
        @property
        def faceCellIDs(self):
            Hids = numerix.zeros((2, self.nx, self.numberOfHorizontalRows), 'l')
            indices = numerix.indices((self.nx, self.numberOfHorizontalRows))
            Hids[1] = indices[0] + indices[1] * self.nx
            Hids[0] = Hids[1] - self.nx
            if self.numberOfHorizontalRows > 0:
                Hids[0, ..., 0] = Hids[1, ..., 0]
                Hids[1, ..., 0] = -1
                Hids[1, ..., -1] = -1

            Vids = numerix.zeros((2, self.numberOfVerticalColumns, self.ny), 'l')
            indices = numerix.indices((self.numberOfVerticalColumns, self.ny))
            Vids[1] = indices[0] + indices[1] * self.nx
            Vids[0] = Vids[1] - 1
            if self.numberOfVerticalColumns > 0:
                Vids[0, 0] = Vids[1, 0]
                Vids[1, 0] = -1
                Vids[1, -1] = -1

            return MA.masked_values(numerix.concatenate((Hids.reshape((2, self.numberOfHorizontalFaces), order='F'),
                                                         Vids.reshape((2, self.numberOfFaces - self.numberOfHorizontalFaces), order='F')), axis=1), value = -1)

    @property
    def _cellVertexIDs(self):
        return self._orderedCellVertexIDs

    @property
    def faceVertexIDs(self):
        Hids = numerix.zeros((2, self.nx, self.numberOfHorizontalRows), 'l')
        indices = numerix.indices((self.nx, self.numberOfHorizontalRows))
        Hids[0] = indices[0] + indices[1] * self.numberOfVerticalColumns
        Hids[1] = Hids[0] + 1

        Vids = numerix.zeros((2, self.numberOfVerticalColumns, self.ny), 'l')
        indices = numerix.indices((self.numberOfVerticalColumns, self.ny))
        Vids[0] = indices[0] + indices[1] * self.numberOfVerticalColumns
        Vids[1] = Vids[0] + self.numberOfVerticalColumns

        return numerix.concatenate((Hids.reshape((2, self.numberOfHorizontalFaces), order='F'),
                                    Vids.reshape((2, self.numberOfFaces - self.numberOfHorizontalFaces),
                                                 order='F')),
                                   axis=1)

    def _calcOrderedCellVertexIDs(self):
        """Correct ordering for VTK_PIXEL"""
        ids = numerix.zeros((4, self.nx, self.ny), 'l')
        indices = numerix.indices((self.nx, self.ny))
        ids[2] = indices[0] + (indices[1] + 1) * self.numberOfVerticalColumns
        ids[1] = ids[2] + 1
        ids[3] = indices[0] + indices[1] * self.numberOfVerticalColumns
        ids[0] = ids[3] + 1

        return ids.reshape((4, self.numberOfCells), order='F')

    def _getNearestCellID(self, points):
        """
        Test cases

           >>> from fipy import *
           >>> m = Grid2D(nx=3, ny=2)
           >>> eps = numerix.array([[1e-5, 1e-5]])
           >>> print(m._getNearestCellID(((0., .9, 3.), (0., 2., 2.))))
           [0 3 5]
           >>> print(m._getNearestCellID(([1.1], [1.5])))
           [4]
           >>> m0 = Grid2D(nx=2, ny=2, dx=1., dy=1.)
           >>> m1 = Grid2D(nx=4, ny=4, dx=.5, dy=.5)
           >>> print(m0._getNearestCellID(m1.cellCenters.globalValue))
           [0 0 1 1 0 0 1 1 2 2 3 3 2 2 3 3]

        """
        nx = self.args['nx']
        ny = self.args['ny']

        if nx == 0 or ny == 0:
            return numerix.arange(0)

        x0, y0 = self.cellCenters.globalValue[..., 0]
        xi, yi = points
        dx, dy = self.dx, self.dy

        i = numerix.array(numerix.rint(((xi - x0) / dx)), 'l')
        i[i < 0] = 0
        i[i > nx - 1] = nx - 1

        j = numerix.array(numerix.rint(((yi - y0) / dy)), 'l')
        j[j < 0] = 0
        j[j > ny - 1]  = ny - 1

        return j * nx + i

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.

            >>> dx = 0.5
            >>> dy = 2.
            >>> nx = 3
            >>> ny = 2

            >>> mesh = UniformGrid2D(nx = nx, ny = ny, dx = dx, dy = dy)

            >>> vertices = numerix.array(((0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.),
            ...                           (0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.)))
            >>> vertices *= numerix.array(((dx,), (dy,)))
            >>> print(numerix.allequal(vertices,
            ...                        mesh.vertexCoords)) # doctest: +PROCESSOR_0
            True

            >>> faces = numerix.array([[0, 1, 2, 4, 5, 6, 8, 9, 10, 0, 1, 2, 3,
            ...                         4, 5, 6, 7],
            ...                        [1, 2, 3, 5, 6, 7, 9, 10, 11, 4, 5, 6, 7,
            ...                         8, 9, 10, 11]])
            >>> print(numerix.allequal(faces,
            ...                        mesh.faceVertexIDs)) # doctest: +PROCESSOR_0
            True

            >>> cells = numerix.array(((0, 1, 2, 3, 4, 5),
            ...                        (10, 11, 12, 14, 15, 16),
            ...                        (3, 4, 5, 6, 7, 8),
            ...                        (9, 10, 11, 13, 14, 15)))
            >>> print(numerix.allequal(cells,
            ...                        mesh.cellFaceIDs)) # doctest: +PROCESSOR_0
            True

            >>> externalFaces = numerix.array((0, 1, 2, 6, 7, 8, 9, 12, 13, 16))
            >>> print(numerix.allequal(externalFaces,
            ...                        numerix.nonzero(mesh.exteriorFaces))) # doctest: +PROCESSOR_0
            True

            >>> internalFaces = numerix.array((3, 4, 5, 10, 11, 14, 15))
            >>> print(numerix.allequal(internalFaces,
            ...                        numerix.nonzero(mesh.interiorFaces))) # doctest: +PROCESSOR_0
            True

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values(((0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 0, 1, 2, 3, 3, 4, 5),
            ...                                 (-1, -1, -1, 3, 4, 5, -1, -1, -1, -1, 1, 2, -1, -1, 4, 5, -1)), -1)
            >>> print(numerix.allequal(faceCellIds, mesh.faceCellIDs)) # doctest: +PROCESSOR_0
            True

            >>> faceAreas = numerix.array((dx, dx, dx, dx, dx, dx, dx, dx, dx,
            ...                            dy, dy, dy, dy, dy, dy, dy, dy))
            >>> print(numerix.allclose(faceAreas, mesh._faceAreas, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> faceCoords = numerix.take(vertices, faces, axis=1)
            >>> faceCenters = (faceCoords[..., 0,:] + faceCoords[..., 1,:]) / 2.
            >>> print(numerix.allclose(faceCenters, mesh.faceCenters, atol = 1e-10, rtol = 1e-10))
            True

            >>> faceNormals = numerix.array(((0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 1., 1., 1., -1., 1., 1., 1.),
            ...                              (-1., -1., -1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0.)))
            >>> print(numerix.allclose(faceNormals, mesh.faceNormals, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> cellToFaceOrientations = numerix.array(((1,  1,  1, -1, -1, -1),
            ...                                         (1,  1,  1,  1,  1,  1),
            ...                                         (1,  1,  1,  1,  1,  1),
            ...                                         (1, -1, -1,  1, -1, -1)))
            >>> print(numerix.allequal(cellToFaceOrientations, mesh._cellToFaceOrientations)) # doctest: +PROCESSOR_0
            True

            >>> cellVolumes = numerix.array((dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy))
            >>> print(numerix.allclose(cellVolumes, mesh.cellVolumes, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> cellCenters = numerix.array(((dx/2., 3.*dx/2., 5.*dx/2., dx/2., 3.*dx/2., 5.*dx/2.),
            ...                              (dy/2., dy/2., dy/2., 3.*dy/2., 3.*dy/2., 3.*dy/2.)))
            >>> print(numerix.allclose(cellCenters, mesh.cellCenters, atol = 1e-10, rtol = 1e-10))
            True

            >>> faceToCellDistances = MA.masked_values(((dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2.),
            ...                                         (-1, -1, -1, dy / 2., dy / 2., dy / 2., -1, -1, -1, -1, dx / 2., dx / 2., -1, -1, dx / 2., dx / 2., -1)), -1)
            >>> print(numerix.allclose(faceToCellDistances, mesh._faceToCellDistances, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> cellDistances = numerix.array((dy / 2., dy / 2., dy / 2.,
            ...                                dy, dy, dy,
            ...                                dy / 2., dy / 2., dy / 2.,
            ...                                dx / 2., dx, dx,
            ...                                dx / 2.,
            ...                                dx / 2., dx, dx,
            ...                                dx / 2.))
            >>> print(numerix.allclose(cellDistances, mesh._cellDistances, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> print(numerix.allclose(faceToCellDistanceRatios, mesh._faceToCellDistanceRatio, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> areaProjections = faceNormals * faceAreas
            >>> print(numerix.allclose(areaProjections, mesh._areaProjections, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> tangents1 = numerix.array(((1., 1., 1., -1., -1., -1., -1., -1., -1., 0., 0., 0., 0., 0., 0., 0., 0.),
            ...                            (0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 1., 1., 1., -1., 1., 1., 1.)))
            >>> print(numerix.allclose(tangents1, mesh._faceTangents1, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> tangents2 = numerix.zeros((2, 17), 'd')
            >>> print(numerix.allclose(tangents2, mesh._faceTangents2, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> cellToCellIDs = MA.masked_values(((-1, -1, -1, 0, 1, 2),
            ...                                   (1, 2, -1, 4, 5, -1),
            ...                                   (3, 4, 5, -1, -1, -1),
            ...                                   (-1, 0, 1, -1, 3, 4)), -1)
            >>> print(numerix.allequal(cellToCellIDs, mesh._cellToCellIDs)) # doctest: +PROCESSOR_0
            True

            >>> cellToCellDistances = MA.masked_values(((dy / 2., dy / 2., dy / 2.,      dy,      dy,      dy),
            ...                                         (     dx,      dx, dx / 2.,      dx,      dx, dx / 2.),
            ...                                         (     dy,      dy,      dy, dy / 2., dy / 2., dy / 2.),
            ...                                         (dx / 2.,      dx,      dx, dx / 2.,      dx,      dx)), -1)
            >>> print(numerix.allclose(cellToCellDistances, mesh._cellToCellDistances, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> cellNormals = numerix.array(((( 0,  0,  0,  0,  0,  0),
            ...                               ( 1,  1,  1,  1,  1,  1),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               (-1, -1, -1, -1, -1, -1)),
            ...                              ((-1, -1, -1, -1, -1, -1),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 1,  1,  1,  1,  1,  1),
            ...                               ( 0,  0,  0,  0,  0,  0))))
            >>> print(numerix.allclose(cellNormals, mesh._cellNormals, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> cellAreaProjections = numerix.array((((  0,  0,  0,  0,  0,  0),
            ...                                       ( dy, dy, dy, dy, dy, dy),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       (-dy, -dy, -dy, -dy, -dy, -dy)),
            ...                                      ((-dx, -dx, -dx, -dx, -dx, -dx),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       ( dx, dx, dx, dx, dx, dx),
            ...                                       (  0,  0,  0,  0,  0,  0))))
            >>> print(numerix.allclose(cellAreaProjections, mesh._cellAreaProjections, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> cellVertexIDs = MA.masked_array(((5, 6, 7, 9, 10, 11),
            ...                                  (4, 5, 6, 8, 9, 10),
            ...                                  (1, 2, 3, 5, 6, 7),
            ...                                  (0, 1, 2, 4, 5, 6)), -1000)

            >>> print(numerix.allclose(mesh._cellVertexIDs, cellVertexIDs)) # doctest: +PROCESSOR_0
            True

            >>> from fipy.tools import dump
            >>> (f, filename) = dump.write(mesh, extension = '.gz')
            >>> unpickledMesh = dump.read(filename, f)

            >>> print(numerix.allequal(mesh.cellCenters, unpickledMesh.cellCenters))
            True

            >>> faceVertexIDs = [[ 0, 1, 2, 4, 5, 6, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7],
            ...                  [ 1, 2, 3, 5, 6, 7, 9, 10, 11, 4, 5, 6, 7, 8, 9, 10, 11]]
            >>> print(numerix.allequal(mesh.faceVertexIDs, faceVertexIDs)) # doctest: +PROCESSOR_0
            True

            >>> mesh = UniformGrid2D(nx=3)
            >>> print(numerix.allequal(mesh._adjacentCellIDs[0],
            ...                        [0, 1, 2, 0, 1, 2, 0, 0, 1, 2])) # doctest: +PROCESSOR_0
            True
            >>> print(numerix.allequal(mesh._adjacentCellIDs[1],
            ...                        [0, 1, 2, 0, 1, 2, 0, 1, 2, 2])) # doctest: +PROCESSOR_0
            True

            >>> faceCellIDs = [[0, 1, 2, 0, 1, 2, 0, 0, 1, 2],
            ...                [-1, -1, -1, -1, -1, -1, -1, 1, 2, -1]]
            >>> print(numerix.allequal(mesh.faceCellIDs.filled(-1),
            ...                        faceCellIDs)) # doctest: +PROCESSOR_0
            True

            >>> mesh = UniformGrid2D(ny=3)
            >>> print(numerix.allequal(mesh._adjacentCellIDs[0],
            ...                        [0, 0, 1, 2, 0, 0, 1, 1, 2, 2])) # doctest: +PROCESSOR_0
            True
            >>> print(numerix.allequal(mesh._adjacentCellIDs[1],
            ...                        [0, 1, 2, 2, 0, 0, 1, 1, 2, 2])) # doctest: +PROCESSOR_0
            True
            >>> faceCellIDs = [[0, 0, 1, 2, 0, 0, 1, 1, 2, 2],
            ...                [-1, 1, 2, -1, -1, -1, -1, -1, -1, -1]]
            >>> print(numerix.allequal(mesh.faceCellIDs.filled(-1),
            ...                        faceCellIDs)) # doctest: +PROCESSOR_0
            True

        Following test added to change `nx`, `ny` argument to integer when its a float to prevent
        warnings from the solver.

            >>> from fipy import *
            >>> mesh = UniformGrid2D(nx=3, ny=3, dx=1., dy=1.)
            >>> var = CellVariable(mesh=mesh)
            >>> DiffusionTerm().solve(var)

        Ensure that ghost faces are excluded from accumulating operations
        (#856).  Four exterior surfaces of :math:`10\times 10` square mesh
        should each have a total area of 10, regardless of partitioning.

            >>> square = UniformGrid2D(nx=10, dx=1., ny=10, dy=1.)

            >>> area = (square._faceAreas * square.facesBottom).sum()
            >>> print(numerix.allclose(area, 10))
            True

            >>> area = (square._faceAreas * square.facesTop).sum()
            >>> print(numerix.allclose(area, 10))
            True

            >>> area = (square._faceAreas * square.facesLeft).sum()
            >>> print(numerix.allclose(area, 10))
            True

            >>> area = (square._faceAreas * square.facesRight).sum()
            >>> print(numerix.allclose(area, 10))
            True

        Ensure no double-counting of faces on boundary between partitions.

            >>> area = (square._faceAreas * (square.faceCenters[1] == 5.)).sum()
            >>> print(numerix.allclose(area, 10)) # doctest: +PARALLEL_2
            True

        Size of global value should not depend on number of processors (#400)

            >>> print(square.cellCenters.globalValue.shape)
            (2, 100)
            >>> print(square.faceCenters.globalValue.shape)
            (2, 220)
        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


