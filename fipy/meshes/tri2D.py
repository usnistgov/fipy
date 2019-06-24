from __future__ import unicode_literals
__docformat__ = "restructuredtext"

from fipy.tools import numerix
from fipy.tools import vector
from fipy.tools.dimensions.physicalField import PhysicalField

from fipy.meshes.mesh2D import Mesh2D
from fipy.meshes.representations.gridRepresentation import _Grid2DRepresentation
from fipy.meshes.topologies.meshTopology import _Mesh2DTopology

__all__ = ["Tri2D"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Tri2D(Mesh2D):
    """
    This class creates a mesh made out of triangles.  It does this by
    starting with a standard Cartesian mesh (`Grid2D`) and dividing each cell
    in that mesh (hereafter referred to as a "box") into four equal
    parts with the dividing lines being the diagonals.
    """

    def __init__(self, dx = 1., dy = 1., nx = 1, ny = 1,
                 _RepresentationClass=_Grid2DRepresentation, _TopologyClass=_Mesh2DTopology):
        """
        Creates a 2D triangular mesh with horizontal faces numbered first then
        vertical faces, then diagonal faces.  Vertices are numbered starting
        with the vertices at the corners of boxes and then the vertices at the
        centers of boxes.  Cells on the right of boxes are numbered first, then
        cells on the top of boxes, then cells on the left of boxes, then cells
        on the bottom of boxes.  Within each of the "sub-categories" in the
        above, the vertices, cells and faces are numbered in the usual way.

        Parameters
        ----------
        float dx, dy : float
            The X and Y dimensions of each "box".
            If `dx` <> `dy`, the line segments connecting the cell
            centers will not be orthogonal to the faces.
        int nx, ny : int
            The number of boxes in the X direction and the Y direction.
            The total number of boxes will be equal to `nx * ny`, and the total
            number of cells will be equal to `4 * nx * ny`.
        """

        self.args = {
            'dx': dx,
            'dy': dy,
            'nx': nx,
            'ny': ny
        }

        self.nx = nx
        self.ny = ny

        self.numberOfHorizontalFaces   = self.nx * (self.ny + 1)
        self.numberOfVerticalFaces     = self.ny * (self.nx + 1)
        self.numberOfEachDiagonalFaces = self.nx * self.ny

        self.dx  = PhysicalField(value = dx)
        scale    = PhysicalField(value = 1, unit = self.dx.unit)
        self.dx /= scale

        self.dy = PhysicalField(value = dy)
        if self.dy.unit.isDimensionless():
            self.dy = dy
        else:
            self.dy /= scale

        self.numberOfCornerVertices = (self.nx + 1) * (self. ny + 1)
        self.numberOfCenterVertices = self.nx * self.ny
        self.numberOfTotalVertices  = self.numberOfCornerVertices + self.numberOfCenterVertices

        self.offset = (0, 0)

        vertices = self._createVertices()
        faces    = self._createFaces()

        cells = self._createCells()
        cells = numerix.sort(cells, axis=0)

        Mesh2D.__init__(self, vertices, faces, cells,
                        _RepresentationClass=_RepresentationClass, _TopologyClass=_TopologyClass)

        self.scale = scale

    def _createVertices(self):

        x = numerix.arange(self.nx + 1) * self.dx
        y = numerix.arange(self.ny + 1) * self.dy
        x = numerix.resize(x, (self.numberOfCornerVertices,))
        y = numerix.repeat(y, self.nx + 1)
        boxCorners = numerix.array((x, y))
        x = numerix.arange(0.5, self.nx + 0.5) * self.dx
        y = numerix.arange(0.5, self.ny + 0.5) * self.dy
        x = numerix.resize(x, (self.numberOfCenterVertices,))
        y = numerix.repeat(y, self.nx)
        boxCenters = numerix.array((x, y))
        return numerix.concatenate((boxCorners, boxCenters), axis=1)

    def _createFaces(self):
        ## v1, v2 refer to the cells.
        ## Horizontal faces are first
        v1 = numerix.arange(self.numberOfCornerVertices)
        v2 = v1 + 1
        horizontalFaces = vector.prune(numerix.array((v1, v2)), self.nx + 1, self.nx, axis=1)
        v1 = numerix.arange(self.numberOfCornerVertices - (self.nx + 1))
        v2 = v1 + self.nx + 1
        verticalFaces =  numerix.array((v1, v2))

        ## reverse some of the face orientations to obtain the correct normals

        tmp = horizontalFaces.copy()
        horizontalFaces[0, :self.nx] = tmp[1, :self.nx]
        horizontalFaces[1, :self.nx] = tmp[0, :self.nx]

        tmp = verticalFaces.copy()
        verticalFaces[0] = tmp[1]
        verticalFaces[1] = tmp[0]
        verticalFaces[0, ::(self.nx + 1)] = tmp[0, ::(self.nx + 1)]
        verticalFaces[1, ::(self.nx + 1)] = tmp[1, ::(self.nx + 1)]

        ## do the center ones now

        cellCenters = numerix.arange(self.numberOfCornerVertices, self.numberOfTotalVertices)
        lowerLefts = vector.prune(numerix.arange(self.numberOfCornerVertices - (self.nx + 1)), self.nx + 1, self.nx)
        lowerRights = lowerLefts + 1
        upperLefts = lowerLefts + self.nx + 1
        upperRights = lowerLefts + self.nx + 2
        lowerLeftFaces = numerix.array((cellCenters, lowerLefts))
        lowerRightFaces = numerix.array((lowerRights, cellCenters))
        upperLeftFaces = numerix.array((cellCenters, upperLefts))
        upperRightFaces = numerix.array((cellCenters, upperRights))
        return numerix.concatenate((horizontalFaces, verticalFaces, lowerLeftFaces, lowerRightFaces, upperLeftFaces, upperRightFaces), axis=1)

    def _createCells(self):
        ## cells = (f1, f2, f3, f4) going anticlockwise.
        ## f1 etc. refer to the faces
        bottomFaces = numerix.arange(0, self.numberOfHorizontalFaces - self.nx)
        topFaces = numerix.arange(self.nx, self.numberOfHorizontalFaces)
        leftFaces = vector.prune(numerix.arange(self.numberOfHorizontalFaces, self.numberOfHorizontalFaces + self.numberOfVerticalFaces), self.nx + 1, self.nx)
        rightFaces = vector.prune(numerix.arange(self.numberOfHorizontalFaces, self.numberOfHorizontalFaces + self.numberOfVerticalFaces), self.nx + 1, 0)
        lowerLeftDiagonalFaces = numerix.arange(self.numberOfHorizontalFaces + self.numberOfVerticalFaces, self.numberOfHorizontalFaces + self.numberOfVerticalFaces + self.numberOfEachDiagonalFaces)
        lowerRightDiagonalFaces = lowerLeftDiagonalFaces + self.numberOfEachDiagonalFaces
        upperLeftDiagonalFaces = lowerRightDiagonalFaces + self.numberOfEachDiagonalFaces
        upperRightDiagonalFaces = upperLeftDiagonalFaces + self.numberOfEachDiagonalFaces
        ##faces in arrays, now get the cells
        bottomOfBoxCells = numerix.array([bottomFaces, lowerRightDiagonalFaces, lowerLeftDiagonalFaces])
        rightOfBoxCells = numerix.array([rightFaces, upperRightDiagonalFaces, lowerRightDiagonalFaces])
        topOfBoxCells = numerix.array([topFaces, upperLeftDiagonalFaces, upperRightDiagonalFaces])
        leftOfBoxCells = numerix.array([leftFaces, lowerLeftDiagonalFaces, upperLeftDiagonalFaces])
        return numerix.concatenate((rightOfBoxCells, topOfBoxCells, leftOfBoxCells, bottomOfBoxCells), axis=1)

    @property
    def physicalShape(self):
        """Return physical dimensions of `Grid2D`.
        """
        return PhysicalField(value = (self.nx * self.dx * self.scale,
                                      self.ny * self.dy * self.scale))

    @property
    def _meshSpacing(self):
        return numerix.array((self.dx, self.dy))[..., numerix.newaxis]

    @property
    def shape(self):
        return (self.nx, self.ny)

    @property
    def _isOrthogonal(self):
        return True

## pickling

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.

            >>> dx = 0.5
            >>> dy = 2.
            >>> nx = 3
            >>> ny = 2

            >>> mesh = Tri2D(nx = nx, ny = ny, dx = dx, dy = dy)

            >>> vertices = numerix.array(((0.0, 0.5, 1.0, 1.5, 0.0, 0.5, 1.0, 1.5, 0.0, 0.5, 1.0, 1.5, 0.25, 0.75, 1.25, 0.25, 0.75, 1.25),
            ...                           (0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 4.0, 4.0, 4.0, 4.0, 1.0,  1.0,  1.0,  3.0,  3.0,  3.0)))

            >>> from fipy.tools import numerix
            >>> numerix.allequal(vertices, mesh._createVertices())
            1

            >>> faces = numerix.array(((1, 2, 3, 4, 5, 6, 8,  9, 10, 0, 5, 6, 7, 4, 9, 10, 11, 12, 13, 14, 15, 16, 17,  1,  2,  3,  5,  6,  7, 12, 13, 14, 15, 16, 17, 12, 13, 14, 15, 16, 17),
            ...                        (0, 1, 2, 5, 6, 7, 9, 10, 11, 4, 1, 2, 3, 8, 5,  6,  7,  0,  1,  2,  4,  5,  6, 12, 13, 14, 15, 16, 17,  4,  5,  6,  8,  9, 10,  5,  6,  7,  9, 10, 11)))
            >>> numerix.allequal(faces, mesh._createFaces())
            1

            >>> cells = numerix.array(((10, 11, 12, 14, 15, 16,  3,  4,  5,  6,  7,  8,  9, 10, 11, 13, 14, 15,  0,  1,  2,  3,  4,  5),
            ...                        (35, 36, 37, 38, 39, 40, 29, 30, 31, 32, 33, 34, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28),
            ...                        (23, 24, 25, 26, 27, 28, 35, 36, 37, 38, 39, 40, 29, 30, 31, 32, 33, 34, 17, 18, 19, 20, 21, 22)))
            >>> numerix.allequal(cells, mesh._createCells())
            1

            >>> externalFaces = numerix.array((0, 1, 2, 6, 7, 8, 9, 12, 13, 16))
            >>> print(numerix.allequal(externalFaces,
            ...                        numerix.nonzero(mesh.exteriorFaces)))
            1

            >>> internalFaces = numerix.array((3, 4, 5, 10, 11, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40))
            >>> print(numerix.allequal(internalFaces,
            ...                        numerix.nonzero(mesh.interiorFaces)))
            1

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values(((18, 19, 20,  6,  7,  8,  9, 10, 11, 12,  0,  1,  2, 15,  3,  4,  5, 12, 13, 14, 15, 16, 17,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,  0,  1,  2,  3,  4,  5),
            ...                                 (-1, -1, -1, 21, 22, 23, -1, -1, -1, -1, 13, 14, -1, -1, 16, 17, -1, 18, 19, 20, 21, 22, 23, 18, 19, 20, 21, 22, 23, 12, 13, 14, 15, 16, 17,  6,  7,  8,  9, 10, 11)), -1)
            >>> numerix.allequal(faceCellIds, mesh.faceCellIDs)
            1

            >>> d = (numerix.sqrt((dx*dx)+(dy*dy))) / 2.0 ## length of diagonal edges
            >>> faceAreas = numerix.array((0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
            ...                            2, 2, 2, 2, 2, 2, 2, 2,
            ...                            d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d))
            >>> numerix.allclose(faceAreas, mesh._faceAreas, atol = 1e-10, rtol = 1e-10)
            1

            >>> faceCoords = numerix.take(vertices, faces, axis=1)
            >>> faceCenters = (faceCoords[..., 0,:] + faceCoords[..., 1,:]) / 2.
            >>> print(numerix.allclose(faceCenters, mesh.faceCenters, atol = 1e-10, rtol = 1e-10))
            True

            >>> xc = dy  / numerix.sqrt((dx * dx) + (dy * dy))
            >>> yc = dx  / numerix.sqrt((dx * dx) + (dy * dy))

            >>> faceNormals = numerix.array((( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, xc, xc, xc, xc, xc, xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc),
            ...                              (-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, yc, yc, yc, yc, yc, yc)))
            >>> numerix.allclose(faceNormals, mesh.faceNormals, atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = numerix.array((( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, 1, -1, -1, 1, 1, 1, -1, -1, -1),
            ...                                         ( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1),
            ...                                         ( 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)))
            >>> numerix.allequal(cellToFaceOrientations, mesh._cellToFaceOrientations)
            1

            >>> cellVolumes = numerix.array((0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
            ...                              0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
            ...                              0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25))
            >>> numerix.allclose(cellVolumes, mesh.cellVolumes, atol = 1e-10, rtol = 1e-10)
            1

            >>> sixth = 1.0 / 6.0
            >>> cellCenters = numerix.array(((  5*sixth, 11*sixth, 17*sixth,  5*sixth, 11*sixth, 17*sixth,      0.5,      1.5,      2.5,      0.5,      1.5,      2.5,  1*sixth,  7*sixth, 13*sixth,  1*sixth,  7*sixth, 13*sixth,      0.5,      1.5,      2.5,      0.5,      1.5,      2.5),
            ...                              (      0.5,      0.5,      0.5,      1.5,      1.5,      1.5,  5*sixth,  5*sixth,  5*sixth, 11*sixth, 11*sixth, 11*sixth,      0.5,      0.5,      0.5,      1.5,      1.5,      1.5,  1*sixth,  1*sixth,  1*sixth,  7*sixth,  7*sixth,  7*sixth)))
            >>> cellCenters *= numerix.array([[dx], [dy]])
            >>> print(numerix.allclose(cellCenters, mesh.cellCenters, atol = 1e-10, rtol = 1e-10))
            True

            >>> yd = numerix.sqrt(((dx/12.0)*(dx/12.0)) + ((dy/ 4.0)*(dy/ 4.0)))
            >>> xd = numerix.sqrt(((dx/ 4.0)*(dx/ 4.0)) + ((dy/12.0)*(dy/12.0)))
            >>> faceToCellDistances = MA.masked_values(((dy/6.0, dy/6.0, dy/6.0, dy/6.0, dy/6.0, dy/6.0,  dy/6.0, dy/6.0, dy/6.0, dx/6.0, dx/6.0, dx/6.0, dx/6.0, dx/6.0, dx/6.0, dx/6.0, dx/6.0,     yd,     yd,     yd,     yd,     yd,     yd,     yd,     yd,     yd,     yd,     yd,     yd,     xd,     xd,     xd,     xd,     xd,     xd,     yd,     yd,     yd,     yd,     yd,     yd),
            ...                                         (    -1,     -1,     -1, dy/6.0, dy/6.0, dy/6.0,      -1,     -1,     -1,     -1, dx/6.0, dx/6.0,     -1,     -1, dx/6.0, dx/6.0,     -1,     xd,     xd,     xd,     xd,     xd,     xd,     xd,     xd,     xd,     xd,     xd,     xd,     yd,     yd,     yd,     yd,     yd,     yd,     xd,     xd,     xd,     xd,     xd,     xd)), -1)
            >>> numerix.allclose(faceToCellDistances, mesh._faceToCellDistances, atol = 1e-10, rtol = 1e-10)
            1

            >>> dd = numerix.sqrt((dx*dx)+(dy*dy)) / 3.0
            >>> cellDistances = numerix.array((dy/6.0, dy/6.0, dy/6.0, dy/3.0, dy/3.0, dy/3.0, dy/6.0, dy/6.0, dy/6.0,
            ...                                dx/6.0, dx/3.0, dx/3.0, dx/6.0, dx/6.0, dx/3.0, dx/3.0, dx/6.0,
            ...                                dd, dd, dd, dd, dd, dd, dd, dd, dd, dd, dd, dd,
            ...                                dd, dd, dd, dd, dd, dd, dd, dd, dd, dd, dd, dd))
            >>> numerix.allclose(cellDistances, mesh._cellDistances, atol = 1e-10, rtol = 1e-10)
            1

            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> numerix.allclose(faceToCellDistanceRatios, mesh._faceToCellDistanceRatio, atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas
            >>> numerix.allclose(areaProjections, mesh._areaProjections, atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents1 = numerix.array(((1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  yc,  yc,  yc,  yc,  yc,  yc,  yc,  yc,  yc,  yc,  yc,  yc,  yc,  yc,  yc,  yc,  yc,  yc, -yc, -yc, -yc, -yc, -yc, -yc),
            ...                            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0,  xc,  xc,  xc,  xc,  xc,  xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc)))
            >>> numerix.allclose(tangents1, mesh._faceTangents1, atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents2 = numerix.zeros((2, 41), 'l')
            >>> numerix.allclose(tangents2, mesh._faceTangents2, atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((13, 14, -1, 16, 17, -1, 21, 22, 23, -1, -1, -1, -1,  0,  1, -1,  3,  4, -1, -1, -1,  6,  7,  8),
            ...                                   (18, 19, 20, 21, 22, 23, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 12, 13, 14, 15, 16, 17),
            ...                                   ( 6,  7,  8,  9, 10, 11,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,  0,  1,  2,  3,  4,  5)), -1)
            >>> numerix.allequal(cellToCellIDs, mesh._cellToCellIDs)
            1

            >>> cellToCellDistances = numerix.array(((dx/3.0, dx/3.0, dx/6.0, dx/3.0, dx/3.0, dx/6.0, dy/3.0, dy/3.0, dy/3.0, dy/6.0, dy/6.0, dy/6.0, dx/6.0, dx/3.0, dx/3.0, dx/6.0, dx/3.0, dx/3.0, dy/6.0, dy/6.0, dy/6.0, dy/3.0, dy/3.0, dy/3.0),
            ...                                      (    dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd),
            ...                                      (    dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd,     dd)))
            >>> numerix.allclose(cellToCellDistances, mesh._cellToCellDistances, atol = 1e-10, rtol = 1e-10)
            1

            >>> interiorCellIDs = numerix.array((0, 1, 3, 4, 6, 7, 8, 13, 14, 16, 17, 21, 22, 23))
            >>> numerix.allequal(interiorCellIDs, mesh._interiorCellIDs)
            1

            >>> exteriorCellIDs = numerix.array((2, 5, 9, 10, 11, 12, 15, 18, 19, 20))
            >>> numerix.allequal(exteriorCellIDs, mesh._exteriorCellIDs)
            1

            >>> cellNormals = numerix.array((((  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0),
            ...                               (-xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, -xc, xc, xc, xc, xc, xc, xc, -xc, -xc, -xc, -xc, -xc, -xc),
            ...                               (-xc, -xc, -xc, -xc, -xc, -xc, xc, xc, xc, xc, xc, xc, xc, xc, xc, xc, xc, xc, xc, xc, xc, xc, xc, xc)),
            ...                              ((  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1),
            ...                               (-yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, -yc, yc, yc, yc, yc, yc, yc),
            ...                               ( yc, yc, yc, yc, yc, yc, -yc, -yc, -yc, -yc, -yc, -yc, yc, yc, yc, yc, yc, yc, yc, yc, yc, yc, yc, yc))))
            >>> numerix.allclose(cellNormals, mesh._cellNormals, atol = 1e-10, rtol = 1e-10)
            1

            >>> cellAreaProjections = numerix.array((((     dy,     dy,     dy,     dy,     dy,     dy,      0,      0,      0,      0,      0,      0,    -dy,    -dy,    -dy,    -dy,    -dy,    -dy,      0,      0,      0,      0,      0,      0),
            ...                                       ( -dy/2., -dy/2., -dy/2., -dy/2., -dy/2., -dy/2., -dy/2., -dy/2., -dy/2., -dy/2., -dy/2., -dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2., -dy/2., -dy/2., -dy/2., -dy/2., -dy/2., -dy/2.),
            ...                                       ( -dy/2., -dy/2., -dy/2., -dy/2., -dy/2., -dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.,  dy/2.)),
            ...                                      ((      0,      0,      0,      0,      0,      0,     dx,     dx,     dx,     dx,     dx,     dx,      0,      0,      0,      0,      0,      0,    -dx,    -dx,    -dx,    -dx,    -dx,    -dx),
            ...                                       ( -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2.),
            ...                                       (  dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2., -dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2.,  dx/2.))))
            >>> numerix.allclose(cellAreaProjections, mesh._cellAreaProjections, atol = 1e-10, rtol = 1e-10)
            1

            >>> tmp1 = numerix.array((12, 5, 1))
            >>> tmp2 = numerix.array((12, 5, 4))
            >>> tmp3 = numerix.array((12, 4, 0))
            >>> tmp4 = numerix.array((12, 1, 0))
            >>> tmp5 = numerix.array((0, 1, 1))
            >>> cellVertexIDs = numerix.array((tmp1, tmp1 + 1, tmp1 + 2, tmp1 + 3 + tmp5, tmp1 + 4 + tmp5, tmp1 + 5 + tmp5,
            ...                                tmp2, tmp2 + 1, tmp2 + 2, tmp2 + 3 + tmp5, tmp2 + 4 + tmp5, tmp2 + 5 + tmp5,
            ...                                tmp3, tmp3 + 1, tmp3 + 2, tmp3 + 3 + tmp5, tmp3 + 4 + tmp5, tmp3 + 5 + tmp5,
            ...                                tmp4, tmp4 + 1, tmp4 + 2, tmp4 + 3 + tmp5, tmp4 + 4 + tmp5, tmp4 + 5 + tmp5))
            >>> cellVertexIDs = cellVertexIDs.swapaxes(0, 1)
            >>> numerix.allclose(mesh._cellVertexIDs, cellVertexIDs)
            1


            >>> from fipy.tools import dump
            >>> (f, filename) = dump.write(mesh, extension = '.gz')
            >>> unpickledMesh = dump.read(filename, f)

            >>> print(numerix.allequal(mesh.cellCenters, unpickledMesh.cellCenters))
            True
        """

## test test test

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


