from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools import parallelComm

from fipy.meshes.mesh import Mesh
from fipy.meshes.builders import _NonuniformGrid3DBuilder
from fipy.meshes.representations.gridRepresentation import _Grid3DRepresentation
from fipy.meshes.topologies.gridTopology import _Grid3DTopology

__all__ = ["NonUniformGrid3D"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class NonUniformGrid3D(Mesh):
    """
    3D rectangular-prism Mesh

    X axis runs from left to right.
    Y axis runs from bottom to top.
    Z axis runs from front to back.

    Numbering System:

    Vertices: Numbered in the usual way. X coordinate changes most quickly, then Y, then Z.

    Cells: Same numbering system as vertices.

    Faces: XY faces numbered first, then XZ faces, then YZ faces. Within each subcategory, it is numbered in the usual way.
    """
    def __init__(self, dx = 1., dy = 1., dz = 1., nx = None, ny = None, nz = None, overlap=2, communicator=parallelComm,
                 _RepresentationClass=_Grid3DRepresentation, _TopologyClass=_Grid3DTopology):

        builder = _NonuniformGrid3DBuilder()

        self.args = {
            'dx': dx,
            'dy': dy,
            'dz': dz,
            'nx': nx,
            'ny': ny,
            'nz': nz,
            'overlap': overlap,
        }

        if self.args['nx'] is None:
            self.args['nx'] = len(self.args['dx'])

        if self.args['ny'] is None:
            self.args['ny'] = len(self.args['dy'])

        if self.args['nz'] is None:
            self.args['nz'] = len(self.args['dz'])

        builder.buildGridData([dx, dy, dz], [nx, ny, nz], overlap,
                              communicator)

        ([self.dx, self.dy, self.dz],
         [self.nx, self.ny, self.nz],
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
         self.numberOfXYFaces,
         self.numberOfXZFaces,
         self.numberOfYZFaces,
         self.numberOfHorizontalRows,
         self.numberOfVerticalColumns,
         self.numberOfLayersDeep,
         vertices,
         faces,
         cells,
         self.Xoffset, self.Yoffset, self.Zoffset) = builder.gridData

        Mesh.__init__(self, vertices, faces, cells, communicator=communicator,
                      _RepresentationClass=_RepresentationClass, _TopologyClass=_TopologyClass)

        self._setScale(scaleLength = scale)

    def _calcScaleArea(self):
        return self.scale['length']**2

    def _calcScaleVolume(self):
        return self.scale['length']**3

    def _calcFaceNormals(self):
        XYFaceNormals = numerix.zeros((3, self.numberOfXYFaces), 'l')
        XYFaceNormals[2, (self.nx * self.ny):] = 1
        XYFaceNormals[2, :(self.nx * self.ny)] = -1
        XZFaceNormals = numerix.zeros((3, self.numberOfXZFaces), 'l')
        xzd = numerix.arange(self.numberOfXZFaces)
        xzd = xzd % (self.nx * (self.ny + 1))
        xzd = (xzd < self.nx)
        xzd = 1 - (2 * xzd)
        XZFaceNormals[1,:] = xzd
        YZFaceNormals = numerix.zeros((3, self.numberOfYZFaces), 'l')
        YZFaceNormals[0,:] = 1
        YZFaceNormals[0, ::self.nx + 1] = -1
        return numerix.concatenate((XYFaceNormals,
                                    XZFaceNormals,
                                    YZFaceNormals),
                                   axis=-1)

    def _calcFaceTangents(self):
        ## need to see whether order matters.
        faceTangents1 = numerix.zeros((3, self.numberOfFaces), 'd')
        faceTangents2 = numerix.zeros((3, self.numberOfFaces), 'd')
        ## XY faces
        faceTangents1[0, :self.numberOfXYFaces] = 1.
        faceTangents2[1, :self.numberOfXYFaces] = 1.
        ## XZ faces
        faceTangents1[0, self.numberOfXYFaces:self.numberOfXYFaces + self.numberOfXZFaces] = 1.
        faceTangents2[2, self.numberOfXYFaces:self.numberOfXYFaces + self.numberOfXZFaces] = 1.
        ## YZ faces
        faceTangents1[1, self.numberOfXYFaces + self.numberOfXZFaces:] = 1.
        faceTangents2[2, self.numberOfXYFaces + self.numberOfXZFaces:] = 1.
        return faceTangents1, faceTangents2

## The following method is broken when dx, dy or dz are not scalar. Simpler to use the generic
## _calcFaceAreas rather than do the required type checking, resizing and outer product.
##
##     def _calcFaceAreas(self):
##         XYFaceAreas = numerix.ones(self.numberOfXYFaces)
##         XYFaceAreas = XYFaceAreas * self.dx * self.dy
##         XZFaceAreas = numerix.ones(self.numberOfXZFaces)
##         XZFaceAreas = XZFaceAreas * self.dx * self.dz
##         YZFaceAreas = numerix.ones(self.numberOfYZFaces)
##         YZFaceAreas = YZFaceAreas * self.dy * self.dz
##         self.faceAreas =  numerix.concatenate((XYFaceAreas, XZFaceAreas, YZFaceAreas))

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.

            >>> dx = 0.5
            >>> dy = 2.
            >>> dz = 4.
            >>> nx = 3
            >>> ny = 2
            >>> nz = 1

            >>> mesh = NonUniformGrid3D(nx = nx, ny = ny, nz = nz, dx = dx, dy = dy, dz = dz)

            >>> adjacentCellIDs = (numerix.array([0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 0,
            ...                               1, 2, 3, 4, 5, 0, 0, 1, 2, 3, 3, 4, 5]),
            ...                numerix.array([0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3,
            ...                               4, 5, 3, 4, 5, 0, 1, 2, 2, 3, 4, 5, 5]))
            >>> print(numerix.allequal(mesh._adjacentCellIDs, adjacentCellIDs)) # doctest: +PROCESSOR_0
            True

            >>> vertices = numerix.array(((0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.,
            ...                            0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.),
            ...                           (0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.,
            ...                            0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.),
            ...                           (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            ...                            1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.)))
            >>> vertices *= numerix.array([[dx], [dy], [dz]])

            >>> print(numerix.allequal(vertices, mesh.vertexCoords)) # doctest: +PROCESSOR_0
            True

            >>> faces = numerix.array(((0, 1, 2, 4,  5,  6, 12, 13, 14, 16, 17, 18,  0,  1,  2,  4,  5,  6,  8,  9, 10,  0,  1,  2,  3,  4,  5,  6,  7),
            ...                        (1, 2, 3, 5,  6,  7, 13, 14, 15, 17, 18, 19,  1,  2,  3,  5,  6,  7,  9, 10, 11,  4,  5,  6,  7,  8,  9, 10, 11),
            ...                        (5, 6, 7, 9, 10, 11, 17, 18, 19, 21, 22, 23, 13, 14, 15, 17, 18, 19, 21, 22, 23, 16, 17, 18, 19, 20, 21, 22, 23),
            ...                        (4, 5, 6, 8,  9, 10, 16, 17, 18, 20, 21, 22, 12, 13, 14, 16, 17, 18, 20, 21, 22, 12, 13, 14, 15, 16, 17, 18, 19)))
            >>> print(numerix.allequal(faces, mesh.faceVertexIDs)) # doctest: +PROCESSOR_0
            True

            >>> cells = numerix.array(((21, 22, 23, 25, 26, 27),
            ...                        (22, 23, 24, 26, 27, 28),
            ...                        (12, 13, 14, 15, 16, 17),
            ...                        (15, 16, 17, 18, 19, 20),
            ...                        ( 0,  1,  2,  3,  4,  5),
            ...                        ( 6,  7,  8,  9, 10, 11)))
            >>> print(numerix.allequal(cells, mesh.cellFaceIDs)) # doctest: +PROCESSOR_0
            True

            >>> externalFaces = numerix.array((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 24, 25, 28))
            >>> print(numerix.allequal(externalFaces,
            ...                        numerix.nonzero(mesh.exteriorFaces))) # doctest: +PROCESSOR_0
            True

            >>> internalFaces = numerix.array((15, 16, 17, 22, 23, 26, 27))
            >>> print(numerix.allequal(internalFaces,
            ...                        numerix.nonzero(mesh.interiorFaces))) # doctest: +PROCESSOR_0
            True

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values((( 0,  1,  2,  3,  4,  5,  0,  1,  2,  3,  4,  5,  0,  1,  2,  0, 1, 2,  3,  4,  5,  0, 0, 1,  2,  3, 3, 4,  5),
            ...                                 (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  3, 4, 5, -1, -1, -1, -1, 1, 2, -1, -1, 4, 5, -1)), -1)
            >>> print(numerix.allequal(faceCellIds, mesh.faceCellIDs)) # doctest: +PROCESSOR_0
            True

            >>> xy = dx * dy
            >>> xz = dx * dz
            >>> yz = dy * dz
            >>> faceAreas = numerix.array((xy, xy, xy, xy, xy, xy, xy, xy, xy, xy, xy, xy,
            ...                            xz, xz, xz, xz, xz, xz, xz, xz, xz,
            ...                            yz, yz, yz, yz, yz, yz, yz, yz))
            >>> print(numerix.allclose(faceAreas, mesh._faceAreas, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> faceCoords = numerix.take(vertices, faces, axis=1)
            >>> faceCenters = (faceCoords[..., 0,:] + faceCoords[..., 1,:] + faceCoords[..., 2,:] + faceCoords[..., 3,:]) / 4.
            >>> print(numerix.allclose(faceCenters, mesh.faceCenters, atol = 1e-10, rtol = 1e-10))
            True

            >>> faceNormals = numerix.array((( 0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0, -1, 1, 1, 1, -1, 1, 1, 1),
            ...                              ( 0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, -1, -1, -1, 1, 1, 1, 1, 1, 1,  0, 0, 0, 0,  0, 0, 0, 0),
            ...                              (-1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1,  0,  0,  0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0)))
            >>> print(numerix.allclose(faceNormals, mesh.faceNormals, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> cellToFaceOrientations = numerix.array(((1, -1, -1, 1, -1, -1),
            ...                                         (1, 1, 1, 1, 1, 1),
            ...                                         (1, 1, 1, -1, -1, -1),
            ...                                         (1, 1, 1, 1, 1, 1),
            ...                                         (1, 1, 1, 1, 1, 1),
            ...                                         (1, 1, 1, 1, 1, 1)))
            >>> print(numerix.allequal(cellToFaceOrientations, mesh._cellToFaceOrientations)) # doctest: +PROCESSOR_0
            True

            >>> cellVolumes = numerix.array((dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz))
            >>> print(numerix.allclose(cellVolumes, mesh.cellVolumes, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> cellCenters = numerix.array(((dx/2., 3.*dx/2., 5.*dx/2.,    dx/2., 3.*dx/2., 5.*dx/2.),
            ...                              (dy/2.,    dy/2.,    dy/2., 3.*dy/2., 3.*dy/2., 3.*dy/2.),
            ...                              (dz/2.,    dz/2.,    dz/2.,    dz/2.,    dz/2.,    dz/2.)))
            >>> print(numerix.allclose(cellCenters, mesh.cellCenters, atol = 1e-10, rtol = 1e-10))
            True

            >>> faceToCellDistances = MA.masked_values(((dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dy/2, dx/2, dx/2, dx/2, dx/2, dx/2, dx/2, dx/2, dx/2),
            ...                                         (  -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1, dy/2, dy/2, dy/2,   -1,   -1,   -1,   -1, dx/2, dx/2,   -1,   -1, dx/2,   -1)), -1)
            >>> print(numerix.allclose(faceToCellDistances, mesh._faceToCellDistances, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            1

            >>> cellDistances = numerix.array((dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2,
            ...                                dy/2, dy/2, dy/2, dy, dy, dy, dy/2, dy/2, dy/2,
            ...                                dx/2, dx, dx, dx/2, dx/2, dx, dx, dx/2))
            >>> print(numerix.allclose(cellDistances, mesh._cellDistances, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> print(numerix.allclose(faceToCellDistanceRatios, mesh._faceToCellDistanceRatio, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> areaProjections = faceNormals * faceAreas
            >>> print(numerix.allclose(areaProjections, mesh._areaProjections, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> tangents1 = numerix.array(((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
            ...                            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1),
            ...                            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)))
            >>> print(numerix.allclose(tangents1, mesh._faceTangents1, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> tangents2 = numerix.array(((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            ...                            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            ...                            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)))
            >>> print(numerix.allclose(tangents2, mesh._faceTangents2, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> cellToCellIDs = MA.masked_values(((-1, 0, 1, -1, 3, 4),
            ...                                   (1, 2, -1, 4, 5, -1),
            ...                                   (-1, -1, -1, 0, 1, 2),
            ...                                   (3, 4, 5, -1, -1, -1),
            ...                                   (-1, -1, -1, -1, -1, -1),
            ...                                   (-1, -1, -1, -1, -1, -1)), -1)
            >>> print(numerix.allequal(mesh._cellToCellIDs, cellToCellIDs)) # doctest: +PROCESSOR_0
            True

            >>> cellToCellIDsFilled = numerix.array(((0, 0, 1, 3, 3, 4),
            ...                                      (1, 2, 2, 4, 5, 5),
            ...                                      (0, 1, 2, 0, 1, 2),
            ...                                      (3, 4, 5, 3, 4, 5),
            ...                                      (0, 1, 2, 3, 4, 5),
            ...                                      (0, 1, 2, 3, 4, 5)))
            >>> print(numerix.allequal(mesh._cellToCellIDsFilled, cellToCellIDsFilled)) # doctest: +PROCESSOR_0
            True

            >>> cellToCellDistances = numerix.take(cellDistances, cells)
            >>> print(numerix.allclose(cellToCellDistances, mesh._cellToCellDistances, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> interiorCellIDs = numerix.array(())
            >>> numerix.allequal(interiorCellIDs, mesh._interiorCellIDs)
            True

            >>> exteriorCellIDs = numerix.array((0, 1, 2, 3, 4, 5))
            >>> print(numerix.allequal(exteriorCellIDs, mesh._exteriorCellIDs)) # doctest: +PROCESSOR_0
            True

            >>> cellNormals = numerix.array((((-1, -1, -1, -1, -1, -1),
            ...                               ( 1,  1,  1,  1,  1,  1),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0)),
            ...                              (( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               (-1, -1, -1, -1, -1, -1),
            ...                               ( 1,  1,  1,  1,  1,  1),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0)),
            ...                              (( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               (-1, -1, -1, -1, -1, -1),
            ...                               ( 1,  1,  1,  1,  1,  1))))
            >>> print(numerix.allclose(cellNormals, mesh._cellNormals, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> cellAreaProjections = numerix.array((((-yz, -yz, -yz, -yz, -yz, -yz),
            ...                                       ( yz,  yz,  yz,  yz,  yz,  yz),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0)),
            ...                                      ((  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (-xz, -xz, -xz, -xz, -xz, -xz),
            ...                                       ( xz,  xz,  xz,  xz,  xz,  xz),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0)),
            ...                                      ((  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (  0,   0,   0,   0,   0,   0),
            ...                                       (-xy, -xy, -xy, -xy, -xy, -xy),
            ...                                       ( xy,  xy,  xy,  xy,  xy,  xy))))
            >>> print(numerix.allclose(cellAreaProjections, mesh._cellAreaProjections, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> cellVertexIDs = numerix.array((17, 16, 13, 12, 5, 4, 1, 0))
            >>> cellVertexIDs = numerix.array((cellVertexIDs, cellVertexIDs + 1, cellVertexIDs + 2,
            ...                                cellVertexIDs + 4, cellVertexIDs + 5, cellVertexIDs + 6))
            >>> cellVertexIDs = cellVertexIDs.swapaxes(0, 1)


            >>> print(numerix.allclose(mesh._cellVertexIDs, cellVertexIDs)) # doctest: +PROCESSOR_0
            True

            >>> from fipy.tools import dump
            >>> (f, filename) = dump.write(mesh, extension = '.gz')
            >>> unpickledMesh = dump.read(filename, f)

            >>> print(numerix.allclose(mesh.cellCenters, unpickledMesh.cellCenters))
            True

            The following test was for a bug when `dx`, `dy` or `dz` are arrays.
            The `_calcFaceAreas()` method was commented out to fix this.

            >>> NonUniformGrid3D(nx=2, ny=2, nz=2, dx=(1., 2.), dy=(1., 2.), dz=(1., 2.))
            NonUniformGrid3D(dx=(1.0, 2.0), nx=2, dy=(1.0, 2.0), ny=2, dz=(1.0, 2.0), nz=2)

        Test for https://github.com/usnistgov/fipy/issues/364.

            >>> from fipy.meshes.nonUniformGrid3D import NonUniformGrid3D
            >>> m = NonUniformGrid3D(nx=1, ny=1, nz=9, overlap=1)
            >>> print(min(m.z) == 0.5) # doctest: +SERIAL
            True
            >>> print(min(m.z) == 3.5) # doctest: +PROCESSOR_1_OF_2
            True
            >>> print(min(m.z) == 5.5) # doctest: +PROCESSOR_2_OF_3
            True

        Ensure that ghost faces are excluded from accumulating operations
        (#856).  Four exterior surfaces of :math:`10\times 10\times 10`
        cube mesh should each have a total area of 100, regardless of
        partitioning.

            >>> cube = NonUniformGrid3D(nx=10, dx=1., ny=10, dy=1., nz=10, dz=1.)

            >>> area = (cube._faceAreas * cube.facesBottom).sum()
            >>> print(numerix.allclose(area, 100))
            True

            >>> area = (cube._faceAreas * cube.facesTop).sum()
            >>> print(numerix.allclose(area, 100))
            True

            >>> area = (cube._faceAreas * cube.facesLeft).sum()
            >>> print(numerix.allclose(area, 100))
            True

            >>> area = (cube._faceAreas * cube.facesRight).sum()
            >>> print(numerix.allclose(area, 100))
            True

            >>> area = (cube._faceAreas * cube.facesFront).sum()
            >>> print(numerix.allclose(area, 100))
            True

            >>> area = (cube._faceAreas * cube.facesBack).sum()
            >>> print(numerix.allclose(area, 100))
            True
        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


