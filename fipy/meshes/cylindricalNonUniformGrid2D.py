"""
2D rectangular Mesh
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'


from fipy.tools import numerix

from fipy.meshes.nonUniformGrid2D import NonUniformGrid2D
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import parallelComm

__all__ = ["CylindricalNonUniformGrid2D"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class CylindricalNonUniformGrid2D(NonUniformGrid2D):
    """
    Creates a 2D cylindrical grid mesh with horizontal faces numbered
    first and then vertical faces.
    """
    def __init__(self, dx=1., dy=1., nx=None, ny=None,
                 origin=((0.,), (0.,)), overlap=2, communicator=parallelComm, *args, **kwargs):
        scale = PhysicalField(value=1, unit=PhysicalField(value=dx).unit)
        self.origin = PhysicalField(value=origin)
        self.origin /= scale

        super(CylindricalNonUniformGrid2D, self).__init__(dx=dx, dy=dy, nx=nx, ny=ny, overlap=overlap,
                        communicator=communicator, *args, **kwargs)

        self._faceAreas *= self.faceCenters[0].value

        self._scaledFaceAreas = self._scale['area'] * self._faceAreas
        self._areaProjections = self.faceNormals * self._faceAreas
        self._orientedAreaProjections = self._calcOrientedAreaProjections()
        self._faceAspectRatios = self._calcFaceAspectRatios()

        self._cellAreas = self._calcCellAreas()
        self._cellNormals = self._calcCellNormals()

        self.vertexCoords += self.origin
        self.args['origin'] = self.origin

    def _calcFaceCenters(self):
        return super(CylindricalNonUniformGrid2D, self)._calcFaceCenters() + self.origin

    def _calcCellVolumes(self):
        return super(CylindricalNonUniformGrid2D, self)._calcCellVolumes() \
          * self._calcCellCenters()[0]

    def _translate(self, vector):
        return CylindricalNonUniformGrid2D(dx=self.args['dx'], nx=self.args['nx'],
                                           dy=self.args['dy'], ny=self.args['ny'],
                                           origin=self.args['origin'] + vector,
                                           overlap=self.args['overlap'])

    def __mul__(self, factor):
        if len(numerix.shape(factor)) == 0:
            factor = numerix.resize(factor, (2, 1))

        return CylindricalNonUniformGrid2D(dx=self.args['dx'] * numerix.array(factor[0]), nx=self.args['nx'],
                                           dy=self.args['dy'] * numerix.array(factor[1]), ny=self.args['ny'],
                                           origin=self.args['origin'] * factor,
                                           overlap=self.args['overlap'])

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.

            >>> import fipy as fp

            >>> dx = 0.5
            >>> dy = 2.
            >>> nx = 3
            >>> ny = 2

            >>> mesh = CylindricalNonUniformGrid2D(nx = nx, ny = ny, dx = dx, dy = dy)

            >>> vertices = numerix.array(((0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.),
            ...                           (0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.)))
            >>> vertices *= numerix.array(((dx,), (dy,)))
            >>> print(numerix.allequal(vertices,
            ...                        mesh.vertexCoords)) # doctest: +PROCESSOR_0
            True

            >>> faces = numerix.array(((1, 2, 3, 4, 5, 6, 8, 9, 10, 0, 5, 6, 7, 4, 9, 10, 11),
            ...                        (0, 1, 2, 5, 6, 7, 9, 10, 11, 4, 1, 2, 3, 8, 5, 6, 7)))
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
            >>> faceAreas = faceAreas * mesh.faceCenters[0] 
            >>> print(numerix.allclose(faceAreas, mesh._faceAreas, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True
            >>> ignore = numerix.allclose(faceAreas, mesh._faceAreas).value # doctest: +PROCESSOR_NOT_0

            >>> faceCoords = numerix.take(vertices, faces, axis=1)
            >>> faceCenters = (faceCoords[..., 0,:] + faceCoords[..., 1,:]) / 2.
            >>> print(numerix.allclose(faceCenters, mesh.faceCenters, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True
            >>> ignore = numerix.allclose(faceCenters, mesh.faceCenters).value # doctest: +PROCESSOR_NOT_0

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

            >>> type(mesh.cellCenters)
            <class 'fipy.variables.cellVariable.CellVariable'>

           >>> testCellVolumes = mesh.cellCenters[0].globalValue * numerix.array((dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy))

            >>> print(isinstance(mesh.cellVolumes, numerix.ndarray))
            True

            >>> globalValue = fp.CellVariable(mesh=mesh, value=mesh.cellVolumes).globalValue
            >>> print(numerix.allclose(testCellVolumes, globalValue, atol = 1e-10, rtol = 1e-10))
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
            >>> ignore = numerix.allclose(areaProjections, mesh._areaProjections).value # doctest: +PROCESSOR_NOT_0

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

            >>> interiorCellIDs = numerix.array(())
            >>> print(numerix.allequal(interiorCellIDs, mesh._interiorCellIDs)) # doctest: +PROCESSOR_0
            True

            >>> exteriorCellIDs = numerix.array((0, 1, 2, 3, 4, 5))
            >>> print(numerix.allequal(exteriorCellIDs, mesh._exteriorCellIDs)) # doctest: +PROCESSOR_0
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
            >>> cellAreaProjections[:, 0] = cellAreaProjections[:, 0] * mesh.cellCenters[0] # doctest: +PROCESSOR_0
            >>> cellAreaProjections[:, 1] = cellAreaProjections[:, 1] * (mesh.cellCenters[0] + mesh.dx / 2.) # doctest: +PROCESSOR_0
            >>> cellAreaProjections[:, 2] = cellAreaProjections[:, 2] * mesh.cellCenters[0] # doctest: +PROCESSOR_0
            >>> cellAreaProjections[:, 3] = cellAreaProjections[:, 3] * (mesh.cellCenters[0] - mesh.dx / 2.) # doctest: +PROCESSOR_0
            >>> print(numerix.allclose(cellAreaProjections, mesh._cellAreaProjections, atol = 1e-10, rtol = 1e-10)) # doctest: +PROCESSOR_0
            True

            >>> cellVertexIDs = MA.masked_values(((5, 6, 7, 9, 10, 11),
            ...                                   (4, 5, 6, 8,  9, 10),
            ...                                   (1, 2, 3, 5,  6,  7),
            ...                                   (0, 1, 2, 4,  5,  6)), -1000)

            >>> print(numerix.allclose(mesh._cellVertexIDs, cellVertexIDs)) # doctest: +PROCESSOR_0
            True

            >>> from fipy.tools import dump
            >>> (f, filename) = dump.write(mesh, extension = '.gz')
            >>> unpickledMesh = dump.read(filename, f)

            >>> print(numerix.allclose(mesh.cellCenters, unpickledMesh.cellCenters))
            True

            >>> mesh = CylindricalNonUniformGrid2D(dx=(1., 2.), dy=(1.,)) + ((1.,), (0.,))
            >>> print(mesh.cellCenters)
            [[ 1.5  3. ]
             [ 0.5  0.5]]
            >>> print(fp.CellVariable(mesh=mesh, value=mesh.cellVolumes).globalValue)
            [ 1.5  6. ]

        This test is for https://github.com/usnistgov/fipy/issues/372. Cell
        volumes were being returned as `binOps` rather than arrays.

            >>> m = CylindricalNonUniformGrid2D(dx=(1., 2.), dy=(1., 2.))
            >>> print(isinstance(m.cellVolumes, numerix.ndarray))
            True
            >>> print(isinstance(m._faceAreas, numerix.ndarray))
            True

        If the above types aren't correct, the divergence operator's value can be a `binOp`

            >>> print(isinstance(fp.CellVariable(mesh=m).arithmeticFaceValue.divergence.value, numerix.ndarray))
            True

        Test for https://github.com/usnistgov/fipy/issues/393. `exteriorFaces` were
        `ndarrays` rather than `FaceVariables`.

            >>> print(isinstance(m.facesTop, fp.FaceVariable))
            True


        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


