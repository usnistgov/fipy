from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.meshes.builders.abstractGridBuilder import _AbstractGridBuilder
from fipy.tools import numerix
from fipy.tools import vector
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.meshes.builders.utilityClasses import (_UniformNumPts,
                                                 _DOffsets,
                                                 _UniformOrigin,
                                                 _NonuniformNumPts)

class _Grid3DBuilder(_AbstractGridBuilder):

    def buildGridData(self, *args, **kwargs):
        super(_Grid3DBuilder, self).buildGridData(*args, **kwargs)

        self.numberOfHorizontalRows = self.spatialDict["numHorizontalRows"]
        self.numberOfVerticalColumns = self.spatialDict["numVerticalCols"]
        self.numberOfLayersDeep = self.spatialDict["numLayersDeep"]

    def _calcShape(self):
        return (self.ns[0], self.ns[1], self.ns[2])

    def _calcPhysicalShape(self):
        """Return physical dimensions of `Grid3D`
        """
        from fipy.tools.dimensions.physicalField import PhysicalField
        return PhysicalField(value = (self.ns[0] * self.ds[0] * self.scale,
                                      self.ns[1] * self.ds[1] * self.scale,
                                      self.ns[2] * self.ds[2] * self.scale))

    def _calcMeshSpacing(self):
        return numerix.array((self.ds[0], self.ds[1], self.ds[2]))[..., numerix.newaxis]

    @property
    def _specificGridData(self):
        return [self.numberOfXYFaces,
                self.numberOfXZFaces,
                self.numberOfYZFaces,
                self.numberOfHorizontalRows,
                self.numberOfVerticalColumns,
                self.numberOfLayersDeep]


    @staticmethod
    def createVertices(dx, dy, dz, nx, ny, nz,
                       numVertices, numHorizRows, numVertCols):
        x = _AbstractGridBuilder.calcVertexCoordinates(dx, nx)
        x = numerix.resize(x, (numVertices,))

        y = _AbstractGridBuilder.calcVertexCoordinates(dy, ny)
        y = numerix.repeat(y, numVertCols)
        y = numerix.resize(y, (numVertices,))

        z = _AbstractGridBuilder.calcVertexCoordinates(dz, nz)
        z = numerix.repeat(z, numHorizRows * numVertCols)
        z = numerix.resize(z, (numVertices,))

        return numerix.array((x, y, z))

    @staticmethod
    def createFaces(nx, ny, nz):
        """
        XY faces are first, then XZ faces, then YZ faces
        """
        ## do the XY faces
        v1 = numerix.arange((nx + 1) * (ny))
        v1 = vector.prune(v1, nx + 1, nx)
        v1 = _Grid3DBuilder._repeatWithOffset(v1, (nx + 1) * (ny + 1), nz + 1)
        v2 = v1 + 1
        v3 = v1 + (nx + 2)
        v4 = v1 + (nx + 1)
        XYFaces = numerix.array((v1, v2, v3, v4))

        ## do the XZ faces
        v1 = numerix.arange((nx + 1) * (ny + 1))
        v1 = vector.prune(v1, nx + 1, nx)
        v1 = _Grid3DBuilder._repeatWithOffset(v1, (nx + 1) * (ny + 1), nz)
        v2 = v1 + 1
        v3 = v1 + ((nx + 1)*(ny + 1)) + 1
        v4 = v1 + ((nx + 1)*(ny + 1))
        XZFaces = numerix.array((v1, v2, v3, v4))

        ## do the YZ faces
        v1 = numerix.arange((nx + 1) * ny)
        v1 = _Grid3DBuilder._repeatWithOffset(v1, (nx + 1) * (ny + 1), nz)
        v2 = v1 + (nx + 1)
        v3 = v1 + ((nx + 1)*(ny + 1)) + (nx + 1)
        v4 = v1 + ((nx + 1)*(ny + 1))
        YZFaces = numerix.array((v1, v2, v3, v4))

        numberOfXYFaces = (nx * ny * (nz + 1))
        numberOfXZFaces = (nx * (ny + 1) * nz)
        numberOfYZFaces = ((nx + 1) * ny * nz)
        numberOfFaces = numberOfXYFaces + numberOfXZFaces + numberOfYZFaces

        return ([numberOfXYFaces, numberOfXZFaces, numberOfYZFaces, numberOfFaces],
                numerix.concatenate((XYFaces, XZFaces, YZFaces), axis=1))

    @staticmethod
    def createCells(nx, ny, nz, numXYFaces, numXZFaces, numYZFaces):
        """
        cells = (front face, back face, left face, right face, bottom face, top face)
        front and back faces are YZ faces
        left and right faces are XZ faces
        top and bottom faces are XY faces
        """
        ## front and back faces
        frontFaces = numerix.arange(numYZFaces)
        frontFaces = vector.prune(frontFaces, nx + 1, nx)
        frontFaces = frontFaces + numXYFaces + numXZFaces
        backFaces = frontFaces + 1

        ## left and right faces
        leftFaces = numerix.arange(nx * ny)
        leftFaces = _Grid3DBuilder._repeatWithOffset(leftFaces, nx * (ny + 1), nz)
        leftFaces = numerix.ravel(leftFaces)
        leftFaces = leftFaces + numXYFaces
        rightFaces = leftFaces + nx

        ## bottom and top faces
        bottomFaces = numerix.arange(nx * ny * nz)
        topFaces = bottomFaces + (nx * ny)

        return numerix.array((frontFaces, backFaces, leftFaces,
                              rightFaces, bottomFaces, topFaces))

    @staticmethod
    def _repeatWithOffset(array, offset, reps):
        a = numerix.fromfunction(lambda rnum, x: array + (offset * rnum),
                                 (reps, numerix.size(array))).astype('l')
        return numerix.ravel(a)


    def _packOverlap(self, first, second):
        return {'left': 0, 'right': 0, 'bottom' : 0, 'top' : 0,
                'front': first, 'back': second}

    def _packOffset(self, arg):
        return (0, 0, arg)

class _NonuniformGrid3DBuilder(_Grid3DBuilder):

    def __init__(self):
        self.NumPtsCalcClass = _NonuniformNumPts

        super(_NonuniformGrid3DBuilder, self).__init__()

    def buildGridData(self, *args, **kwargs):
        super(_NonuniformGrid3DBuilder, self).buildGridData(*args, **kwargs)

        ([self.Xoffset, self.Yoffset, self.Zoffset],
         self.ds) = _DOffsets.calcDOffsets(self.ds,
                                           self.ns,
                                           self.offset)

        self.vertices = _Grid3DBuilder.createVertices(self.ds[0], self.ds[1],
                                                      self.ds[2],
                                                      self.ns[0], self.ns[1],
                                                      self.ns[2],
                                                      self.numberOfVertices,
                                                      self.numberOfHorizontalRows,
                                                      self.numberOfVerticalColumns) \
                         + ((self.Xoffset,), (self.Yoffset,), (self.Zoffset,))

        numFacesList, self.faces = _Grid3DBuilder.createFaces(self.ns[0],
                                                              self.ns[1],
                                                              self.ns[2])

        self.numberOfXYFaces = numFacesList[0]
        self.numberOfXZFaces = numFacesList[1]
        self.numberOfYZFaces = numFacesList[2]
        self.numberOfFaces   = numFacesList[3]

        self.cells = _Grid3DBuilder.createCells(self.ns[0],
                                                self.ns[1],
                                                self.ns[2],
                                                self.numberOfXYFaces,
                                                self.numberOfXZFaces,
                                                self.numberOfYZFaces)

    @property
    def _specificGridData(self):
        return super(_NonuniformGrid3DBuilder, self)._specificGridData \
                + [self.vertices,
                   self.faces,
                   self.cells,
                   self.Xoffset, self.Yoffset, self.Zoffset]

class _UniformGrid3DBuilder(_Grid3DBuilder):

    def __init__(self):
        self.NumPtsCalcClass = _UniformNumPts

        super(_UniformGrid3DBuilder, self).__init__()

    def buildGridData(self, ds, ns, overlap, communicator, origin):
        super(_UniformGrid3DBuilder, self).buildGridData(ds, ns, overlap,
                                                        communicator)

        self.origin = _UniformOrigin.calcOrigin(origin,
                                                self.offset, self.ds, self.scale)

        self.numberOfXYFaces = self.ns[0] * self.ns[1] * (self.ns[2] + 1)
        self.numberOfXZFaces = self.ns[0] * (self.ns[1] + 1) * self.ns[2]
        self.numberOfYZFaces = (self.ns[0] + 1) * self.ns[1] * self.ns[2]
        self.numberOfFaces = self.numberOfXYFaces + self.numberOfXZFaces \
                              + self.numberOfYZFaces

    @property
    def _specificGridData(self):
        return super(_UniformGrid3DBuilder, self)._specificGridData \
                + [self.origin]
