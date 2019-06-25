from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.meshes.builders.abstractGridBuilder import _AbstractGridBuilder
from fipy.tools import inline
from fipy.tools import numerix
from fipy.tools import vector
from fipy.meshes.builders.utilityClasses import (_UniformNumPts,
                                                 _DOffsets,
                                                 _UniformOrigin,
                                                 _NonuniformNumPts)

class _Grid2DBuilder(_AbstractGridBuilder):

    def buildGridData(self, *args, **kwargs):
        # call super for side-effects
        super(_Grid2DBuilder, self).buildGridData(*args, **kwargs)

        self.numberOfVerticalColumns = self.spatialDict["numVerticalCols"]
        self.numberOfHorizontalRows = self.spatialDict["numHorizontalRows"]

    def _dsUniformLen(self):
        """
        Return True if all entries in `self.ds` are the same length, False
        otherwise.

        Exists to get around the fact that `_calcPhysicalShape` and
        `_calcMeshSpacing` don't work for cylindrical grids.
        """

        # if the first entry in `ds` is non-scalar
        if numerix.shape(self.ds[0]) != ():
            lenDs = len(self.ds[0])

            for d in self.ds[1:]:
                if numerix.shape(d) == () or len(d) != lenDs:
                    return False

        # if any other entry in `ds` is non-scalar and first isn't
        elif True in [numerix.shape(d) != () for d in self.ds[1:]]:
            return False

        return True

    def _calcShape(self):
        return (self.ns[0], self.ns[1])

    def _calcPhysicalShape(self):
        """Return physical dimensions of `Grid2D`
        """
        from fipy.tools.dimensions.physicalField import PhysicalField

        if self._dsUniformLen():
            return PhysicalField(value = (self.ns[0] * self.ds[0] * self.scale,
                                          self.ns[1] * self.ds[1] * self.scale))
        else:
            return None

    def _calcMeshSpacing(self):
        if self._dsUniformLen():
            return numerix.array((self.ds[0], self.ds[1]))[..., numerix.newaxis]
        else:
            return None

    @property
    def _specificGridData(self):
        return [self.numberOfHorizontalRows,
                self.numberOfVerticalColumns,
                self.numberOfHorizontalFaces]

    @staticmethod
    def createVertices(nx, ny, dx, dy, numVerts, numVertCols):
        x = _AbstractGridBuilder.calcVertexCoordinates(dx, nx)
        x = numerix.resize(x, (numVerts,))

        y = _AbstractGridBuilder.calcVertexCoordinates(dy, ny)
        y = numerix.repeat(y, numVertCols)

        return numerix.array((x, y))

    @staticmethod
    def createFaces(nx, numVerts, numVertCols):
        """
        `v1`, `v2` refer to the vertices.
        Horizontal faces are first

        Ugly return to avoid side-effects.
        """
        v1 = numerix.arange(numVerts)
        v2 = v1 + 1

        horizontalFaces = vector.prune(numerix.array((v1, v2)), numVertCols, nx, axis=1)

        v1 = numerix.arange(numVerts - numVertCols)
        v2 = v1 + numVertCols
        verticalFaces =  numerix.array((v1, v2))

        ## The cell normals must point out of the cell.
        ## The left and bottom faces have only one neighboring cell,
        ## in the 2nd neighbor position (there is nothing in the 1st).
        ##
        ## reverse some of the face orientations to obtain the correct normals

        tmp = horizontalFaces.copy()
        horizontalFaces[0, :nx] = tmp[1, :nx]
        horizontalFaces[1, :nx] = tmp[0, :nx]

        numberOfHorizontalFaces = horizontalFaces.shape[-1]

        tmp = verticalFaces.copy()
        verticalFaces[0,:] = tmp[1,:]
        verticalFaces[1,:] = tmp[0,:]
        if numVertCols > 0:
            verticalFaces[0, ::numVertCols] = tmp[0, ::numVertCols]
            verticalFaces[1, ::numVertCols] = tmp[1, ::numVertCols]

        return (numerix.concatenate((horizontalFaces, verticalFaces), axis=1),
                numberOfHorizontalFaces)

    if inline.doInline:
        @staticmethod
        def createCells(nx, ny, numFaces, numHorizFaces, numVertCols):
            """
            `cells = (f1, f2, f3, f4)` going anticlockwise.
            `f1` etc. refer to the faces
            """
            cellFaceIDs = numerix.zeros((4, nx * ny), 'l')

            inline._runInline("""
                int ID = j * ni + i;
                int NCELLS = ni * nj;
                cellFaceIDs[ID + 0 * NCELLS] = ID;
                cellFaceIDs[ID + 2 * NCELLS] = cellFaceIDs[ID + 0 * NCELLS] + ni;
                cellFaceIDs[ID + 3 * NCELLS] = horizontalFaces + ID + j;
                cellFaceIDs[ID + 1 * NCELLS] = cellFaceIDs[ID + 3 * NCELLS] + 1;
            """,
            horizontalFaces=numHorizFaces,
            cellFaceIDs=cellFaceIDs,
            ni=nx,
            nj=ny)

            return cellFaceIDs

    else:
        @staticmethod
        def createCells(nx, ny, numFaces, numHorizFaces, numVertCols):
            """
            `cells = (f1, f2, f3, f4)` going anticlockwise.
            `f1` etc. refer to the faces
            """
            cellFaceIDs = numerix.zeros((4, nx * ny), 'l')
            faceIDs = numerix.arange(numFaces)
            if numFaces > 0:
                cellFaceIDs[0,:] = faceIDs[:numHorizFaces - nx]
                cellFaceIDs[2,:] = cellFaceIDs[0,:] + nx
                cellFaceIDs[1,:] = vector.prune(faceIDs[numHorizFaces:],
                                                numVertCols)
                cellFaceIDs[3,:] = cellFaceIDs[1,:] - 1
            return cellFaceIDs

    def _packOverlap(self, first, second):
        return {'left': 0, 'right': 0, 'bottom': first, 'top': second}

    def _packOffset(self, arg):
        return (0, arg)

class _NonuniformGrid2DBuilder(_Grid2DBuilder):

    def __init__(self):
        self.NumPtsCalcClass = _NonuniformNumPts

        super(_NonuniformGrid2DBuilder, self).__init__()

    def buildGridData(self, *args, **kwargs):
        # call super for side-effects
        super(_NonuniformGrid2DBuilder, self).buildGridData(*args, **kwargs)

        (self.offsets,
         self.ds) = _DOffsets.calcDOffsets(self.ds, self.ns, self.offset)

        self.vertices = _Grid2DBuilder.createVertices(self.ns[0], self.ns[1],
                                        self.ds[0], self.ds[1],
                                        self.numberOfVertices,
                                        self.numberOfVerticalColumns) \
                          + ((self.offsets[0],), (self.offsets[1],))

        (self.faces,
         self.numberOfHorizontalFaces) = _Grid2DBuilder.createFaces(self.ns[0],
                                          self.numberOfVertices,
                                          self.numberOfVerticalColumns)
        self.numberOfFaces = len(self.faces[0])
        self.cells = _Grid2DBuilder.createCells(self.ns[0], self.ns[1],
                                               self.numberOfFaces,
                                               self.numberOfHorizontalFaces,
                                               self.numberOfVerticalColumns)

    @property
    def _specificGridData(self):
        return super(_NonuniformGrid2DBuilder, self)._specificGridData \
                 + [self.vertices,
                    self.faces,
                    self.cells,
                    self.offsets]

class _UniformGrid2DBuilder(_Grid2DBuilder):

    def __init__(self):
        self.NumPtsCalcClass = _UniformNumPts

        super(_UniformGrid2DBuilder, self).__init__()

    def buildGridData(self, ds, ns, overlap, communicator, origin):
        # call super for side-effects
        super(_UniformGrid2DBuilder, self).buildGridData(ds, ns, overlap,
                                                        communicator)

        self.origin = _UniformOrigin.calcOrigin(origin,
                                                self.offset, self.ds, self.scale)

        self.numberOfHorizontalFaces = self.ns[0] * self.numberOfHorizontalRows
        self.numberOfVerticalFaces = self.numberOfVerticalColumns * self.ns[1]
        self.numberOfFaces = self.numberOfHorizontalFaces \
                               + self.numberOfVerticalFaces

    @property
    def _specificGridData(self):
        return super(_UniformGrid2DBuilder, self)._specificGridData \
                + [self.numberOfVerticalFaces,
                   self.origin]
