"""
1D Mesh
"""
from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools import parallelComm

from fipy.meshes.uniformGrid import UniformGrid
from fipy.meshes.builders import _UniformGrid1DBuilder
from fipy.meshes.builders import _Grid1DBuilder
from fipy.meshes.representations.gridRepresentation import _Grid1DRepresentation
from fipy.meshes.topologies.gridTopology import _Grid1DTopology

__all__ = ["UniformGrid1D"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class UniformGrid1D(UniformGrid):
    """
    Creates a 1D grid mesh.

        >>> mesh = UniformGrid1D(nx = 3)
        >>> print(mesh.cellCenters)
        [[ 0.5  1.5  2.5]]

    """
    def __init__(self, dx=1., nx=1, origin=(0,), overlap=2,
                 communicator=parallelComm,
                 _RepresentationClass=_Grid1DRepresentation,
                 _TopologyClass=_Grid1DTopology):

        super(UniformGrid1D, self).__init__(communicator=communicator,
                                            _RepresentationClass=_RepresentationClass,
                                            _TopologyClass=_TopologyClass)

        builder = _UniformGrid1DBuilder()

        origin = numerix.array(origin)

        self.args = {
            'dx': dx,
            'nx': nx,
            'origin': origin,
            'overlap': overlap
        }

        self._scale = {
            'length': 1.,
            'area': 1.,
            'volume': 1.
        }

        builder.buildGridData([dx], [nx], overlap, communicator, origin)

        ([self.dx],
         [self.nx],
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
         self.occupiedNodes,
         self.origin) = builder.gridData

        self._setTopology()

    """
    Topology set and calculate
    """

    def _setTopology(self):
        self._exteriorFaces = self.facesLeft | self.facesRight

    @property
    def _interiorFaces(self):
        from fipy.variables.faceVariable import FaceVariable
        interiorFaces = FaceVariable(mesh=self, value=False)
        interiorFaces[numerix.arange(self.numberOfFaces-2) + 1] = True
        return interiorFaces

    @property
    def _cellToFaceOrientations(self):
        orientations = numerix.ones((2, self.numberOfCells), 'l')
        if self.numberOfCells > 0:
            orientations[0] *= -1
            orientations[0, 0] = 1
        return orientations

    @property
    def _adjacentCellIDs(self):
        c1 = numerix.arange(self.numberOfFaces)
        ids = numerix.array((c1 - 1, c1))
        if self.numberOfFaces > 0:
            ids[0, 0] = ids[1, 0]
            ids[1, -1] = ids[0, -1]
        return ids[0], ids[1]

    @property
    def _cellToCellIDs(self):
        c1 = numerix.arange(self.numberOfCells)
        ids = MA.array((c1 - 1, c1 + 1))
        if self.numberOfCells > 0:
            ids[0, 0] = MA.masked
            ids[1, -1] = MA.masked
        return ids

    @property
    def _cellToCellIDsFilled(self):
        ids = self._cellToCellIDs.filled()
        if self.numberOfCells > 0:
            ids[0, 0] = 0
            ids[1, -1] = self.numberOfCells - 1
        return ids

    def _getExteriorFaces(self):
        return self._exteriorFaces

    def _setExteriorFaces(self, e):
        self._exteriorFaces = e

    exteriorFaces = property(_getExteriorFaces, _setExteriorFaces)

    """
    Geometry set and calc
    """

    @property
    def _faceAreas(self):
        return numerix.ones(self.numberOfFaces, 'd')

    @property
    def _faceCenters(self):
        return numerix.arange(self.numberOfFaces)[numerix.NewAxis, ...] * self.dx + self.origin

    @property
    def faceNormals(self):
        faceNormals = numerix.ones((1, self.numberOfFaces), 'd')
        # The left-most face has neighboring cells None and the left-most cell.
        # We must reverse the normal to make fluxes work correctly.
        if self.numberOfFaces > 0:
            faceNormals[..., 0] *= -1
        return faceNormals

    @property
    def _orientedFaceNormals(self):
        return self.faceNormals

    @property
    def _cellVolumes(self):
        return numerix.ones(self.numberOfCells, 'd') * self.dx

    @property
    def _cellCenters(self):
        ccs = ((numerix.arange(self.numberOfCells)[numerix.NewAxis, ...] + 0.5) \
               * self.dx + self.origin) * self.scale['length']
        return ccs

    @property
    def _cellDistances(self):
        distances = numerix.ones(self.numberOfFaces, 'd')
        distances *= self.dx
        if len(distances) > 0:
            distances[0] = self.dx / 2.
            distances[-1] = self.dx / 2.
        return distances

    @property
    def _faceTangents1(self):
        return numerix.zeros(self.numberOfFaces, 'd')[numerix.NewAxis, ...]

    @property
    def _faceTangents2(self):
        return numerix.zeros(self.numberOfFaces, 'd')[numerix.NewAxis, ...]

    @property
    def _cellToCellDistances(self):
        distances = MA.zeros((2, self.numberOfCells), 'd')
        distances[:] = self.dx
        if self.numberOfCells > 0:
            distances[0, 0] = self.dx / 2.
            distances[1, -1] = self.dx / 2.
        return distances

    @property
    def _cellNormals(self):
        normals = numerix.ones((1, 2, self.numberOfCells), 'd')
        if self.numberOfCells > 0:
            normals[:, 0] = -1
        return normals

    @property
    def _cellAreas(self):
        return numerix.ones((2, self.numberOfCells), 'd')

    @property
    def _cellAreaProjections(self):
        return MA.array(self._cellNormals)

    @property
    def _cellVolumes(self):
        return numerix.ones(self.numberOfCells, 'd') * self.dx

    """
    Scaled geometry set and calculate
    """

    @property
    def _faceToCellDistanceRatio(self):
        """how far face is from first to second cell
        
        distance from center of face to center of first cell divided by distance
        between cell centers
        """
        distances = numerix.ones(self.numberOfFaces, 'd')
        distances *= 0.5
        if len(distances) > 0:
            distances[0] = 1
            distances[-1] = 1
        return distances

    @property
    def _areaProjections(self):
        return self.faceNormals

    @property
    def _orientedAreaProjections(self):
        return self._areaProjections

    @property
    def _getFaceAspectRatios(self):
        return 1. / self._cellDistances


    def _translate(self, vector):
        return UniformGrid1D(dx=self.dx,
                             nx=self.args['nx'],
                             origin=self.args['origin'] + numerix.array(vector),
                             overlap=self.args['overlap'])

    def __mul__(self, factor):
        dx = numerix.array([[self.args['dx']]])
        dx *= factor

        return UniformGrid1D(dx=dx[0,0],
                             nx=self.args['nx'],
                             origin=self.args['origin'] * factor,
                             overlap=self.args['overlap'])

    @property
    def _concatenableMesh(self):
        from fipy.meshes.mesh1D import Mesh1D
        return Mesh1D(vertexCoords = self.vertexCoords,
                      faceVertexIDs = _Grid1DBuilder.createFaces(self.numberOfVertices),
                      cellFaceIDs = _Grid1DBuilder.createCells(self.nx))

    @property
    def _cellFaceIDs(self):
        return MA.array(_Grid1DBuilder.createCells(self.nx))

    @property
    def _maxFacesPerCell(self):
        return 2

    @property
    def vertexCoords(self):
        return numerix.array(self.faceCenters)

    @property
    def faceCellIDs(self):
        c1 = numerix.arange(self.numberOfFaces)
        ids = MA.array((c1 - 1, c1))
        if self.numberOfFaces > 0:
            ids[0, 0] = ids[1, 0]
            ids[1, 0] = MA.masked
            ids[1, -1] = MA.masked
        return ids

    @property
    def _cellVertexIDs(self):
        c1 = numerix.arange(self.numberOfCells)
        return numerix.array((c1 + 1, c1))

    def _getNearestCellID(self, points):
        """
        Test cases

           >>> from fipy import *
           >>> m = Grid1D(nx=3)
           >>> print(m._getNearestCellID(([0., .9, 3.],)))
           [0 0 2]
           >>> print(m._getNearestCellID(([1.1],)))
           [1]
           >>> m0 = Grid1D(nx=2, dx=1.)
           >>> m1 = Grid1D(nx=4, dx=.5)
           >>> print(m0._getNearestCellID(m1.cellCenters.globalValue))
           [0 0 1 1]

        """
        nx = self.globalNumberOfCells

        if nx == 0:
            return numerix.arange(0)

        x0, = self.cellCenters.globalValue[..., 0]
        xi, = points
        dx = self.dx

        i = numerix.array(numerix.rint(((xi - x0) / dx)), 'l')
        i[i < 0] = 0
        i[i > nx - 1] = nx - 1

        return i

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected. The following was broken, now fixed.

            >>> from fipy import *
            >>> mesh = Grid1D(nx=3., dx=1.)
            >>> var = CellVariable(mesh=mesh)
            >>> DiffusionTerm().solve(var, solver=DummySolver())

        Size of global value should not depend on number of processors (#400)

            >>> mesh = Grid1D(nx=10)
            >>> print(mesh.cellCenters.globalValue.shape)
            (1, 10)
            >>> print(mesh.faceCenters.globalValue.shape)
            (1, 11)
        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()

