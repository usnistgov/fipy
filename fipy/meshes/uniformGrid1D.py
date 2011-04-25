#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "uniformGrid1D.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 # ###################################################################
 ##

"""
1D Mesh
"""
__docformat__ = 'restructuredtext'

from fipy.tools.numerix import MA

from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import numerix
from fipy.tools import parallel

from fipy.meshes.builders import UniformGrid1DBuilder
from fipy.meshes.builders import Grid1DBuilder
from fipy.meshes.gridlike import Gridlike1D
from fipy.meshes.uniformGrid import UniformGrid

class UniformGrid1D(UniformGrid):
    """
    Creates a 1D grid mesh.
    
        >>> mesh = UniformGrid1D(nx = 3)
        >>> print mesh.cellCenters
        [[ 0.5  1.5  2.5]]
         
    """
    def __init__(self, dx=1., nx=1, origin=(0,), overlap=2,
                       communicator=parallel):
        builder = UniformGrid1DBuilder()

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

        self.communicator = communicator

        self._setTopology()
 
    def __getstate__(self):
        return Gridlike1D.__getstate__(self)

    def __setstate__(self, dict):
        return Gridlike1D.__setstate__(self, dict)

    def __repr__(self):
        return Gridlike1D.__repr__(self)

    def _isOrthogonal(self):
        return Gridlike1D._isOrthogonal(self)

    @property
    def _concatenatedClass(self):
        return Gridlike1D._concatenatedClass
                                                                
    @property
    def _globalNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Does not include the IDs of boundary cells.

        E.g., would return [0, 1, 4, 5] for mesh A

            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return Gridlike1D._globalNonOverlappingCellIDs(self)

    @property
    def _globalOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 4, 5, 6] for mesh A

            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return Gridlike1D._globalOverlappingCellIDs(self)

    @property
    def _localNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Does not include the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

            A        B
        ------------------
        | 3 | 4 || 4 | 5 |
        ------------------
        | 0 | 1 || 1 | 2 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return Gridlike1D._localNonOverlappingCellIDs(self)

    @property
    def _localOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3, 4, 5] for mesh A

            A        B
        ------------------
        | 3 | 4 || 5 |   |
        ------------------
        | 0 | 1 || 2 |   |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return Gridlike1D._localOverlappingCellIDs(self)

    @property
    def _globalNonOverlappingFaceIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Does not include the IDs of boundary cells.

        E.g., would return [0, 1, 4, 5, 8, 9, 12, 13, 14, 17, 18, 19]
        for mesh A

            A   ||   B
        --8---9---10--11--
       17   18  19  20   21
        --4---5----6---7--
       12   13  14  15   16
        --0---1----2---3--
                ||
                
        .. note:: Trivial except for parallel meshes
        """
        return Gridlike1D._globalNonOverlappingFaceIDs(self)

    @property
    def _globalOverlappingFaceIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 
        14, 15, 17, 18, 19, 20] for mesh A

            A   ||   B
        --8---9---10--11--
       17   18  19  20   21
        --4---5----6---7--
       12   13  14  15   16
        --0---1----2---3--
                ||
                
        .. note:: Trivial except for parallel meshes
        """
        return Gridlike1D._globalOverlappingFaceIDs(self)

    @property
    def _localNonOverlappingFaceIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Does not include the IDs of boundary cells.
        
        E.g., would return [0, 1, 3, 4, 6, 7, 9, 10, 11, 13, 14, 15]
        for mesh A

            A   ||   B
        --6---7-----7---8--
       13   14 15/14 15   16
        --3---4-----4---5--
        9   10 11/10 11   12
        --0---1-----1---2--
                ||
        
        .. note:: Trivial except for parallel meshes
        """
        return Gridlike1D._localNonOverlappingFaceIDs(self)

    @property
    def _localOverlappingFaceIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
        12, 13, 14, 15, 16] for mesh A

            A   ||   B
        --6---7----8------
       13   14  15  16   |
        --3---4----5------
        9   10  11  12   |
        --0---1----2------
                ||
        
        .. note:: Trivial except for parallel meshes
        """
        return Gridlike1D._localOverlappingFaceIDs(self)
     
    """
    Topology set and calc
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
        orientations = numerix.ones((2, self.numberOfCells))
        if self.numberOfCells > 0:
            orientations[0] *= -1
            orientations[0,0] = 1
        return orientations

    @property
    def _adjacentCellIDs(self):
        c1 = numerix.arange(self.numberOfFaces)
        ids = numerix.array((c1 - 1, c1))
        if self.numberOfFaces > 0:
            ids[0,0] = ids[1,0]
            ids[1,-1] = ids[0,-1]
        return ids[0], ids[1]

    @property
    def _cellToCellIDs(self):
        c1 = numerix.arange(self.numberOfCells)
        ids = MA.array((c1 - 1, c1 + 1))
        if self.numberOfCells > 0:
            ids[0,0] = MA.masked
            ids[1,-1] = MA.masked
        return ids
        
    @property
    def _cellToCellIDsFilled(self):
        ids = self._cellToCellIDs.filled()
        if self.numberOfCells > 0:
            ids[0,0] = 0
            ids[1,-1] = self.numberOfCells - 1
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
        return numerix.ones(self.numberOfFaces,'d')

    @property
    def _faceCenters(self):
        return numerix.arange(self.numberOfFaces)[numerix.NewAxis, ...] * self.dx + self.origin

    @property
    def _faceNormals(self):
        faceNormals = numerix.ones((1, self.numberOfFaces), 'd')
        # The left-most face has neighboring cells None and the left-most cell.
        # We must reverse the normal to make fluxes work correctly.
        if self.numberOfFaces > 0:
            faceNormals[...,0] *= -1
        return faceNormals

    @property
    def _orientedFaceNormals(self):
        return self._faceNormals

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
            distances[0,0] = self.dx / 2.
            distances[1,-1] = self.dx / 2.
        return distances

    @property
    def _cellNormals(self):
        normals = numerix.ones((1, 2, self.numberOfCells), 'd')
        if self.numberOfCells > 0:
            normals[:,0] = -1
        return normals
        
    @property
    def _cellAreas(self):
        return numerix.ones((2, self.numberOfCells), 'd')

    @property
    def _cellAreaProjections(self):
        return MA.array(self._cellNormals)

    @property
    def _faceCenters(self):
        return numerix.arange(self.numberOfFaces)[numerix.NewAxis, ...] * self.dx + self.origin
     
    @property
    def _cellVolumes(self):
        return numerix.ones(self.numberOfCells, 'd') * self.dx
                                                          
    """
    Scaled geometry set and calc
    """

    @property
    def _faceToCellDistanceRatio(self):
        distances = numerix.ones(self.numberOfFaces, 'd')
        distances *= 0.5
        if len(distances) > 0:
            distances[0] = 1
            distances[-1] = 1
        return distances

    @property
    def _areaProjections(self):
        return self._faceNormals
        
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
        return UniformGrid1D(dx=self.dx * factor,
                             nx=self.args['nx'],
                             origin=self.args['origin'] * factor,
                             overlap=self.args['overlap'])

    @property
    def _concatenableMesh(self):
        from mesh1D import Mesh1D
        return Mesh1D(vertexCoords = self.vertexCoords, 
                      faceVertexIDs = Grid1DBuilder.createFaces(self.numberOfVertices), 
                      cellFaceIDs = Grid1DBuilder.createCells(self.nx))
                      
    @property
    def _cellFaceIDs(self):
        return MA.array(Grid1DBuilder.createCells(self.nx))

    @property
    def _maxFacesPerCell(self):
        return 2
    
    @property
    def vertexCoords(self):
        return self.faceCenters

    @property
    def faceCellIDs(self):
        c1 = numerix.arange(self.numberOfFaces)
        ids = MA.array((c1 - 1, c1))
        if self.numberOfFaces > 0:
            ids[0,0] = ids[1,0]
            ids[1,0] = MA.masked
            ids[1,-1] = MA.masked
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
           >>> print m._getNearestCellID(([0., .9, 3.],))
           [0 0 2]
           >>> print m._getNearestCellID(([1.1],))
           [1]
           >>> m0 = Grid1D(nx=2, dx=1.)
           >>> m1 = Grid1D(nx=4, dx=.5)
           >>> print m0._getNearestCellID(m1.cellCenters.globalValue)
           [0 0 1 1]
           
        """
        nx = self.globalNumberOfCells
        
        if nx == 0:
            return numerix.arange(0)
            
        x0, = self.cellCenters.globalValue[...,0]        
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

        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
