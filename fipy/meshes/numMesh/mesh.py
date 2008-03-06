#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
 #                                    created: 11/10/03 {2:44:42 PM} 
 #                                last update: 3/6/08 {9:28:54 AM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Alexander Mont <alexander.mont@nist.gov>
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools.numerix import MA

from fipy.meshes.common.mesh import Mesh as _CommonMesh

from fipy.meshes.meshIterator import FaceIterator
from fipy.meshes.numMesh.cell import Cell

from fipy.tools.dimensions.physicalField import PhysicalField

class MeshAdditionError(Exception):
    pass
    
class Mesh(_CommonMesh):
    """Generic mesh class using numerix to do the calculations

        Meshes contain cells, faces, and vertices.

        This is built for a non-mixed element mesh.
    """

    def __init__(self, vertexCoords, faceVertexIDs, cellFaceIDs):
        """faceVertexIds and cellFacesIds must be padded with minus ones."""

        from fipy.variables.vertexVariable import _VertexVariable
        from fipy.variables.faceVariable import FaceVariable
        from fipy.variables.cellVariable import CellVariable
        
        self.vertexCoords = _VertexVariable(mesh=self, value=vertexCoords, 
                                            _bootstrap=True)
        self.faceVertexIDs = FaceVariable(mesh=self, 
                                          elementshape=faceVertexIDs.shape[:-1],
                                          value=MA.masked_values(faceVertexIDs, -1), 
                                          _bootstrap=True)
        self.cellFaceIDs = CellVariable(mesh=self, 
                                        elementshape=cellFaceIDs.shape[:-1], 
                                        value=MA.masked_values(cellFaceIDs, -1), 
                                        _bootstrap=True)

        _CommonMesh.__init__(self)
        
    """Topology methods"""

    def __add__(self, other):
        if(isinstance(other, Mesh)):
            return self._concatenate(other, smallNumber = 1e-15)
        else:
            return self._translate(other)

    __radd__ = __add__
    
    def __mul__(self, factor):
        newCoords = self.vertexCoords * factor
        newmesh = Mesh(vertexCoords=newCoords, 
                       faceVertexIDs=numerix.array(self.faceVertexIDs), 
                       cellFaceIDs=numerix.array(self.cellFaceIDs))
        return newmesh

    __rmul__ = __mul__

    def _concatenate(self, other, smallNumber):
        return Mesh(**self._getAddedMeshValues(other, smallNumber))

    def _connectFaces(self, faces0, faces1):
        """
        
        Merge faces on the same mesh. This is used to create periodic
        meshes. The first list of faces, `faces1`, will be the faces
        that are used to add to the matrix diagonals. The faces in
        `faces2` will not be used. They aren't deleted but their
        adjacent cells are made to point at `faces1`. The list
        `faces2` are not altered, they still remain as members of
        exterior faces.

           >>> from fipy.meshes.numMesh.grid2D import Grid2D
           >>> mesh = Grid2D(nx = 2, ny = 2, dx = 1., dy = 1.)

           >>> print mesh._getCellFaceIDs()
           [[0 1 2 3]
            [7 8 10 11]
            [2 3 4 5]
            [6 7 9 10]]
            
           >>> mesh._connectFaces(mesh.getFacesLeft(), mesh.getFacesRight())

           >>> print mesh._getCellFaceIDs()
           [[0 1 2 3]
            [7 6 10 9]
            [2 3 4 5]
            [6 7 9 10]]

        """

        ## check for errors

        ## check that faces are members of exterior faces
        from sets import Set
        assert Set(faces0).union(Set(faces1)).issubset(Set(self.getExteriorFaces()))

        ## following assert checks number of faces are equal, normals are opposite and areas are the same
        assert numerix.alltrue(numerix.take(self.areaProjections, faces0, axis=1) 
                               == numerix.take(-self.areaProjections, faces1, axis=1))

        ## extract the adjacent cells for both sets of faces
        self.faceCellIDs = self.faceCellIDs.copy()
        ## set the new adjacent cells for `faces0`
        newFaces0 = self.faceCellIDs[0].take(faces0)
        newFaces1 = self.faceCellIDs[0].take(faces1)
        
        self.faceCellIDs[1].put(faces0, newFaces0)
        self.faceCellIDs[0].put(faces0, newFaces1)
        
        ## extract the face to cell distances for both sets of faces
        self.faceToCellDistances = self.faceToCellDistances.copy()
        ## set the new faceToCellDistances for `faces0`
        newDistances0 = self.faceToCellDistances[0].take(faces0)
        newDistances1 = self.faceToCellDistances[0].take(faces1)
        
        self.faceToCellDistances[1].put(faces0, newDistances0)
        self.faceToCellDistances[0].put(faces0, newDistances1)

        ## calculate new cell distances and add them to faces0
        self.cellDistances = self.cellDistances.copy()
        self.cellDistances.put(faces0, (self.faceToCellDistances[0] 
                                        + self.faceToCellDistances[1]).take(faces0))

        ## change the direction of the face normals for faces0
        self.faceNormals = self.faceNormals.copy()
        for dim in range(self.getDim()):
            faceNormals = self.faceNormals[dim].copy()
            numerix.put(faceNormals, faces0, faceNormals.take(faces1))
            self.faceNormals[dim] = faceNormals

        ## Cells that are adjacent to faces1 are changed to point at faces0
        ## get the cells adjacent to faces1
        faceCellIDs = self.faceCellIDs[0].take(faces1)
        ## get all the adjacent faces for those particular cells
        self.cellFaceIDs = self.cellFaceIDs.copy()
        cellFaceIDs = self.cellFaceIDs.take(faceCellIDs, axis=1).copy()
        
        for i in range(cellFaceIDs.shape[0]):
            ## if the faces is a member of faces1 then change the face to point at
            ## faces0
            facesInFaces1 = (cellFaceIDs[i] == faces1)
            cellFaceIDs[i] = (facesInFaces1 * faces0
                              + ~facesInFaces1 * cellFaceIDs[i])
            ## add those faces back to the main self.cellFaceIDs
            self.cellFaceIDs[i].put(faceCellIDs, cellFaceIDs[i])

        ## calculate new topology
        _CommonMesh._calcTopology(self)

        ## calculate new geometry
        self._calcOrientedFaceNormals()
        self._calcFaceToCellDistanceRatio()
        self._calcCellToCellDistances()
        self._calcScaledGeometry()
        self._calcFaceAspectRatios()
        
    def _getConcatenableMesh(self):
        return self
        
    def _getAddedMeshValues(self, other, smallNumber):
        """
        Returns a `dictionary` with 3 elements: the new mesh `vertexCoords`,
        `faceVertexIDs`, and `cellFaceIDs`.
        """

        other = other._getConcatenableMesh()

        ## check dimensions
        if self.getDim() != other.getDim():
            raise MeshAdditionError, "Dimensions do not match"
        ## compute vertex correlates
        vertexCorrelates = {}
        
        for i in range(other._getNumberOfVertices()):
            diff = self.vertexCoords() - other.vertexCoords()[...,i,numerix.newaxis]
            j = numerix.array(numerix.dot(diff, diff) < smallNumber).nonzero()[0]
            if len(j) == 1:
                vertexCorrelates[i] = int(j)
            elif len(j) > 1:
                raise MeshAdditionError, "Vertices are indistinguishable"

        if vertexCorrelates == {}:
            raise MeshAdditionError, "Vertices are not aligned"

        
         
        ## compute face correlates
        faceCorrelates = {}
        for i in range(other._getNumberOfFaces()):
            currFace = numerix.array(other.faceVertexIDs()[...,i])
            keepGoing = 1
            currIndex = 0 
            for item in currFace:
                if vertexCorrelates.has_key(item):
                    currFace[currIndex] = vertexCorrelates[item]
                    currIndex = currIndex + 1
                else:
                    keepGoing = 0
            if keepGoing == 1:
                for j in range(self._getNumberOfFaces()):
                    if self._equalExceptOrder(currFace, self.faceVertexIDs()[...,j]):
                        faceCorrelates[i] = j
        if faceCorrelates == {}:
            raise MeshAdditionError, "Faces are not aligned"
        
        faceIndicesToAdd = ()
        for i in range(other._getNumberOfFaces()):
            if not faceCorrelates.has_key(i):
                faceIndicesToAdd = faceIndicesToAdd + (i,)
        vertexIndicesToAdd = ()
        for i in range(other._getNumberOfVertices()):
            if not vertexCorrelates.has_key(i):
                vertexIndicesToAdd = vertexIndicesToAdd + (i,)

        ##compute the full face and vertex correlation list
        a = self._getNumberOfFaces()
        for i in faceIndicesToAdd:
            faceCorrelates[i] = a
            a = a + 1
        b = self._getNumberOfVertices()
        for i in vertexIndicesToAdd:
            vertexCorrelates[i] = b
            b = b + 1

        ## compute what the cells are that we need to add
        cellsToAdd = numerix.ones((self.cellFaceIDs.shape[0], other.cellFaceIDs.shape[-1]))
        cellsToAdd = -1 * cellsToAdd

        for j in range(other.cellFaceIDs.shape[-1]):
            for i in range(other.cellFaceIDs.shape[0]):
                cellsToAdd[i, j] = faceCorrelates[int(other.cellFaceIDs[i, j])]

        cellsToAdd = MA.masked_values(cellsToAdd, -1)


        ## compute what the faces are that we need to add
        facesToAdd = numerix.take(other.faceVertexIDs, faceIndicesToAdd, axis=1).getValue()

        for j in range(facesToAdd.shape[-1]):
            for i in range(facesToAdd.shape[0]):
                facesToAdd[i, j] = vertexCorrelates[int(facesToAdd[i, j])]

        ## compute what the vertices are that we need to add
        verticesToAdd = numerix.take(other.vertexCoords, vertexIndicesToAdd, axis=1)

        return {
            'vertexCoords': numerix.concatenate((self.vertexCoords, verticesToAdd), axis=1), 
            'faceVertexIDs': numerix.concatenate((self.faceVertexIDs, facesToAdd), axis=1), 
            'cellFaceIDs': MA.concatenate((self.cellFaceIDs.getValue(), cellsToAdd), axis=1)
            }

    def _equalExceptOrder(self, first, second):
        """Determines if two lists contain the same set of elements, although they may be in different orders. Does not work if one list contains duplicates of an element.
        """
        res = 0
 
        if (len(first) == len(second)):
            res = 1
        for i in first:
            isthisin = 0
            for j in second:
                if (i == j):
                    isthisin = 1
            if(isthisin == 0):
                res = 0
        return res
    

    def _translate(self, vector):
        newCoords = self.vertexCoords + vector
        newmesh = Mesh(vertexCoords=newCoords,
                       faceVertexIDs=self.faceVertexIDs,
                       cellFaceIDs=self.cellFaceIDs)
        return newmesh

    def _calcTopology(self):
        self.dim = self.vertexCoords.shape[0]
        self.numberOfFaces = self.faceVertexIDs.shape[-1]
        self.numberOfCells = self.cellFaceIDs.shape[-1]
        self._calcFaceCellIDs()
        
        _CommonMesh._calcTopology(self)


    """calc Topology methods"""

    def _calcFaceCellIDs(self):
        from fipy.variables.faceVariable import FaceVariable
        class _FaceCellIDsVariable(FaceVariable):
            def __init__(self, mesh):
                FaceVariable.__init__(self, mesh=mesh,
                                      elementshape=(2,))
                self._requires(mesh.cellFaceIDs)
                self._requires(mesh.numberOfFaces)
                
            def _calcValue(self):
                array = MA.array(MA.indices(self.mesh.cellFaceIDs.shape, 'l')[1], 
                                 mask=self.mesh.cellFaceIDs.getMask())
                faceCellIDs = MA.zeros((2, self.mesh.numberOfFaces), 'l')

                numerix.put(faceCellIDs[0], self.mesh.cellFaceIDs[::-1,::-1], array[::-1,::-1])
                numerix.put(faceCellIDs[1], self.mesh.cellFaceIDs, array)
                return MA.sort(MA.array(faceCellIDs,
                                        mask = ((False,) * self.mesh.numberOfFaces, 
                                                (faceCellIDs[0] == faceCellIDs[1]))),
                               axis=0)
                               
        self.faceCellIDs = _FaceCellIDsVariable(mesh=self)

    def _calcVertexFaceIDs(self):
        from fipy.variables.vertexVariable import _VertexVariable
        class _VertexFaceIDsVariable(_VertexVariable):
            def __init__(self, mesh):
                _VertexVariable.__init__(self, mesh=mesh,
                                         elementshape=(2,))
                self._requires(mesh.vertexFaceIDs)
##                 self._requires(mesh.numberOfFaces)
                
            def _calcValue(self):
                array = MA.array(MA.indices(self.mesh.faceVertexIDs.shape, 'l')[1], 
                                 mask=self.mesh.faceVertexIDs.getMask())
                vertexFaceIDs = MA.zeros((2, self.mesh.numberOfVertices), 'l')

                numerix.put(vertexFaceIDs[0], self.mesh.cellFaceIDs[::-1,::-1], array[::-1,::-1])
                numerix.put(faceCellIDs[1], self.mesh.cellFaceIDs, array)
                return MA.sort(MA.array(faceCellIDs,
                                        mask = ((False,) * self.mesh.numberOfFaces, 
                                                (faceCellIDs[0] == faceCellIDs[1]))),
                               axis=0)
                               
        self.vertexFaceIDs = _VertexFaceIDsVariable(mesh=self)


    def _calcInteriorAndExteriorFaceIDs(self):
        self.exteriorFaces = FaceIterator(mesh=self, 
                                          ids=self.faceCellIDs[1].getMask().nonzero())
        self.interiorFaces = FaceIterator(mesh=self, 
                                          ids=(~self.faceCellIDs[1].getMask()).nonzero())

    def _calcInteriorAndExteriorCellIDs(self):
        ids = numerix.take(self.faceCellIDs[0], self.getExteriorFaces(), axis=-1).filled().sorted()
        self.exteriorCellIDs = ids[(ids[:-1] != ids[1:]).append([True] * (len(ids) - len(ids[:-1])))]
        
        from fipy.variables.cellVariable import CellVariable
        self.interiorCellIDs = CellVariable(mesh=self, value=numerix.arange(self.numberOfCells)).delete(self.exteriorCellIDs)
    
##         try:
##             import sets
##             self.exteriorCellIDs = sets.Set(numerix.take(self.faceCellIDs[0], self.getExteriorFaces()))
##             self.interiorCellIDs = list(sets.Set(range(self.numberOfCells)) - self.exteriorCellIDs)
##             self.exteriorCellIDs = list(self.exteriorCellIDs)
##         except:
##             self.exteriorCellIDs = numerix.take(self.faceCellIDs[0], self.getExteriorFaces())
##             tmp = numerix.zeros(self.numberOfCells)
##             numerix.put(tmp, self.exteriorCellIDs, numerix.ones(len(self.exteriorCellIDs)))
##             self.exteriorCellIDs = numerix.nonzero(tmp)            
##             self.interiorCellIDs = numerix.nonzero(numerix.logical_not(tmp))
            
    def _calcCellToFaceOrientations(self):
        tmp = numerix.take(self.faceCellIDs[0], self.cellFaceIDs, axis=-1)
        self.cellToFaceOrientations = (tmp == MA.indices(tmp.shape)[-1]) * 2 - 1

    def _calcAdjacentCellIDs(self):
        mask = self.faceCellIDs[1].getMask()
        self.adjacentCellIDs = (self.faceCellIDs[0].filled(),
                                (mask * self.faceCellIDs[0].filled(0) 
                                + ~mask * self.faceCellIDs[1].filled(0)))


    def _calcCellToCellIDs(self):    
        self.cellToCellIDs = numerix.take(self.faceCellIDs, self.cellFaceIDs, axis=1)
        self.cellToCellIDs = ((self.cellToFaceOrientations == 1) * self.cellToCellIDs[1] 
                              + (self.cellToFaceOrientations != 1) * self.cellToCellIDs[0])
        
    def _calcNumPts(self, d, n = None, axis = "x"):
        """
        Calculate the number of cells along the specified axis, based
        on either the specified number or on the number elements in the
        cell  `d` spacings.
        
        Used by the `Grid` meshes.
        """
        try:
            lend = len(d)
        except TypeError, e:
            return int(n or 1)
            
        n = int(n or lend)
        if n != lend and lend != 1:
            raise IndexError, "n%s != len(d%s)" % (axis, axis)
                
        return n

    def _calcVertexCoordinates(self, d, n):
        """
        Calculate the positions of the vertices along an axis, based on the 
        specified `Cell` `d` spacing or list of `d` spacings.
        
        Used by the `Grid` meshes.
        """
        x = numerix.zeros((n + 1), 'd')
        x[1:] = d
        return numerix.add.accumulate(x)

    """get Topology methods"""

    def getVertexCoords(self):
        return self.vertexCoords

    def getExteriorFaces(self):
        return self.exteriorFaces
            
    def getInteriorFaces(self):
        return self.interiorFaces
        
    def getFaceCellIDs(self):
        return self.faceCellIDs

    def _getFaces(self):
        return FaceIterator(mesh=self,
                            ids=numerix.arange(self.numberOfFaces),
                            checkIDs=False)

    def _getCellsByID(self, ids = None):
        if ids is None:
            ids = range(self.numberOfCells) 
        return [Cell(self, ID) for ID in ids]
    
    def _getMaxFacesPerCell(self):
        return self.cellFaceIDs.shape[0]

    """Geometry methods"""

    def _calcGeometry(self):
        self._calcFaceCenters()
        _CommonMesh._calcGeometry(self)
        self._calcCellNormals()
        
    """calc geometry methods"""

    def _calcFaceAreas(self):
        faceVertexIDs = self.faceVertexIDs.filled(-1)
        substitute = numerix.repeat(faceVertexIDs[numerix.newaxis, 0], 
                                    faceVertexIDs.shape[0], axis=0)
        faceVertexIDs = numerix.where(self.faceVertexIDs.getMaskArray(), substitute, faceVertexIDs)
        faceVertexCoords = numerix.take(self.vertexCoords, faceVertexIDs, axis=1)
        faceOrigins = numerix.repeat(faceVertexCoords[:,0], faceVertexIDs.shape[0], axis=0)
        faceOrigins = numerix.reshape(faceOrigins, MA.shape(faceVertexCoords))
        faceVertexCoords = faceVertexCoords - faceOrigins
        left = range(faceVertexIDs.shape[0])
        right = left[1:] + [left[0]]
        cross = numerix.sum(numerix.crossProd(faceVertexCoords, 
                                              numerix.take(faceVertexCoords, right, axis=1)), 1)
        self.faceAreas = numerix.sqrtDot(cross, cross) / 2.

    def _calcFaceCenters(self):
        faceVertexCoords = numerix.take(self.vertexCoords, self.faceVertexIDs, axis=1)
## 
##         self.faceCenters = MA.filled(MA.average(faceVertexCoords, axis=1))
        
        self.faceCenters = faceVertexCoords.mean(axis=1).filled()

        

    def _calcFaceNormals(self):
        faceVertexIDs = self.faceVertexIDs.filled(0)
        faceVertexCoords = numerix.take(self.vertexCoords, faceVertexIDs, axis=1)
        t1 = faceVertexCoords[:,1,:] - faceVertexCoords[:,0,:]
        t2 = faceVertexCoords[:,2,:] - faceVertexCoords[:,1,:]
        norm = numerix.crossProd(t1, t2)
        ## reordering norm's internal memory for inlining
        norm = norm.copy()
        norm = norm / numerix.sqrtDot(norm, norm)
        
        self.faceNormals = -norm
        
        orientation = 1 - 2 * (numerix.dot(self.faceNormals, self.cellDistanceVectors) < 0)
        self.faceNormals *= orientation

    def _calcFaceCellToCellNormals(self):
        faceCellCentersUp = numerix.take(self.cellCenters, self.getFaceCellIDs()[1], axis=1)
        faceCellCentersDown = numerix.take(self.cellCenters, self.getFaceCellIDs()[0], axis=1)
        faceCellCentersUp = numerix.where(MA.getmaskarray(faceCellCentersUp),
                                          self.getFaceCenters(),
                                          faceCellCentersUp)

        diff = faceCellCentersDown - faceCellCentersUp
        mag = numerix.sqrtDot(diff, diff)
        self.faceCellToCellNormals = diff / mag[numerix.NewAxis, ...]

        orientation = 1 - 2 * (numerix.dot(self.faceNormals, self.faceCellToCellNormals) < 0)
        self.faceCellToCellNormals *= orientation

    def _calcOrientedFaceNormals(self):
        self.orientedFaceNormals = self.faceNormals
        
    def _calcCellVolumes(self):
        tmp = self.faceCenters[0] * self.faceAreas * self.faceNormals[0]
        tmp = numerix.take(tmp, self.cellFaceIDs, axis=-1) * self.cellToFaceOrientations
        self.cellVolumes = tmp.sum(0).filled()
        self.cellVolumes.name = self.__class__.__name__ + ".cellVolumes"

    def _calcCellCenters(self):
        tmp = numerix.take(self.faceCenters, self.cellFaceIDs, axis=1)
        self.cellCenters = tmp.mean(axis=1).filled()
        self.cellCenters.name = self.__class__.__name__ + ".cellCenters"
        
    def _calcFaceToCellDistances(self):
        tmp = self.faceCenters[...,numerix.newaxis,:]
        tmp -= numerix.take(self.cellCenters, self.faceCellIDs, axis=1)
        self.cellToFaceDistanceVectors = tmp
        self.faceToCellDistances = (tmp * tmp).sum(axis=0).sqrt()
        self.faceToCellDistances.name = self.__class__.__name__ + ".faceToCellDistances"

    def _calcCellDistances(self):
        tmp = numerix.take(self.cellCenters, self.faceCellIDs, axis=1)
        tmp = tmp[...,1,:] - tmp[...,0,:]
        self.cellDistanceVectors = (tmp.getMask() * self.cellToFaceDistanceVectors[:,0].filled() 
                                    + ~tmp.getMask() * tmp.filled())
        self.cellDistances = self.cellDistanceVectors.getMag()
        self.cellDistances.name = self.__class__.__name__ + ".cellDistances"

    def _calcFaceToCellDistanceRatio(self):
        dAP = self._getCellDistances()
        dFP = self._getFaceToCellDistances()[0]
        
        self.faceToCellDistanceRatio = (dFP / dAP).filled()
        self.faceToCellDistanceRatio.name = self.__class__.__name__ + ".faceToCellDistanceRatio"

    def _calcAreaProjections(self):
        self.areaProjections = self._getFaceNormals() * self._getFaceAreas()
        self.areaProjections.name = self.__class__.__name__ + ".areaProjections"
        
    def _calcOrientedAreaProjections(self):
        self.orientedAreaProjections = self.areaProjections

    def _calcFaceTangents(self):
        faceVertexCoord = numerix.take(self.vertexCoords, 
                                       self.faceVertexIDs[0], 
                                       axis=1)
        tmp = self.faceCenters - faceVertexCoord
        self.faceTangents1 = tmp / numerix.sqrtDot(tmp, tmp)
        tmp = numerix.crossProd(self.faceTangents1, self.faceNormals)
        self.faceTangents2 = tmp / numerix.sqrtDot(tmp, tmp)
        self.faceTangents1.name = self.__class__.__name__ + ".faceTangents1"
        self.faceTangents2.name = self.__class__.__name__ + ".faceTangents2"
        
    def _calcCellToCellDistances(self):
        self.cellToCellDistances = numerix.take(self.cellDistances, self._getCellFaceIDs(), axis=-1)
        self.cellToCellDistances.name = self.__class__.__name__ + ".cellToCellDistances"

    def _calcCellNormals(self):
        cellNormals = numerix.take(self._getFaceNormals(), self._getCellFaceIDs(), axis=1)
        cellFaceCellIDs = numerix.take(self.faceCellIDs[0], self.cellFaceIDs, axis=-1)
        cellIDs = numerix.repeat(numerix.arange(self.getNumberOfCells())[numerix.newaxis,...], 
                                 self._getMaxFacesPerCell(), 
                                 axis=0)
        direction = (cellFaceCellIDs == cellIDs) * 2 - 1
        self.cellNormals =  direction[numerix.newaxis, ...] * cellNormals
        self.cellNormals.name = self.__class__.__name__ + ".cellNormals"
                         
    """get geometry methods"""

    def getFaceCenters(self):
        return self.faceCenters

    def _getOrderedCellVertexIDs(self):
        return self._getCellVertexIDs()

    def _getCellDistanceNormals(self):
        return self.getCellDistanceVectors() / self.getCellDistances()
        
    def _getCellVertexIDs(self):

        ## Get all the vertices from all the faces for each cell
        cellFaceVertices = numerix.take(self.faceVertexIDs, self.cellFaceIDs, axis=1)

        ## get a sorted list of vertices for each cell 
        cellVertexIDs = cellFaceVertices.reshape((-1, self.getNumberOfCells()))
        cellVertexIDs = cellVertexIDs.sorted(axis=0, fill_value=-1)

        tmp = cellVertexIDs[:-1].masked(cellVertexIDs[:-1] == cellVertexIDs[1:])
        cellVertexIDs = tmp.append(cellVertexIDs[-1, numerix.newaxis], axis=0)
        cellVertexIDs = cellVertexIDs.sorted(axis=0, fill_value=-1)

        ## resize the array to remove extra masked values
        length = cellVertexIDs.getMask().sum(axis=0).min()
        return cellVertexIDs[length:][::-1]

## Below is an ordered version of _getCellVertexIDs()
## It works for the test case in this file (other than the ordering, obviously)
## I've left it in as it may be useful when we need ordered vertices for cells
    
##    def _getOrderedCellVertexIDs(self):

##        ## Get all the vertices from all the faces for each cell
##        from fipy.tools.numerix import take
##        cellFaceVertices = take(self.faceVertexIDs, self.cellFaceIDs)

##        ## get a sorted list of vertices for each cell
##        NCells = self.getNumberOfCells()
##        cellVertexIDs = MA.reshape(cellFaceVertices.flat, (NCells, -1))
##        newmask = MA.getmaskarray(cellVertexIDs).copy()

##        for i in range(len(cellVertexIDs[0]) - 1):
##            for j in range(len(cellVertexIDs[0]))[i + 1:]:

##                newmask[:,j] = MA.where(newmask[:,j],
##                                        newmask[:,j],
##                                        MA.where(MA.filled(cellVertexIDs)[:,i] == MA.filled(cellVertexIDs)[:,j],
##                                                 1,
##                                                 newmask[:,j]))


##        cellVertexIDs = MA.masked_array(cellVertexIDs, newmask)

##        for i in range(len(cellVertexIDs[0]) - 1):
##            j = i + 1
##            while j < len(cellVertexIDs[0]):
##                tmp = cellVertexIDs[:]
##                tmp[:, i] = MA.where(MA.getmaskarray(cellVertexIDs[:,i]),
##                                     MA.where(MA.getmaskarray(cellVertexIDs[:, j]),
##                                              cellVertexIDs[:, i],
##                                              cellVertexIDs[:, j]),
##                                     cellVertexIDs[:, i])
                                                        
##                tmp[:, j] = MA.where(MA.getmaskarray(cellVertexIDs[:,i]),
##                                     MA.where(MA.getmaskarray(cellVertexIDs[:, j]),
##                                              cellVertexIDs[:,j],
##                                              cellVertexIDs[:,i]),
##                                     cellVertexIDs[:, j])

##                cellVertexIDs = tmp[:]

##                j += 1


##        length = len(cellVertexIDs[0]) - min(numerix.sum(MA.getmaskarray(cellVertexIDs), axis = 1))
##        return cellVertexIDs[:, :length]


    """scaling"""

## ##     def setScale(self, value = 1.):
## ##         self.scale = 1.
    
##    def setScale(self):
##	self.scale = PhysicalField(value = 1.)
##        self.faceAreas = self.faceAreas * self.scale**2
##        self.faceCenters = self.faceCenters * self.scale
##        self.cellVolumes = self.cellVolumes * self.scale**3
##        self.cellCenters = self.cellCenters * self.scale
##        self.faceToCellDistances = self.faceToCellDistances * self.scale
##        self.cellDistances = self.cellDistances * self.scale
##        self.areaProjections = self.areaProjections * self.scale**2

## pickling

##    self.__getinitargs__(self):
##        return (self.vertexCoords, self.faceVertexIDs, self.cellFaceIDs)
    

    def __getstate__(self):
        dict = {
            'vertexCoords' : self.vertexCoords.getValue(),            
            'faceVertexIDs' : self.faceVertexIDs.getValue(),
            'cellFaceIDs' : self.cellFaceIDs.getValue() }
        return dict

    def __setstate__(self, dict):
        Mesh.__init__(self, **dict)
     
    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.
        
            >>> dx = 2.
            >>> dy = 1.23456
            >>> dz = 1.e-1
            
            >>> vertices = numerix.array(((0., 1., 1., 0., 0., 1., 1., 0., 2., 2.),
            ...                           (0., 0., 1., 1., 0., 0., 1., 1., 0., 0.),
            ...                           (0., 0., 0., 0., 1., 1., 1., 1., 0., 1.)))
            >>> vertices *= numerix.array([[dx], [dy], [dz]])
            >>> faces = MA.masked_values(((0, 7, 3, 5, 1, 3, 1, 9, 8, 8),
            ...                           (1, 6, 7, 6, 0, 2, 8, 5, 1, 9),
            ...                           (2, 5, 4, 2, 4, 6, 2, 6, 5, 6),
            ...                           (3, 4, 0, 1, 5, 7, -1, -1, 9, 2)), -1)
            >>> cells = MA.masked_values(((0, 3),
            ...                           (1, 6), 
            ...                           (2, 7),
            ...                           (3, 8),
            ...                           (4, 9),
            ...                           (5, -1)), -1)

            >>> mesh = Mesh(vertexCoords=vertices, faceVertexIDs=faces, cellFaceIDs=cells)

            >>> externalFaces = numerix.array((0, 1, 2, 4, 5, 6, 7, 8, 9))
            >>> numerix.allequal(externalFaces, mesh.getExteriorFaces())
            1

            >>> internalFaces = numerix.array((3,))
            >>> numerix.allequal(internalFaces, mesh.getInteriorFaces())
            1

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values((( 0,  0,  0, 0,  0,  0,  1,  1,  1,  1),
            ...                                 (-1, -1, -1, 1, -1, -1, -1, -1, -1, -1)), -1)
            >>> print numerix.allequal(faceCellIds, mesh.getFaceCellIDs())
            1
            
            >>> dxdy = dx * dy
            >>> dxdz = dx * dz
            >>> dydz = dy * dz
            >>> faceAreas = numerix.array((dxdy, dxdy, dydz, dydz, dxdz, dxdz,
            ...                            dxdy/2., dxdy/2., dxdz, numerix.sqrt(dx**2 + dy**2) * dz))
            >>> print numerix.allclose(faceAreas, mesh._getFaceAreas(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = numerix.take(vertices, MA.filled(faces, 0), axis=1)
            >>> faceCenters = faceCoords[...,0,:] + faceCoords[...,1,:] + faceCoords[...,2,:] + faceCoords[...,3,:]
            >>> numVex = numerix.array((4., 4., 4., 4., 4., 4., 3., 3., 4., 4.))
            >>> faceCenters /= numVex
            >>> print numerix.allclose(faceCenters, mesh.getFaceCenters(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> faceNormals = numerix.array((( 0., 0., -1., 1.,  0., 0.,  0., 0.,  0., dy / numerix.sqrt(dy**2 + dx**2)),
            ...                              ( 0., 0.,  0., 0., -1., 1.,  0., 0., -1., dx / numerix.sqrt(dy**2 + dx**2)),
            ...                              (-1., 1.,  0., 0.,  0., 0., -1., 1.,  0., 0.)))
            >>> print numerix.allclose(faceNormals, mesh._getFaceNormals(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = MA.masked_values(((1, -1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, -2)), -2)
            >>> print numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            1
                                             
            >>> cellVolumes = numerix.array((dx*dy*dz, dx*dy*dz / 2.))
            >>> print numerix.allclose(cellVolumes, mesh.getCellVolumes(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> cellCenters = numerix.array(((dx/2., dx+dx/3.),
            ...                              (dy/2.,    dy/3.),
            ...                              (dz/2.,    dz/2.)))
            >>> print numerix.allclose(cellCenters, mesh.getCellCenters(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> d1 = numerix.sqrt((dx / 3.)**2 + (dy / 6.)**2)
            >>> d2 = numerix.sqrt((dx / 6.)**2 + (dy / 3.)**2)
            >>> d3 = numerix.sqrt((dx / 6.)**2 + (dy / 6.)**2)
            >>> d4 = numerix.sqrt((5 * dx / 6.)**2 + (dy / 6.)**2)
            >>> faceToCellDistances = MA.masked_values(((dz / 2., dz / 2., dx / 2., dx / 2., dy / 2., dy / 2., dz / 2., dz / 2., d2, d3),
            ...                                         (     -1,      -1,      -1,      d1,      -1,      -1,      -1,      -1, -1, -1)), -1)
            >>> print numerix.allclose(faceToCellDistances, 
            ...                        mesh._getFaceToCellDistances(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> cellDistances = numerix.array((dz / 2., dz / 2., dx / 2.,
            ...                                d4,
            ...                                dy / 2., dy / 2., dz / 2., dz / 2.,
            ...                                d2,
            ...                                d3))
            >>> print numerix.allclose(cellDistances, mesh._getCellDistances(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> print numerix.allclose(faceToCellDistanceRatios, 
            ...                        mesh._getFaceToCellDistanceRatio(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas
            >>> numerix.allclose(areaProjections, mesh._getAreaProjections(), 
            ...                  atol = 1e-10, rtol = 1e-10)
            1

            >>> v1 = numerix.take(vertices, numerix.array(faces[0]), axis=1)
            >>> tmp = faceCenters - v1
            >>> tangents1 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> print numerix.allclose(tangents1, mesh._getFaceTangents1(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> tmp = numerix.crossProd(tangents1, faceNormals)
            >>> tangents2 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> print numerix.allclose(tangents2, mesh._getFaceTangents2(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((-1,  0),
            ...                                   (-1, -1),
            ...                                   (-1, -1), 
            ...                                   ( 1, -1), 
            ...                                   (-1, -1),
            ...                                   (-1, -1)), -1)
            >>> print numerix.allequal(cellToCellIDs, mesh._getCellToCellIDs())
            1

            >>> cellToCellDistances = MA.masked_values(((dz / 2., d4),
            ...                                         (dz / 2., -1),
            ...                                         (dx / 2., -1),
            ...                                         (     d4, -1),
            ...                                         (dy / 2., -1),
            ...                                         (dy / 2., -1)), -1)
            >>> print numerix.allclose(cellToCellDistances, 
            ...                        mesh._getCellToCellDistances(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> interiorCellIDs = numerix.array(())
            >>> print numerix.allequal(interiorCellIDs, mesh._getInteriorCellIDs())
            1

            >>> exteriorCellIDs = numerix.array((0, 1))
            >>> print numerix.allequal(exteriorCellIDs, mesh._getExteriorCellIDs())
            1

            >>> cellNormals = MA.masked_values((((0, -1),
            ...                                  (0, 0),
            ...                                  (-1, 0),
            ...                                  (1, 0),
            ...                                  (0, dy / numerix.sqrt(dx**2 + dy**2)),
            ...                                  (0, -1000)),
            ...                                 ((0, 0),
            ...                                  (0, 0),
            ...                                  (0, 0),
            ...                                  (0, -1),
            ...                                  (-1, dx / numerix.sqrt(dx**2 + dy**2)),
            ...                                  (1, -1000)),
            ...                                 ((-1, 0),
            ...                                  (1, -1),
            ...                                  (0, 1),
            ...                                  (0, 0),
            ...                                  (0, 0),
            ...                                  (0, -1000))), -1000)
            >>> print numerix.allclose(cellNormals, mesh._getCellNormals(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> cellAreaProjections = MA.masked_values((((0, -dy * dz),
            ...                                          (0, 0),
            ...                                          (-dy * dz, 0),
            ...                                          ( dy * dz, 0),
            ...                                          (0, dy * dz),
            ...                                          (0, -1000)),
            ...                                         ((0, 0),
            ...                                          (0, 0),
            ...                                          (0, 0),
            ...                                          (0, -dx * dz),
            ...                                          (-dx * dz, dx * dz),
            ...                                          ( dx * dz, -1000)),
            ...                                         ((-dx * dy, 0),
            ...                                          ( dx * dy, -dx * dy / 2.),
            ...                                          (0, dx * dy / 2.),
            ...                                          (0, 0),
            ...                                          (0, 0),
            ...                                          (0, -1000))), -1000)
            >>> print numerix.allclose(cellAreaProjections, 
            ...                        mesh._getCellAreaProjections(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> cellVertexIDs = MA.masked_values(((7, 6, 5, 4, 3, 2, 1, 0), 
            ...                                   (9, 8, 6, 5, 2, 1, -1000, -1000)), 
            ...                                  -1000)
            >>> cellVertexIDs = MA.masked_values(((7, 9),
            ...                                   (6, 8),
            ...                                   (5, 6),
            ...                                   (4, 5),
            ...                                   (3, 2),
            ...                                   (2, 1),
            ...                                   (1, -1000),
            ...                                   (0, -1000)), -1000)
            >>> print numerix.allclose(cellVertexIDs, mesh._getCellVertexIDs())
            1


            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')
            >>> unpickledMesh = dump.read(filename, f)

            >>> print numerix.allequal(mesh.getCellCenters(), unpickledMesh.getCellCenters())
            1

            >>> dx = 1.
            >>> dy = 1.
            >>> nx = 10
            >>> ny = 2
            >>> from fipy.meshes.grid2D import Grid2D
            >>> gridMesh = Grid2D(dx, dy, nx, ny)
            >>> from fipy.meshes.tri2D import Tri2D
            >>> triMesh = Tri2D(dx, dy, nx, 1) + [[dx*nx], [0]]
            >>> bigMesh = gridMesh + triMesh
            >>> volumes = numerix.ones(bigMesh.getNumberOfCells(), 'd')
            >>> volumes[20:] = 0.25
            >>> print numerix.allclose(bigMesh.getCellVolumes(), volumes)
            1
            
        """

                      
### test test test    

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
