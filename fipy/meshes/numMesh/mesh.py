#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
 #                                    created: 11/10/03 {2:44:42 PM} 
 #                                last update: 10/22/04 {4:11:54 PM} 
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

"""Generic mesh class using Numeric to do the calculations

    Meshes contain cells, faces, and vertices.

    This is built for a non-mixed element mesh.

    Test cases:
    
   >>> from fipy.meshes.grid2D import Grid2D
   >>> from fipy.meshes.numMesh.grid3D import Grid3D
   >>> from fipy.meshes.numMesh.tri2D import Tri2D
   >>> basemesh = Grid2D(dx = 1.0, dy = 1.0, nx = 2, ny = 2)
   >>> dilatedMesh = basemesh * (3, 2)
   >>> print dilatedMesh.getVertexCoords().tolist()
   [[0.0, 0.0], [3.0, 0.0], [6.0, 0.0], [0.0, 2.0], [3.0, 2.0], [6.0, 2.0], [0.0, 4.0], [3.0, 4.0], [6.0, 4.0]]

   >>> translatedMesh = basemesh + (5, 10)
   >>> print translatedMesh.getVertexCoords().tolist()
   [[5.0, 10.0], [6.0, 10.0], [7.0, 10.0], [5.0, 11.0], [6.0, 11.0], [7.0, 11.0], [5.0, 12.0], [6.0, 12.0], [7.0, 12.0]]

   >>> addedMesh = basemesh + (basemesh + (2, 0))
   >>> print addedMesh.getVertexCoords().tolist()
   [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 1.0], [1.0, 1.0], [2.0, 1.0], [0.0, 2.0], [1.0, 2.0], [2.0, 2.0], [3.0, 0.0], [4.0, 0.0], [3.0, 1.0], [4.0, 1.0], [3.0, 2.0], [4.0, 2.0]]

   >>> print addedMesh.getCellFaceIDs()
   [[ 0, 7, 2, 6,]
    [ 1, 8, 3, 7,]
    [ 2,10, 4, 9,]
    [ 3,11, 5,10,]
    [12,18,14, 8,]
    [13,19,15,18,]
    [14,20,16,11,]
    [15,21,17,20,]]
    
   >>> print addedMesh.faceVertexIDs
   [[ 1, 0,]
    [ 2, 1,]
    [ 3, 4,]
    [ 4, 5,]
    [ 6, 7,]
    [ 7, 8,]
    [ 0, 3,]
    [ 4, 1,]
    [ 5, 2,]
    [ 3, 6,]
    [ 7, 4,]
    [ 8, 5,]
    [ 9, 2,]
    [10, 9,]
    [ 5,11,]
    [11,12,]
    [ 8,13,]
    [13,14,]
    [11, 9,]
    [12,10,]
    [13,11,]
    [14,12,]]
   
   >>> addedMesh = basemesh + (basemesh + (3, 0))
   Traceback (most recent call last):
   ...
   MeshAdditionError: Vertices are not aligned

   >>> addedMesh = basemesh + (basemesh + (2, 2))
   Traceback (most recent call last):
   ...
   MeshAdditionError: Faces are not aligned

   >>> triMesh = Tri2D(dx = 1.0, dy = 1.0, nx = 2, ny = 1)
   >>> triMesh = triMesh + (2, 0)
   >>> triAddedMesh = basemesh + triMesh
   >>> print triAddedMesh.getVertexCoords().tolist()
   [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 1.0], [1.0, 1.0], [2.0, 1.0], [0.0, 2.0], [1.0, 2.0], [2.0, 2.0], [3.0, 0.0], [4.0, 0.0], [3.0, 1.0], [4.0, 1.0], [2.5, 0.5], [3.5, 0.5]]

   >>> print triAddedMesh.getCellFaceIDs()
   [[0 ,7 ,2 ,6 ,]
    [1 ,8 ,3 ,7 ,]
    [2 ,10 ,4 ,9 ,]
    [3 ,11 ,5 ,10 ,]
    [16 ,20 ,24 ,-- ,]
    [17 ,21 ,25 ,-- ,]
    [14 ,22 ,24 ,-- ,]
    [15 ,23 ,25 ,-- ,]
    [8 ,18 ,22 ,-- ,]
    [16 ,19 ,23 ,-- ,]
    [12 ,18 ,20 ,-- ,]
    [13 ,19 ,21 ,-- ,]]

   >>> print triAddedMesh.faceVertexIDs
   [[ 1, 0,]
    [ 2, 1,]
    [ 3, 4,]
    [ 4, 5,]
    [ 6, 7,]
    [ 7, 8,]
    [ 0, 3,]
    [ 4, 1,]
    [ 5, 2,]
    [ 3, 6,]
    [ 7, 4,]
    [ 8, 5,]
    [ 9, 2,]
    [10, 9,]
    [ 5,11,]
    [11,12,]
    [11, 9,]
    [12,10,]
    [13, 2,]
    [14, 9,]
    [ 9,13,]
    [10,14,]
    [13, 5,]
    [14,11,]
    [13,11,]
    [14,12,]]

   >>> ThreeDBaseMesh = Grid3D(dx = 1.0, dy = 1.0, dz = 1.0, nx = 2, ny = 2, nz = 2)
   >>> ThreeDSecondMesh = Grid3D(dx = 1.0, dy = 1.0, dz = 1.0, nx = 1, ny = 1, nz = 1)
   >>> ThreeDAddedMesh = ThreeDBaseMesh + (ThreeDSecondMesh + (2, 0, 0))
   >>> print ThreeDAddedMesh.getVertexCoords().tolist()
   [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0], [0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [2.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [0.0, 2.0, 1.0], [1.0, 2.0, 1.0], [2.0, 2.0, 1.0], [0.0, 0.0, 2.0], [1.0, 0.0, 2.0], [2.0, 0.0, 2.0], [0.0, 1.0, 2.0], [1.0, 1.0, 2.0], [2.0, 1.0, 2.0], [0.0, 2.0, 2.0], [1.0, 2.0, 2.0], [2.0, 2.0, 2.0], [3.0, 0.0, 0.0], [3.0, 1.0, 0.0], [3.0, 0.0, 1.0], [3.0, 1.0, 1.0]]

   >>> print ThreeDAddedMesh.getCellFaceIDs()
   [[24,25,12,14, 0, 4,]
    [25,26,13,15, 1, 5,]
    [27,28,14,16, 2, 6,]
    [28,29,15,17, 3, 7,]
    [30,31,18,20, 4, 8,]
    [31,32,19,21, 5, 9,]
    [33,34,20,22, 6,10,]
    [34,35,21,23, 7,11,]
    [26,40,38,39,36,37,]]

   >>> print ThreeDAddedMesh.faceVertexIDs
   [[ 0, 1, 4, 3,]
    [ 1, 2, 5, 4,]
    [ 3, 4, 7, 6,]
    [ 4, 5, 8, 7,]
    [ 9,10,13,12,]
    [10,11,14,13,]
    [12,13,16,15,]
    [13,14,17,16,]
    [18,19,22,21,]
    [19,20,23,22,]
    [21,22,25,24,]
    [22,23,26,25,]
    [ 0, 1,10, 9,]
    [ 1, 2,11,10,]
    [ 3, 4,13,12,]
    [ 4, 5,14,13,]
    [ 6, 7,16,15,]
    [ 7, 8,17,16,]
    [ 9,10,19,18,]
    [10,11,20,19,]
    [12,13,22,21,]
    [13,14,23,22,]
    [15,16,25,24,]
    [16,17,26,25,]
    [ 0, 3,12, 9,]
    [ 1, 4,13,10,]
    [ 2, 5,14,11,]
    [ 3, 6,15,12,]
    [ 4, 7,16,13,]
    [ 5, 8,17,14,]
    [ 9,12,21,18,]
    [10,13,22,19,]
    [11,14,23,20,]
    [12,15,24,21,]
    [13,16,25,22,]
    [14,17,26,23,]
    [ 2,27,28, 5,]
    [11,29,30,14,]
    [ 2,27,29,11,]
    [ 5,28,30,14,]
    [27,28,30,29,]]
     
   >>> InvalidMesh = ThreeDBaseMesh + basemesh
   Traceback (most recent call last):
   ...
   MeshAdditionError: Dimensions do not match
"""


__docformat__ = 'restructuredtext'

import Numeric
import MA

import fipy.meshes.common.mesh

from fipy.meshes.numMesh.face import Face
from fipy.meshes.numMesh.cell import Cell

import fipy.tools.array
import fipy.tools.vector as vector
from fipy.tools.dimensions.physicalField import PhysicalField


# Necessary because LLNL hires stupidheads
def MAtake(array, indices, fill = 0, axis = 0):
    tmp = MA.take(array, MA.filled(indices, fill), axis = axis)
    if indices.mask() is not None and tmp.shape != indices.mask().shape:
	mask = MA.repeat(indices.mask()[...,Numeric.NewAxis],tmp.shape[-1],len(tmp.shape)-1)
        if tmp.mask() is not None:
            mask = Numeric.logical_or(tmp.mask(), mask)
    else:
	mask = indices.mask()
    return MA.array(data = tmp, mask = mask)

class Mesh(fipy.meshes.common.mesh.Mesh):
    def __init__(self, vertexCoords, faceVertexIDs, cellFaceIDs):
        """faceVertexIds and cellFacesIds must be padded with minus ones."""

        self.vertexCoords = vertexCoords
        self.faceVertexIDs = MA.array(faceVertexIDs)
        self.cellFaceIDs = MA.array(cellFaceIDs)

	fipy.meshes.common.mesh.Mesh.__init__(self)
	
    """Topology methods"""

    def __add__(self, other, smallNumber = 1e-15):
        if(isinstance(other, Mesh)):
            return self.meshAdd(other, smallNumber)
        else:
            return self.translate(other)

    def __mul__(self, other):
        return self.dilate(other)

    def meshAdd(self, other, smallNumber):
        a = self.getAddedMeshValues(other, smallNumber)
        return Mesh(a[0], a[1], a[2])

    def getAddedMeshValues(self, other, smallNumber):
        """
        Returns a tuple with 3 elements: the new mesh vertexCoords, faceVertexIDs, and cellFaceIDs, in that order.
        """
        
        MeshAdditionError = "MeshAdditionError"
        selfNumFaces = self.faceVertexIDs.shape[0]
        selfNumCells = self.cellFaceIDs.shape[0]
        selfNumVertices = self.vertexCoords.shape[0]
        otherNumFaces = other.faceVertexIDs.shape[0]
        otherNumCells = other.cellFaceIDs.shape[0]
        otherNumVertices = other.vertexCoords.shape[0]
        ## check dimensions
        if(self.vertexCoords.shape[1] != other.vertexCoords.shape[1]):
            raise MeshAdditionError, "Dimensions do not match"
        else:
            dimensions = self.vertexCoords.shape[1]
        ## compute vertex correlates
        vertexCorrelates = {}
        for i in range(selfNumVertices):
            for j in range(otherNumVertices):
                diff = self.vertexCoords[i] - other.vertexCoords[j]
                diff = Numeric.array(diff)
                if (sum(diff ** 2) < smallNumber):
                    vertexCorrelates[j] = i
        if (vertexCorrelates == {}):
            raise MeshAdditionError, "Vertices are not aligned"
        ## compute face correlates
        faceCorrelates = {}
        for i in range(otherNumFaces):
            currFace = other.faceVertexIDs[i]
            keepGoing = 1
            currIndex = 0
            for item in currFace:
                if(vertexCorrelates.has_key(item)):
                    currFace[currIndex] = vertexCorrelates[item]
                    currIndex = currIndex + 1
                else:
                    keepGoing = 0
            if(keepGoing == 1):
                for j in range(selfNumFaces):
                    if (self.equalExceptOrder(currFace, self.faceVertexIDs[j])):
                        faceCorrelates[i] = j
        if(faceCorrelates == {}):
            raise MeshAdditionError, "Faces are not aligned"

        faceIndicesToAdd = ()
        for i in range(otherNumFaces):
            if(not faceCorrelates.has_key(i)):
                faceIndicesToAdd = faceIndicesToAdd + (i,)
        vertexIndicesToAdd = ()
        for i in range(otherNumVertices):
            if(not vertexCorrelates.has_key(i)):
                vertexIndicesToAdd = vertexIndicesToAdd + (i,)

        ##compute the full face and vertex correlation list
        a = selfNumFaces
        for i in faceIndicesToAdd:
            faceCorrelates[i] = a
            a = a + 1
        b = selfNumVertices
        for i in vertexIndicesToAdd:
            vertexCorrelates[i] = b
            b = b + 1

        ## compute what the cells are that we need to add
        cellsToAdd = Numeric.ones((other.cellFaceIDs.shape[0], self.cellFaceIDs.shape[1]))
        cellsToAdd = -1 * cellsToAdd
        
        for i in range(len(other.cellFaceIDs)):
            for j in range(len(other.cellFaceIDs[i])):
                cellsToAdd[i, j] = faceCorrelates[other.cellFaceIDs[i, j]]

        cellsToAdd = MA.masked_values(cellsToAdd, -1)


        ## compute what the faces are that we need to add
        facesToAdd = Numeric.take(other.faceVertexIDs, faceIndicesToAdd)
        for i in range(len(facesToAdd)):
            for j in range(len(facesToAdd[i])):
                facesToAdd[i, j] = vertexCorrelates[facesToAdd[i, j]]

        ## compute what the vertices are that we need to add
        verticesToAdd = Numeric.take(other.vertexCoords, vertexIndicesToAdd)

        newMeshCells = MA.concatenate((self.cellFaceIDs, cellsToAdd))
        newMeshFaces = Numeric.concatenate((self.faceVertexIDs, facesToAdd))
        newMeshVertices = Numeric.concatenate((self.vertexCoords, verticesToAdd))
        return (newMeshVertices, newMeshFaces, newMeshCells)

    def equalExceptOrder(self, first, second):
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
    

    def translate(self, vector):
        newCoords = self.vertexCoords + vector
        newmesh = Mesh(newCoords, Numeric.array(self.faceVertexIDs), Numeric.array(self.cellFaceIDs))
        return newmesh

    def dilate(self, factor):
        newCoords = self.vertexCoords * factor
        newmesh = Mesh(newCoords, Numeric.array(self.faceVertexIDs), Numeric.array(self.cellFaceIDs))
        return newmesh
    
    def calcTopology(self):
        self.dim = len(self.vertexCoords[0])
        self.numberOfFaces = len(self.faceVertexIDs)
        self.numberOfCells = len(self.cellFaceIDs)
        self.calcFaceCellIDs()
	
	fipy.meshes.common.mesh.Mesh.calcTopology(self)


    """calc Topology methods"""

    def calcFaceCellIDs(self):
        array = MA.indices((len(self.cellFaceIDs), len(self.cellFaceIDs[0])))[0]
        array = MA.array(data = array, mask = self.cellFaceIDs.mask()).flat
        cellFaceIDsFlat = MA.ravel(self.cellFaceIDs)
        firstRow = MA.zeros(self.numberOfFaces)
        secondRow = MA.zeros(self.numberOfFaces)
        MA.put(firstRow, cellFaceIDsFlat[::-1], array[::-1])
        MA.put(secondRow, cellFaceIDsFlat, array)
        secondRow = MA.array(data = secondRow, mask = (secondRow == firstRow))
        self.faceCellIDs = MA.zeros((len(firstRow),2))
        self.faceCellIDs[:,0] = firstRow[:]
        self.faceCellIDs[:,1] = secondRow[:]
        

    def calcInteriorAndExteriorFaceIDs(self):
        self.exteriorFaceIDs = Numeric.nonzero(self.faceCellIDs[:,1].mask())
        self.interiorFaceIDs = Numeric.nonzero(Numeric.logical_not(self.faceCellIDs[:,1].mask()))

    def calcInteriorAndExteriorCellIDs(self):
        try:
            import sets
            self.exteriorCellIDs = sets.Set(MA.take(self.faceCellIDs[:,0],self.exteriorFaceIDs))
            self.interiorCellIDs = list(sets.Set(range(self.numberOfCells)) - self.exteriorCellIDs)
            self.exteriorCellIDs = list(self.exteriorCellIDs)
        except:
            self.exteriorCellIDs = Numeric.take(self.faceCellIDs[:,0], self.exteriorFaceIDs)
            tmp = Numeric.zeros(self.numberOfCells)
            Numeric.put(tmp, self.exteriorCellIDs, Numeric.ones(len(self.exteriorCellIDs)))
            self.exteriorCellIDs = Numeric.nonzero(tmp)            
            self.interiorCellIDs = Numeric.nonzero(Numeric.logical_not(tmp))
            
    def calcCellToFaceOrientations(self):
	tmp = MAtake(self.faceCellIDs[:,0], self.cellFaceIDs)
	self.cellToFaceOrientations = (tmp == MA.indices(tmp.shape)[0]) * 2 - 1

    def calcAdjacentCellIDs(self):
        self.adjacentCellIDs = (MA.filled(self.faceCellIDs[:,0]), MA.filled(MA.where(self.faceCellIDs[:,1].mask(), self.faceCellIDs[:,0], self.faceCellIDs[:,1])))

    def calcCellToCellIDs(self):        
        self.cellToCellIDs = MAtake(self.faceCellIDs, self.cellFaceIDs)
        self.cellToCellIDs = MA.where(self.cellToFaceOrientations == 1, self.cellToCellIDs[:,:,1], self.cellToCellIDs[:,:,0])

    """get Topology methods"""

    def getVertexCoords(self):
        return self.vertexCoords

    def getExteriorFaces(self):
	return [Face(self, id) for id in self.exteriorFaceIDs]

    def getInteriorFaces(self):
	return [Face(self, id) for id in self.interiorFaceIDs]
	
    def getFaceCellIDs(self):
        return self.faceCellIDs

    def getFaces(self):
        return [Face(self, id) for id in Numeric.arange(self.numberOfFaces)]

    def getCellsByID(self, ids = None):
	if ids is None:
	    ids = range(self.numberOfCells) 
	return [Cell(self, id) for id in ids]
	
    def getMaxFacesPerCell(self):
        return len(self.cellFaceIDs[0])

    """Geometry methods"""

    def calcGeometry(self):
	self.calcFaceCenters()
	fipy.meshes.common.mesh.Mesh.calcGeometry(self)
        self.calcCellNormals()
        
    """calc geometry methods"""

    def calcFaceAreas(self):
        faceVertexIDs = MA.filled(self.faceVertexIDs, -1)
        substitute = Numeric.reshape(Numeric.repeat(faceVertexIDs[:,0],len(faceVertexIDs[0])), Numeric.shape(faceVertexIDs))
        if (self.faceVertexIDs.mask()):
            faceVertexIDs = Numeric.where(self.faceVertexIDs.mask(), substitute, faceVertexIDs)    
        faceVertexCoords = Numeric.take(self.vertexCoords, faceVertexIDs)
        faceOrigins = Numeric.repeat(faceVertexCoords[:,0], len(faceVertexIDs[0]))
        faceOrigins = Numeric.reshape(faceOrigins, MA.shape(faceVertexCoords))
        faceVertexCoords = faceVertexCoords - faceOrigins
        left = range(len(faceVertexIDs[0]))
        right = left[1:] + [left[0]]
        cross = Numeric.sum(fipy.tools.array.crossProd(faceVertexCoords, Numeric.take(faceVertexCoords, right, 1)), 1)
        self.faceAreas = fipy.tools.array.sqrtDot(cross, cross) / 2.

    def calcFaceCenters(self):
        faceVertexIDs = MA.filled(self.faceVertexIDs, 0)
        faceVertexCoords = Numeric.take(self.vertexCoords, faceVertexIDs)
        if self.faceVertexIDs.mask() == None:
            faceVertexCoordsMask = Numeric.zeros(Numeric.shape(faceVertexCoords))
        else:
            faceVertexCoordsMask = Numeric.reshape(Numeric.repeat(self.faceVertexIDs.mask().flat, self.dim), Numeric.shape(faceVertexCoords))
        faceVertexCoords = MA.array(data = faceVertexCoords, mask = faceVertexCoordsMask)

	self.faceCenters = MA.filled(MA.average(faceVertexCoords, axis = 1))

    def calcFaceNormals(self):
        faceVertexIDs = MA.filled(self.faceVertexIDs, 0)
        faceVertexCoords = Numeric.take(self.vertexCoords, faceVertexIDs)
        t1 = faceVertexCoords[:,1,:] - faceVertexCoords[:,0,:]
        t2 = faceVertexCoords[:,2,:] - faceVertexCoords[:,1,:]
        norm = fipy.tools.array.crossProd(t1, t2)
        sqrtDot = fipy.tools.array.sqrtDot(norm, norm)
        norm[:,0] = norm[:,0] / sqrtDot
        norm[:,1] = norm[:,1] / sqrtDot
        norm[:,2] = norm[:,2] / sqrtDot
        
        self.faceNormals = -norm
	
    def calcOrientedFaceNormals(self):
	self.orientedFaceNormals = self.faceNormals
	
    def calcCellVolumes(self):
	tmp = self.faceCenters[:,0] * self.faceAreas * self.faceNormals[:,0]
	tmp = MAtake(tmp, self.cellFaceIDs) * self.cellToFaceOrientations
        self.cellVolumes = MA.filled(MA.sum(tmp, 1))
        
    def calcCellCenters(self):
	tmp = MAtake(self.faceCenters, self.cellFaceIDs)
	self.cellCenters = MA.filled(MA.average(tmp, 1))
	
    def calcFaceToCellDistances(self):
	tmp = MAtake(self.cellCenters, self.faceCellIDs)
	tmp -= MA.repeat(self.faceCenters[:,Numeric.NewAxis,...], 2, 1)
	self.faceToCellDistances = MA.sqrt(MA.sum(tmp * tmp,2))

    def calcCellDistances(self):
	tmp = MAtake(self.cellCenters, self.faceCellIDs)
	tmp = tmp[:,1] - tmp[:,0]
        self.cellDistanceVectors = tmp
	tmp = MA.sqrt(MA.sum(tmp * tmp,1))
	self.cellDistances = MA.filled(MA.where(tmp.mask(), self.faceToCellDistances[:,0], tmp))

    def calcFaceToCellDistanceRatio(self):
        dAP = self.getCellDistances()
        dFP = self.getFaceToCellDistances()[:,0]
	
	self.faceToCellDistanceRatio = MA.filled(dFP / dAP)

    def calcAreaProjections(self):
        self.areaProjections = self.getFaceNormals() * self.getFaceAreas()[:,Numeric.NewAxis]
	
    def calcOrientedAreaProjections(self):
	self.orientedAreaProjections = self.areaProjections

    def calcFaceTangents(self):
        faceVertexCoord = Numeric.take(self.vertexCoords, self.faceVertexIDs[:,0])
        tmp = self.faceCenters - faceVertexCoord
        self.faceTangents1 = tmp / fipy.tools.array.sqrtDot(tmp, tmp)[:,Numeric.NewAxis]  
        tmp = fipy.tools.array.crossProd(self.faceTangents1, self.faceNormals)
        self.faceTangents2 = tmp / fipy.tools.array.sqrtDot(tmp, tmp)[:,Numeric.NewAxis]
        
    def calcCellToCellDistances(self):
        self.cellToCellDistances = MAtake(self.cellDistances, self.getCellFaceIDs())

    def calcCellNormals(self):
        cellNormals = MAtake(self.getFaceNormals(), self.getCellFaceIDs())
        cellFaceCellIDs = MAtake(self.faceCellIDs[:,0], self.cellFaceIDs)
        cellIDs = Numeric.reshape(Numeric.repeat(Numeric.arange(self.getNumberOfCells()), self.getMaxFacesPerCell()), cellFaceCellIDs.shape)
        direction = (cellFaceCellIDs == cellIDs) * 2 - 1
        self.cellNormals =  direction[:,:,Numeric.NewAxis] * cellNormals
                         
    """get geometry methods"""

    def getFaceCenters(self):
        return self.faceCenters

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
            'vertexCoords' : self.vertexCoords,            
            'faceVertexIDs' : self.faceVertexIDs,
            'cellFaceIDs' : self.cellFaceIDs }
        return dict

    def __setstate__(self, dict):
        self.__init__(dict['vertexCoords'], dict['faceVertexIDs'], dict['cellFaceIDs'])
        
                      
### test test test    

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
