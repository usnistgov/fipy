#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
 #                                    created: 11/10/03 {2:44:42 PM} 
 #                                last update: 7/12/05 {1:12:34 PM} 
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

import Numeric
import MA

import fipy.meshes.common.mesh

from fipy.meshes.numMesh.face import Face
from fipy.meshes.numMesh.cell import Cell

from fipy.tools import numerix
from fipy.tools import vector
from fipy.tools.dimensions.physicalField import PhysicalField

class Mesh(fipy.meshes.common.mesh.Mesh):
    """Generic mesh class using Numeric to do the calculations

        Meshes contain cells, faces, and vertices.

        This is built for a non-mixed element mesh.
    """

    def __init__(self, vertexCoords, faceVertexIDs, cellFaceIDs):
        """faceVertexIds and cellFacesIds must be padded with minus ones."""

        self.vertexCoords = vertexCoords
        self.faceVertexIDs = MA.array(faceVertexIDs)
        self.cellFaceIDs = MA.array(cellFaceIDs)

        fipy.meshes.common.mesh.Mesh.__init__(self)
	
    """Topology methods"""

    def __add__(self, other):
        if(isinstance(other, Mesh)):
            return self._concatenate(other, smallNumber = 1e-15)
        else:
            return self._translate(other)

    def __mul__(self, factor):
        newCoords = self.vertexCoords * factor
        newmesh = Mesh(newCoords, Numeric.array(self.faceVertexIDs), Numeric.array(self.cellFaceIDs))
        return newmesh

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

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(nx = 2, ny = 2, dx = 1., dy = 1.)

           >>> print mesh._getCellFaceIDs()
           [[ 0, 7, 2, 6,]
            [ 1, 8, 3, 7,]
            [ 2,10, 4, 9,]
            [ 3,11, 5,10,]]

           >>> mesh._connectFaces(mesh.getFacesLeft(), mesh.getFacesRight())

           >>> print mesh._getCellFaceIDs()
           [[ 0, 7, 2, 6,]
            [ 1, 6, 3, 7,]
            [ 2,10, 4, 9,]
            [ 3, 9, 5,10,]]
        
        """

        ## extract the face IDs for each set of faces
        ids0 = [face.getID() for face in faces0]
        ids1 = [face.getID() for face in faces1]

        ## check for errors

        ## check that faces are members of exterior faces
        from sets import Set
        assert Set(ids0).union(Set(ids1)).issubset(Set(self.getExteriorFaceIDs()))

        ## following assert checks number of faces are equal, normals are opposite and areas are the same
        assert Numeric.take(self.areaProjections, ids0) == Numeric.take(-self.areaProjections, ids1)

        ## extract the adjacent cells for both sets of faces
        faceCellIDs0 = self.faceCellIDs[:,0]
        faceCellIDs1 = self.faceCellIDs[:,1]
        ## set the new adjacent cells for `faces0`
        MA.put(faceCellIDs1, ids0, MA.take(faceCellIDs0, ids0))
        MA.put(faceCellIDs0, ids0, MA.take(faceCellIDs0, ids1))
        self.faceCellIDs[:,0] = faceCellIDs0
        self.faceCellIDs[:,1] = faceCellIDs1
        
        ## extract the face to cell distances for both sets of faces
        faceToCellDistances0 = self.faceToCellDistances[:,0]
        faceToCellDistances1 = self.faceToCellDistances[:,1]
        ## set the new faceToCellDistances for `faces0`
        MA.put(faceToCellDistances1, ids0, MA.take(faceToCellDistances0, ids0))
        MA.put(faceToCellDistances0, ids0, MA.take(faceToCellDistances0, ids1))
        self.faceToCellDistances[:,0] = faceToCellDistances0
        self.faceToCellDistances[:,1] = faceToCellDistances1

        ## calculate new cell distances and add them to faces0
        Numeric.put(self.cellDistances, ids0, MA.take(faceToCellDistances0 + faceToCellDistances1, ids0))

        ## change the direction of the face normals for faces0
        for dim in range(self.getDim()):
            faceNormals = self.faceNormals[:,dim].copy()
            Numeric.put(faceNormals, ids0, MA.take(faceNormals, ids1))
            self.faceNormals[:,dim] = faceNormals

        ## Cells that are adjacent to faces1 are changed to point at faces0
        ## get the cells adjacent to faces1
        faceCellIDs = MA.take(self.faceCellIDs[:,0], ids1)
        ## get all the adjacent faces for those particular cells
        cellFaceIDs = MA.take(self.cellFaceIDs[:], faceCellIDs)
        for i in range(len(cellFaceIDs[0,:])):
            ## if the faces is a member of faces1 then change the face to point at
            ## faces0
            cellFaceIDs[:,i] = MA.where(cellFaceIDs[:,i] == ids1,
                                        ids0,
                                        cellFaceIDs[:,i])
            ## add those faces back to the main self.cellFaceIDs
            tmp = self.cellFaceIDs[:,i]
            MA.put(tmp, faceCellIDs, cellFaceIDs[:,i])
            self.cellFaceIDs[:,i] = tmp

        ## calculate new topology
        fipy.meshes.common.mesh.Mesh._calcTopology(self)

        ## calculate new geometry
        self._calcFaceToCellDistanceRatio()
        self._calcCellToCellDistances()
        self._calcScaledGeometry()
        self._calcFaceAspectRatios()
        
        
    def _getAddedMeshValues(self, other, smallNumber):
        """
        Returns a `dictionary` with 3 elements: the new mesh vertexCoords, faceVertexIDs, and cellFaceIDs.
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
                    if (self._equalExceptOrder(currFace, self.faceVertexIDs[j])):
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

        return {
            'vertexCoords': Numeric.concatenate((self.vertexCoords, verticesToAdd)), 
            'faceVertexIDs': Numeric.concatenate((self.faceVertexIDs, facesToAdd)), 
            'cellFaceIDs': MA.concatenate((self.cellFaceIDs, cellsToAdd))
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
        newmesh = Mesh(newCoords, Numeric.array(self.faceVertexIDs), Numeric.array(self.cellFaceIDs))
        return newmesh

    def _calcTopology(self):
        self.dim = len(self.vertexCoords[0])
        self.numberOfFaces = len(self.faceVertexIDs)
        self.numberOfCells = len(self.cellFaceIDs)
        self._calcFaceCellIDs()
	
        fipy.meshes.common.mesh.Mesh._calcTopology(self)


    """calc Topology methods"""

    def _calcFaceCellIDs(self):
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
        

    def _calcInteriorAndExteriorFaceIDs(self):
        self.exteriorFaceIDs = Numeric.nonzero(self.faceCellIDs[:,1].mask())
        self.interiorFaceIDs = Numeric.nonzero(Numeric.logical_not(self.faceCellIDs[:,1].mask()))

    def _calcInteriorAndExteriorCellIDs(self):
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
            
    def _calcCellToFaceOrientations(self):
	tmp = numerix.MAtake(self.faceCellIDs[:,0], self.cellFaceIDs)
	self.cellToFaceOrientations = (tmp == MA.indices(tmp.shape)[0]) * 2 - 1

    def _calcAdjacentCellIDs(self):
        self.adjacentCellIDs = (MA.filled(self.faceCellIDs[:,0]), MA.filled(MA.where(self.faceCellIDs[:,1].mask(), self.faceCellIDs[:,0], self.faceCellIDs[:,1])))

    def _calcCellToCellIDs(self):        
        self.cellToCellIDs = numerix.MAtake(self.faceCellIDs, self.cellFaceIDs)
        self.cellToCellIDs = MA.where(self.cellToFaceOrientations == 1, self.cellToCellIDs[:,:,1], self.cellToCellIDs[:,:,0])
        
    def _calcNumPts(self, d, n = None, axis = "x"):
        """
        Calculate the number of cells along the specified axis, based
        on either the specified number or on the number elements in the
        cell  `d` spacings.
        
        Used by the `Grid` meshes.
        """
        if type(d) in [type(1), type(1.)]:
            n = n or 1
        else:
            n = n or len(d)
            if n != len(d):
                raise IndexError, "n%s != len(d%s)" % (axis, axis)
                
        return n

    def _calcVertexCoordinates(self, d, n):
        """
        Calculate the positions of the vertices along an axis, based on the 
        specified `Cell` `d` spacing or list of `d` spacings.
        
        Used by the `Grid` meshes.
        """
        x = Numeric.zeros((n + 1), 'd')
        x[1:] = d
        return Numeric.add.accumulate(x)

    """get Topology methods"""

    def getVertexCoords(self):
        return self.vertexCoords

    def getExteriorFaces(self):
	return [Face(self, id) for id in self.exteriorFaceIDs]

    def _getInteriorFaces(self):
	return [Face(self, id) for id in self.interiorFaceIDs]
	
    def getFaceCellIDs(self):
        return self.faceCellIDs

    def _getFaces(self):
        return [Face(self, id) for id in Numeric.arange(self.numberOfFaces)]

    def _getCellsByID(self, ids = None):
	if ids is None:
	    ids = range(self.numberOfCells) 
	return [Cell(self, id) for id in ids]
    
    def getFacesByID(self, ids = None):
	if ids is None:
	    ids = range(self.numberOfFaces) 
	return [Face(self, id) for id in ids]
	
    def _getMaxFacesPerCell(self):
        return len(self.cellFaceIDs[0])

    """Geometry methods"""

    def _calcGeometry(self):
	self._calcFaceCenters()
        fipy.meshes.common.mesh.Mesh._calcGeometry(self)
        self._calcCellNormals()
        
    """calc geometry methods"""

    def _calcFaceAreas(self):
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
        cross = Numeric.sum(numerix.crossProd(faceVertexCoords, Numeric.take(faceVertexCoords, right, 1)), 1)
        self.faceAreas = numerix.sqrtDot(cross, cross) / 2.

    def _calcFaceCenters(self):
        faceVertexIDs = MA.filled(self.faceVertexIDs, 0)
        faceVertexCoords = Numeric.take(self.vertexCoords, faceVertexIDs)
        if self.faceVertexIDs.mask() == None:
            faceVertexCoordsMask = Numeric.zeros(Numeric.shape(faceVertexCoords))
        else:
            faceVertexCoordsMask = Numeric.reshape(Numeric.repeat(self.faceVertexIDs.mask().flat, self.dim), Numeric.shape(faceVertexCoords))
        faceVertexCoords = MA.array(data = faceVertexCoords, mask = faceVertexCoordsMask)

	self.faceCenters = MA.filled(MA.average(faceVertexCoords, axis = 1))

    def _calcFaceNormals(self):
        faceVertexIDs = MA.filled(self.faceVertexIDs, 0)
        faceVertexCoords = Numeric.take(self.vertexCoords, faceVertexIDs)
        t1 = faceVertexCoords[:,1,:] - faceVertexCoords[:,0,:]
        t2 = faceVertexCoords[:,2,:] - faceVertexCoords[:,1,:]
        norm = numerix.crossProd(t1, t2)
        sqrtDot = numerix.sqrtDot(norm, norm)
        norm[:,0] = norm[:,0] / sqrtDot
        norm[:,1] = norm[:,1] / sqrtDot
        norm[:,2] = norm[:,2] / sqrtDot
        
        self.faceNormals = -norm
	
    def _calcOrientedFaceNormals(self):
	self.orientedFaceNormals = self.faceNormals
	
    def _calcCellVolumes(self):
	tmp = self.faceCenters[:,0] * self.faceAreas * self.faceNormals[:,0]
	tmp = numerix.MAtake(tmp, self.cellFaceIDs) * self.cellToFaceOrientations
        self.cellVolumes = MA.filled(MA.sum(tmp, 1))
        
    def _calcCellCenters(self):
	tmp = numerix.MAtake(self.faceCenters, self.cellFaceIDs)
	self.cellCenters = MA.filled(MA.average(tmp, 1))
	
    def _calcFaceToCellDistances(self):
	tmp = numerix.MAtake(self.cellCenters, self.faceCellIDs)
	tmp -= MA.repeat(self.faceCenters[:,Numeric.NewAxis,...], 2, 1)
	self.faceToCellDistances = MA.sqrt(MA.sum(tmp * tmp,2))

    def _calcCellDistances(self):
	tmp = numerix.MAtake(self.cellCenters, self.faceCellIDs)
	tmp = tmp[:,1] - tmp[:,0]
        self.cellDistanceVectors = tmp
	tmp = MA.sqrt(MA.sum(tmp * tmp,1))
	self.cellDistances = MA.filled(MA.where(tmp.mask(), self.faceToCellDistances[:,0], tmp))

    def _calcFaceToCellDistanceRatio(self):
        dAP = self._getCellDistances()
        dFP = self._getFaceToCellDistances()[:,0]
	
	self.faceToCellDistanceRatio = MA.filled(dFP / dAP)

    def _calcAreaProjections(self):
        self.areaProjections = self._getFaceNormals() * self._getFaceAreas()[:,Numeric.NewAxis]
	
    def _calcOrientedAreaProjections(self):
	self.orientedAreaProjections = self.areaProjections

    def _calcFaceTangents(self):
        faceVertexCoord = Numeric.take(self.vertexCoords, self.faceVertexIDs[:,0])
        tmp = self.faceCenters - faceVertexCoord
        self.faceTangents1 = tmp / numerix.sqrtDot(tmp, tmp)[:,Numeric.NewAxis]  
        tmp = numerix.crossProd(self.faceTangents1, self.faceNormals)
        self.faceTangents2 = tmp / numerix.sqrtDot(tmp, tmp)[:,Numeric.NewAxis]
        
    def _calcCellToCellDistances(self):
        self.cellToCellDistances = numerix.MAtake(self.cellDistances, self._getCellFaceIDs())

    def _calcCellNormals(self):
        cellNormals = numerix.MAtake(self._getFaceNormals(), self._getCellFaceIDs())
        cellFaceCellIDs = numerix.MAtake(self.faceCellIDs[:,0], self.cellFaceIDs)
        cellIDs = Numeric.reshape(Numeric.repeat(Numeric.arange(self.getNumberOfCells()), self._getMaxFacesPerCell()), cellFaceCellIDs.shape)
        direction = (cellFaceCellIDs == cellIDs) * 2 - 1
        self.cellNormals =  direction[:,:,Numeric.NewAxis] * cellNormals
                         
    """get geometry methods"""

    def getFaceCenters(self):
        return self.faceCenters

    def _getCellVertexIDs(self):

        ## Get all the vertices from all the faces for each cell
        from fipy.tools.numerix import MAtake
        NCells = self.getNumberOfCells()
        cellFaceVertices = MAtake(self.faceVertexIDs, self.cellFaceIDs)

        ## get a sorted list of vertices for each cell 
        cellVertexIDs = MA.reshape(cellFaceVertices.flat, (NCells, -1))
        cellVertexIDs = MA.sort(cellVertexIDs, axis = 1, fill_value = -1)
        cellVertexIDs = cellVertexIDs[:,::-1]

        ## get a unique sorted list of vertices for each cell
        newmask = MA.getmaskarray(cellVertexIDs).copy()
        for i in range(len(cellVertexIDs[0]) - 1):
            newmask[:,i + 1] = MA.where(newmask[:,i + 1],
                                        newmask[:,i + 1],
                                        
                                        MA.where(MA.filled(cellVertexIDs)[:,i] == MA.filled(cellVertexIDs)[:,i + 1],
                                                 1,
                                                 newmask[:,i + 1]))

        cellVertexIDs = MA.masked_array(cellVertexIDs, newmask)
        cellVertexIDs = MA.sort(cellVertexIDs, axis = 1, fill_value = -1)
        cellVertexIDs = cellVertexIDs[:,::-1]

        ## resize the array to remove extra masked values
        length = len(cellVertexIDs[0]) - min(Numeric.sum(MA.getmaskarray(cellVertexIDs), axis = 1))
        return cellVertexIDs[:, :length]

## Below is an ordered version of _getCellVertexIDs()
## It works for the test case in this file (other than the ordering, obviously)
## I've left it in as it may be useful when we need ordered vertices for cells
    
##    def _getOrderedCellVertexIDs(self):

##        ## Get all the vertices from all the faces for each cell
##        from fipy.tools.numerix import MAtake
##        cellFaceVertices = MAtake(self.faceVertexIDs, self.cellFaceIDs)

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


##        length = len(cellVertexIDs[0]) - min(Numeric.sum(MA.getmaskarray(cellVertexIDs), axis = 1))
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
            'vertexCoords' : self.vertexCoords,            
            'faceVertexIDs' : self.faceVertexIDs,
            'cellFaceIDs' : self.cellFaceIDs }
        return dict

    def __setstate__(self, dict):
        Mesh.__init__(self, dict['vertexCoords'], dict['faceVertexIDs'], dict['cellFaceIDs'])
##        self.__init__(dict['vertexCoords'], dict['faceVertexIDs'], dict['cellFaceIDs'])
     
    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.
        
            >>> dx = 2.
            >>> dy = 1.23456
            >>> dz = 1.e-1
            
            >>> vertices = Numeric.array(((0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
            ...                           (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
            ...                           (2., 0., 0.), (2., 0., 1.)))
            >>> vertices *= Numeric.array((dx, dy, dz))
            >>> faces = MA.masked_values(((0, 1, 2, 3), (7, 6, 5, 4),
            ...                           (3, 7, 4, 0), (5, 6, 2, 1),
            ...                           (1, 0, 4, 5), (3, 2, 6, 7),
            ...                           (1, 8, 2, -1), (9, 5, 6, -1), (8, 1, 5, 9), (8, 9, 6, 2)),-1)
            >>> cells = MA.masked_values(((0, 1, 2, 3, 4, 5),
            ...                           (3 , 6, 7, 8, 9, -1)), -1)

            >>> mesh = Mesh(vertexCoords = vertices, faceVertexIDs = faces, cellFaceIDs = cells)

            >>> externalFaces = Numeric.array((0, 1, 2, 4, 5, 6, 7, 8, 9))
            >>> numerix.allequal(externalFaces, [face.getID() for face in mesh.getExteriorFaces()])
            1

            >>> internalFaces = Numeric.array((3,))
            >>> numerix.allequal(internalFaces, [face.getID() for face in mesh._getInteriorFaces()])
            1

            >>> import MA
            >>> faceCellIds = MA.masked_values(((0, -1), (0, -1), (0, -1),
            ...                                 (0, 1), (0, -1), (0, -1),
            ...                                 (1, -1), (1, -1), (1, -1), (1, -1)), -1)
            >>> numerix.allequal(faceCellIds, mesh.getFaceCellIDs())
            1
            
            >>> dxdy = dx * dy
            >>> dxdz = dx * dz
            >>> dydz = dy * dz
            >>> faceAreas = Numeric.array((dxdy, dxdy, dydz, dydz, dxdz, dxdz,
            ...                            dxdy/2., dxdy/2., dxdz, Numeric.sqrt(dx**2 + dy**2) * dz))
            >>> numerix.allclose(faceAreas, mesh._getFaceAreas(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = Numeric.take(vertices, MA.filled(faces, 0))
            >>> faceCenters = faceCoords[:,0] + faceCoords[:,1] + faceCoords[:,2] + faceCoords[:,3]
            >>> numVex = Numeric.array((4., 4., 4., 4., 4., 4., 3., 3., 4., 4.))
            >>> faceCenters /= numVex[..., Numeric.NewAxis]
            >>> numerix.allclose(faceCenters, mesh.getFaceCenters(), atol = 1e-10, rtol = 1e-10)
            1

            >>> faceNormals = Numeric.array(((0., 0., -1.),
            ...                              (0., 0., 1.),
            ...                              (-1, 0., 0.),
            ...                              (1. , 0., 0.),
            ...                              (0, -1., 0.),
            ...                              (0, 1., 0.),
            ...                              (0., 0., -1.),
            ...                              (0., 0., 1.),
            ...                              (0., -1., 0.),
            ...                              (dy / Numeric.sqrt(dy**2 + dx**2), dx / Numeric.sqrt(dy**2 + dx**2), 0.)))
            >>> numerix.allclose(faceNormals, mesh._getFaceNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = MA.masked_values(((1, 1, 1, 1, 1, 1),
            ...                                            (-1, 1, 1, 1, 1, -2)), -2)
            >>> numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            1
                                             
            >>> cellVolumes = Numeric.array((dx*dy*dz, dx*dy*dz / 2.))
            >>> numerix.allclose(cellVolumes, mesh.getCellVolumes(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellCenters = Numeric.array(((dx/2.,dy/2.,dz/2.), (dx+dx/3.,dy/3.,dz/2.)))
            >>> numerix.allclose(cellCenters, mesh.getCellCenters(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> d1 = Numeric.sqrt((dx / 3.)**2 + (dy / 6.)**2)
            >>> d2 = Numeric.sqrt((dx / 6.)**2 + (dy / 3.)**2)
            >>> d3 = Numeric.sqrt((dx / 6.)**2 + (dy / 6.)**2)
            >>> d4 = Numeric.sqrt((5 * dx / 6.)**2 + (dy / 6.)**2)
            >>> faceToCellDistances = MA.masked_values(((dz / 2., -1),
            ...                                         (dz / 2., -1),
            ...                                         (dx / 2., -1),
            ...                                         (dx / 2., d1),
            ...                                         (dy / 2., -1),
            ...                                         (dy / 2., -1),
            ...                                         (dz / 2., -1),
            ...                                         (dz / 2., -1),
            ...                                         (d2, -1),
            ...                                         (d3, -1)), -1)
            >>> numerix.allclose(faceToCellDistances, mesh._getFaceToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> cellDistances = Numeric.array((dz / 2., dz / 2., dx / 2.,
            ...                                d4,
            ...                                dy / 2., dy / 2., dz / 2., dz / 2.,
            ...                                d2,
            ...                                d3))
            >>> numerix.allclose(cellDistances, mesh._getCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceToCellDistanceRatios = faceToCellDistances[...,0] / cellDistances
            >>> numerix.allclose(faceToCellDistanceRatios, mesh._getFaceToCellDistanceRatio(), atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas[...,Numeric.NewAxis]
            >>> numerix.allclose(areaProjections, mesh._getAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> v1 = Numeric.take(vertices, faces[...,0])
            >>> tmp = faceCenters - v1
            >>> tangents1 = tmp / numerix.sqrtDot(tmp, tmp)[...,Numeric.NewAxis]
            >>> numerix.allclose(tangents1, mesh._getFaceTangents1(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tmp = numerix.crossProd(tangents1, faceNormals)
            >>> tangents2 = tmp / numerix.sqrtDot(tmp, tmp)[...,Numeric.NewAxis]
            >>> numerix.allclose(tangents2, mesh._getFaceTangents2(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((-1, -1, -1, 1, -1, -1),
            ...                                   (0, -1, -1, -1, -1, -1)), -1)
            >>> numerix.allequal(cellToCellIDs, mesh._getCellToCellIDs())
            1

            >>> cellToCellDistances = MA.masked_values(((dz / 2., dz / 2., dx / 2., d4, dy / 2., dy / 2.),
            ...                                         (d4, -1, -1, -1, -1, -1)), -1)
            >>> numerix.allclose(cellToCellDistances, mesh._getCellToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1

            >>> interiorCellIDs = Numeric.array(())
            >>> numerix.allequal(interiorCellIDs, mesh._getInteriorCellIDs())
            1

            >>> exteriorCellIDs = Numeric.array((0, 1))
            >>> numerix.allequal(exteriorCellIDs, mesh._getExteriorCellIDs())
            1

            >>> xnor = Numeric.array((1, 0, 0))
            >>> ynor = Numeric.array((0, 1, 0))
            >>> znor = Numeric.array((0, 0, 1))
            >>> nor =  Numeric.array((dy, dx, 0)) / Numeric.sqrt(dx**2 + dy**2)
            >>> cellNormals = MA.masked_values(((-znor, znor, -xnor, xnor, -ynor, ynor),
            ...                                 (-xnor, -znor, znor, -ynor, nor, (-1000, -1000, -1000))), -1000)
            >>> numerix.allclose(cellNormals, mesh._getCellNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellAreaProjections = MA.masked_values(((-znor * dx * dy, znor * dx * dy, -xnor * dy * dz,
            ...                                           xnor * dy * dz, -ynor * dx * dz, ynor * dx * dz),
            ...                                         (-xnor * dy * dz, -znor * dx * dy / 2, znor * dx * dy / 2, 
            ...                                          -ynor * dx * dz, nor * Numeric.sqrt(dx**2 + dy**2) * dz, 
            ...                                          (-1000, -1000, -1000))), -1000)
            >>> numerix.allclose(cellAreaProjections, mesh._getCellAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellVertexIDs = MA.masked_values(((7, 6, 5, 4, 3, 2, 1, 0), (9, 8, 6, 5, 2, 1, -1000, -1000)), -1000)
            >>> numerix.allclose(cellVertexIDs, mesh._getCellVertexIDs())
            1

            >>> import tempfile
            >>> import os
            >>> from fipy.tools import dump
            
            >>> (f, filename) = tempfile.mkstemp('.gz')
            >>> pickledMesh = dump.write(mesh, filename)
            
            >>> unpickledMesh = dump.read(filename)
            >>> os.close(f)
            >>> os.remove(filename)

            >>> numerix.allequal(mesh.getCellCenters(), unpickledMesh.getCellCenters())
            1
        """

                      
### test test test    

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
