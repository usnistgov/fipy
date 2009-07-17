#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
 #
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

        self.vertexCoords = vertexCoords
        self.faceVertexIDs = MA.masked_values(faceVertexIDs, -1)
        self.cellFaceIDs = MA.masked_values(cellFaceIDs, -1)

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

           >>> print (mesh._getCellFaceIDs() == [[0, 1, 2, 3],
           ...                                   [7, 8, 10, 11],
           ...                                   [2, 3, 4, 5],
           ...                                   [6, 7, 9, 10]]).flatten().all()
           True

           >>> mesh._connectFaces(numerix.nonzero(mesh.getFacesLeft()), numerix.nonzero(mesh.getFacesRight()))

           >>> print (mesh._getCellFaceIDs() == [[0, 1, 2, 3],
           ...                                   [7, 6, 10, 9],
           ...                                   [2, 3, 4, 5],
           ...                                   [6, 7, 9, 10]]).flatten().all()
           True

        """

        ## check for errors

        ## check that faces are members of exterior faces
        from fipy.variables.faceVariable import FaceVariable
        faces = FaceVariable(mesh=self, value=False)
        faces[faces0] = True
        faces[faces1] = True
        assert (faces | self.getExteriorFaces() == self.getExteriorFaces()).all()

        ## following assert checks number of faces are equal, normals are opposite and areas are the same
        assert numerix.alltrue(numerix.take(self.areaProjections, faces0, axis=1) 
                               == numerix.take(-self.areaProjections, faces1, axis=1))

        ## extract the adjacent cells for both sets of faces
        faceCellIDs0 = self.faceCellIDs[0]
        faceCellIDs1 = self.faceCellIDs[1]
        ## set the new adjacent cells for `faces0`
        MA.put(faceCellIDs1, faces0, MA.take(faceCellIDs0, faces0))
        MA.put(faceCellIDs0, faces0, MA.take(faceCellIDs0, faces1))
        self.faceCellIDs[0] = faceCellIDs0
        self.faceCellIDs[1] = faceCellIDs1
        
        ## extract the face to cell distances for both sets of faces
        faceToCellDistances0 = self.faceToCellDistances[0]
        faceToCellDistances1 = self.faceToCellDistances[1]
        ## set the new faceToCellDistances for `faces0`
        MA.put(faceToCellDistances1, faces0, MA.take(faceToCellDistances0, faces0))
        MA.put(faceToCellDistances0, faces0, MA.take(faceToCellDistances0, faces1))
        self.faceToCellDistances[0] = faceToCellDistances0
        self.faceToCellDistances[1] = faceToCellDistances1

        ## calculate new cell distances and add them to faces0
        numerix.put(self.cellDistances, faces0, MA.take(faceToCellDistances0 + faceToCellDistances1, faces0))

        ## change the direction of the face normals for faces0
        for dim in range(self.getDim()):
            faceNormals = self.faceNormals[dim].copy()
            numerix.put(faceNormals, faces0, MA.take(faceNormals, faces1))
            self.faceNormals[dim] = faceNormals

        ## Cells that are adjacent to faces1 are changed to point at faces0
        ## get the cells adjacent to faces1
        faceCellIDs = MA.take(self.faceCellIDs[0], faces1)
        ## get all the adjacent faces for those particular cells
        cellFaceIDs = numerix.take(self.cellFaceIDs, faceCellIDs, axis=1)
        for i in range(cellFaceIDs.shape[0]):
            ## if the faces is a member of faces1 then change the face to point at
            ## faces0
            cellFaceIDs[i] = MA.where(cellFaceIDs[i] == faces1,
                                      faces0,
                                      cellFaceIDs[i])
            ## add those faces back to the main self.cellFaceIDs
            tmp = self.cellFaceIDs[i]
            numerix.put(tmp, faceCellIDs, cellFaceIDs[i])
            self.cellFaceIDs[i] = tmp

        ## calculate new topology
        _CommonMesh._calcTopology(self)

        ## calculate new geometry
        self._calcFaceToCellDistanceRatio()
        self._calcCellToCellDistances()
        self._calcScaledGeometry()
        self._calcFaceAspectRatios()
        
    def _getConcatenableMesh(self):
        return self
        
    def _getAddedMeshValues(self, other, smallNumber):
        """
        Returns a `dictionary` with 3 elements: the new mesh vertexCoords, faceVertexIDs, and cellFaceIDs.
        """

        other = other._getConcatenableMesh()

        selfNumFaces = self.faceVertexIDs.shape[-1]
        selfNumVertices = self.vertexCoords.shape[-1]
        otherNumFaces = other.faceVertexIDs.shape[-1]
        otherNumVertices = other.vertexCoords.shape[-1]
        ## check dimensions
        if(self.vertexCoords.shape[0] != other.vertexCoords.shape[0]):
            raise MeshAdditionError, "Dimensions do not match"
        ## compute vertex correlates
        vertexCorrelates = {}
        for i in range(selfNumVertices):
            for j in range(otherNumVertices):
                diff = self.vertexCoords[...,i] - other.vertexCoords[...,j]
                diff = numerix.array(diff)
                if (sum(diff ** 2) < smallNumber):
                    vertexCorrelates[j] = i
        if (vertexCorrelates == {}):
            raise MeshAdditionError, "Vertices are not aligned"

        
         
        ## compute face correlates
        faceCorrelates = {}
        for i in range(otherNumFaces):
##          Seems to be overwriting other.faceVertexIDs with new numpy
##            currFace = other.faceVertexIDs[i]
##            currFace = other.faceVertexIDs[...,i].copy()
##          Changed this again as numpy 1.0.4 seems to have no copy method for
##          masked arrays.
            try:
                currFace = other.faceVertexIDs[...,i].copy()
            except:
                currFace = MA.array(other.faceVertexIDs[...,i], mask=MA.getmask(other.faceVertexIDs[...,i]))

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
                    if (self._equalExceptOrder(currFace, self.faceVertexIDs[...,j])):
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
        cellsToAdd = numerix.ones((self.cellFaceIDs.shape[0], other.cellFaceIDs.shape[-1]))
        cellsToAdd = -1 * cellsToAdd

        for j in range(other.cellFaceIDs.shape[-1]):
            for i in range(other.cellFaceIDs.shape[0]):
                cellsToAdd[i, j] = faceCorrelates[other.cellFaceIDs[i, j]]

        cellsToAdd = MA.masked_values(cellsToAdd, -1)


        ## compute what the faces are that we need to add
        facesToAdd = numerix.take(other.faceVertexIDs, faceIndicesToAdd, axis=1)

        for j in range(facesToAdd.shape[-1]):
            for i in range(facesToAdd.shape[0]):
                facesToAdd[i, j] = vertexCorrelates[facesToAdd[i, j]]

        ## compute what the vertices are that we need to add
        verticesToAdd = numerix.take(other.vertexCoords, vertexIndicesToAdd, axis=1)

        return {
            'vertexCoords': numerix.concatenate((self.vertexCoords, verticesToAdd), axis=1), 
            'faceVertexIDs': numerix.concatenate((self.faceVertexIDs, facesToAdd), axis=1), 
            'cellFaceIDs': MA.concatenate((self.cellFaceIDs, cellsToAdd), axis=1)
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
        newmesh = Mesh(newCoords, numerix.array(self.faceVertexIDs), numerix.array(self.cellFaceIDs))
        return newmesh

    def _calcTopology(self):
        self.dim = len(self.vertexCoords[...,0])
        if not hasattr(self, "numberOfFaces"):
            self.numberOfFaces = self.faceVertexIDs.shape[-1]
        if not hasattr(self, "numberOfCells"):
            self.numberOfCells = self.cellFaceIDs.shape[-1]
        if not hasattr(self, "globalNumberOfCells"):
            self.globalNumberOfCells = self.numberOfCells
        if not hasattr(self, "globalNumberOfFaces"):
            self.globalNumberOfFaces = self.numberOfFaces
        self._calcFaceCellIDs()
        
        _CommonMesh._calcTopology(self)


    """calc Topology methods"""

    def _calcFaceCellIDs(self):
        array = MA.array(MA.indices(self.cellFaceIDs.shape, 'l')[1], 
                         mask=MA.getmask(self.cellFaceIDs))
        self.faceCellIDs = MA.zeros((2, self.numberOfFaces), 'l')

        ## Nasty bug: MA.put(arr, ids, values) fills its ids and
        ## values arguments when masked!  This was not the behavior
        ## that was assumed when used below.  It was only working
        ## because the old fill value was 0 and the first element of
        ## the array needed to be 0 since the cell's face was
        ## 0. numerix.put() has been changed to deal with this
        ## properly.

##         MA.put(firstRow, cellFaceIDsFlat[::-1], array[::-1])
##         MA.put(secondRow, cellFaceIDsFlat, array)
        firstRow = self.faceCellIDs[0]
        secondRow = self.faceCellIDs[1]
        numerix.put(firstRow, self.cellFaceIDs[::-1,::-1], array[::-1,::-1])
        numerix.put(secondRow, self.cellFaceIDs, array)
        
        mask = ((False,) * self.numberOfFaces, (firstRow == secondRow))
        self.faceCellIDs = MA.sort(MA.array(self.faceCellIDs, mask = mask),
                                   axis=0)

    def _calcInteriorAndExteriorFaceIDs(self):
        from fipy.variables.faceVariable import FaceVariable
        mask = MA.getmask(self.faceCellIDs[1])
        self.exteriorFaces = FaceVariable(mesh=self, 
                                          value=mask)
        self.interiorFaces = FaceVariable(mesh=self, 
                                          value=numerix.logical_not(mask))

    def _calcInteriorAndExteriorCellIDs(self):
        try:
            import sets
            self.exteriorCellIDs = sets.Set(self.faceCellIDs[0, self.getExteriorFaces().getValue()])
            self.interiorCellIDs = list(sets.Set(range(self.numberOfCells)) - self.exteriorCellIDs)
            self.exteriorCellIDs = list(self.exteriorCellIDs)
        except:
            self.exteriorCellIDs = self.faceCellIDs[0, self.getExteriorFaces().getValue()]
            tmp = numerix.zeros(self.numberOfCells)
            numerix.put(tmp, self.exteriorCellIDs, numerix.ones(len(self.exteriorCellIDs)))
            self.exteriorCellIDs = numerix.nonzero(tmp)            
            self.interiorCellIDs = numerix.nonzero(numerix.logical_not(tmp))
            
    def _calcCellToFaceOrientations(self):
        tmp = numerix.take(self.faceCellIDs[0], self.cellFaceIDs)
        self.cellToFaceOrientations = (tmp == MA.indices(tmp.shape)[-1]) * 2 - 1

    def _calcAdjacentCellIDs(self):
        self.adjacentCellIDs = (MA.filled(self.faceCellIDs[0]), 
                                MA.filled(MA.where(MA.getmaskarray(self.faceCellIDs[1]), 
                                                   self.faceCellIDs[0], 
                                                   self.faceCellIDs[1])))

    def _calcCellToCellIDs(self):    
        self.cellToCellIDs = numerix.take(self.faceCellIDs, self.cellFaceIDs, axis=1)
        self.cellToCellIDs = MA.where(self.cellToFaceOrientations == 1, self.cellToCellIDs[1], self.cellToCellIDs[0])
        
    def _calcNumPts(self, d, n = None, axis = "x"):
        """
        Calculate the number of cells along the specified axis, based
        on either the specified number or on the number elements in the
        cell  `d` spacings.
        
        Used by the `Grid` meshes.

        This tests a bug that was occuring with PeriodicGrid1D when
        using a numpy float as the argument for the grid spacing.

           >>> from fipy.meshes.periodicGrid1D import PeriodicGrid1D
           >>> PeriodicGrid1D(nx=2, dx=numerix.float32(1.))
           PeriodicGrid1D(dx=1.0, nx=2)

        """

        if type(d) in [type(1), type(1.)] or not hasattr(d, '__len__'):
            n = int(n or 1)
        else:
            n = int(n or len(d))
            if n != len(d) and len(d) != 1:
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
        if hasattr(self, 'vertexCoords'):
            return self.vertexCoords
        else:
            return self._createVertices()

    def getExteriorFaces(self):
        """
        Return only the faces that have one neighboring cell.
        """
        return self.exteriorFaces
            
    def getInteriorFaces(self):
        """
        Return only the faces that have two neighboring cells.
        """
        return self.interiorFaces
        
    def getFaceCellIDs(self):
        return self.faceCellIDs

    def _getFaces(self):
        return numerix.arange(self.numberOfFaces)

    def _getCellsByID(self, ids = None):
        if ids is None:
            ids = range(self.numberOfCells) 
        return numerix.array([Cell(self, ID) for ID in ids])
    
    def _getMaxFacesPerCell(self):
        if self.getNumberOfCells() == 0:
            return 0
        else:
            return len(self.cellFaceIDs[...,0])

    """Geometry methods"""

    def _calcGeometry(self):
        self._calcFaceCenters()
        _CommonMesh._calcGeometry(self)
        self._calcCellNormals()
        
    """calc geometry methods"""

    def _calcFaceAreas(self):
        faceVertexIDs = MA.filled(self.faceVertexIDs, -1)
        substitute = numerix.repeat(faceVertexIDs[numerix.newaxis, 0], 
                                    faceVertexIDs.shape[0], axis=0)
        faceVertexIDs = numerix.where(MA.getmaskarray(self.faceVertexIDs), substitute, faceVertexIDs)
        faceVertexCoords = numerix.take(self.vertexCoords, faceVertexIDs, axis=1)
        faceOrigins = numerix.repeat(faceVertexCoords[:,0], faceVertexIDs.shape[0], axis=0)
        faceOrigins = numerix.reshape(faceOrigins, MA.shape(faceVertexCoords))
        faceVertexCoords = faceVertexCoords - faceOrigins
        left = range(faceVertexIDs.shape[0])
        right = left[1:] + [left[0]]
        cross = numerix.sum(numerix.cross(faceVertexCoords, numerix.take(faceVertexCoords, right, 1), axis=0), 1)
        self.faceAreas = numerix.sqrtDot(cross, cross) / 2.

    def _calcFaceCenters(self):
        faceVertexIDs = MA.filled(self.faceVertexIDs, 0)

        faceVertexCoords = numerix.take(self.vertexCoords, faceVertexIDs, axis=1)


        if MA.getmask(self.faceVertexIDs) is False:
            faceVertexCoordsMask = numerix.zeros(numerix.shape(faceVertexCoords))
        else:
            faceVertexCoordsMask = numerix.repeat(MA.getmaskarray(self.faceVertexIDs)[numerix.newaxis,...], self.dim, axis=0)

            
        faceVertexCoords = MA.array(data=faceVertexCoords, mask=faceVertexCoordsMask)

        self.faceCenters = MA.filled(MA.average(faceVertexCoords, axis=1))

        

    def _calcFaceNormals(self):
        faceVertexIDs = MA.filled(self.faceVertexIDs, 0)
        faceVertexCoords = numerix.take(self.vertexCoords, faceVertexIDs, axis=1)
        t1 = faceVertexCoords[:,1,:] - faceVertexCoords[:,0,:]
        t2 = faceVertexCoords[:,2,:] - faceVertexCoords[:,1,:]
        norm = numerix.cross(t1, t2, axis=0)
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
        mag = numerix.sqrt(numerix.sum(diff**2))
        self.faceCellToCellNormals = diff / numerix.resize(mag, (self.dim, len(mag)))

        orientation = 1 - 2 * (numerix.dot(self.faceNormals, self.faceCellToCellNormals) < 0)
        self.faceCellToCellNormals *= orientation

    def _calcOrientedFaceNormals(self):
        self.orientedFaceNormals = self.faceNormals
        
    def _calcCellVolumes(self):
        tmp = self.faceCenters[0] * self.faceAreas * self.faceNormals[0]
        tmp = numerix.take(tmp, self.cellFaceIDs) * self.cellToFaceOrientations
        self.cellVolumes = MA.filled(MA.sum(tmp, 0))

    def _calcCellCenters(self):
        tmp = numerix.take(self.faceCenters, self.cellFaceIDs, axis=1)
        self.cellCenters = MA.filled(MA.average(tmp, 1))
        
    def _calcFaceToCellDistances(self):
        tmp = MA.repeat(self.faceCenters[...,numerix.NewAxis,:], 2, 1)
        # array -= masked_array screws up masking for on numpy 1.1

        tmp = tmp - numerix.take(self.cellCenters, self.faceCellIDs, axis=1)
        self.cellToFaceDistanceVectors = tmp
        self.faceToCellDistances = MA.sqrt(MA.sum(tmp * tmp,0))

    def _calcCellDistances(self):
        tmp = numerix.take(self.cellCenters, self.faceCellIDs, axis=1)
        tmp = tmp[...,1,:] - tmp[...,0,:]
        tmp = MA.filled(MA.where(MA.getmaskarray(tmp), self.cellToFaceDistanceVectors[:,0], tmp))
        self.cellDistanceVectors = tmp
        self.cellDistances = MA.filled(MA.sqrt(MA.sum(tmp * tmp, 0)))

    def _calcFaceToCellDistanceRatio(self):
        dAP = self._getCellDistances()
        dFP = self._getFaceToCellDistances()[0]
        
        self.faceToCellDistanceRatio = MA.filled(dFP / dAP)

    def _calcAreaProjections(self):
        self.areaProjections = self._getFaceNormals() * self._getFaceAreas()
        
    def _calcOrientedAreaProjections(self):
        self.orientedAreaProjections = self.areaProjections

    def _calcFaceTangents(self):
        faceVertexCoord = numerix.array(numerix.take(self.vertexCoords, 
                                                     self.faceVertexIDs[0], 
                                                     axis=1))
        tmp = self.faceCenters - faceVertexCoord
        self.faceTangents1 = tmp / numerix.sqrtDot(tmp, tmp)
        tmp = numerix.cross(self.faceTangents1, self.faceNormals, axis=0)
        self.faceTangents2 = tmp / numerix.sqrtDot(tmp, tmp)
        
    def _calcCellToCellDistances(self):
        self.cellToCellDistances = numerix.take(self.cellDistances, self._getCellFaceIDs())

    def _calcCellNormals(self):
        cellNormals = numerix.take(self._getFaceNormals(), self._getCellFaceIDs(), axis=1)
        cellFaceCellIDs = numerix.take(self.faceCellIDs[0], self.cellFaceIDs)
        cellIDs = numerix.repeat(numerix.arange(self.getNumberOfCells())[numerix.newaxis,...], 
                                 self._getMaxFacesPerCell(), 
                                 axis=0)
        direction = (cellFaceCellIDs == cellIDs) * 2 - 1
        if self._getMaxFacesPerCell() > 0:
            self.cellNormals =  direction[numerix.newaxis, ...] * cellNormals
        else:
            self.cellNormals = cellNormals
                         
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
        cellVertexIDs = MA.sort(cellVertexIDs, axis=0, fill_value=-1)

        cellVertexIDs = MA.sort(MA.concatenate((cellVertexIDs[-1, numerix.newaxis], 
                                                MA.masked_where(cellVertexIDs[:-1] 
                                                                == cellVertexIDs[1:], 
                                                                cellVertexIDs[:-1]))), 
                                axis=0, fill_value=-1)
        
        ## resize the array to remove extra masked values
        length = min(numerix.sum(MA.getmaskarray(cellVertexIDs), axis=0))
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
            >>> print numerix.allequal(externalFaces, 
            ...                        numerix.nonzero(mesh.getExteriorFaces()))
            1

            >>> internalFaces = numerix.array((3,))
            >>> print numerix.allequal(internalFaces, 
            ...                        numerix.nonzero(mesh.getInteriorFaces()))
            1

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values((( 0,  0,  0, 0,  0,  0,  1,  1,  1,  1),
            ...                                 (-1, -1, -1, 1, -1, -1, -1, -1, -1, -1)), -1)
            >>> numerix.allequal(faceCellIds, mesh.getFaceCellIDs())
            1
            
            >>> dxdy = dx * dy
            >>> dxdz = dx * dz
            >>> dydz = dy * dz
            >>> faceAreas = numerix.array((dxdy, dxdy, dydz, dydz, dxdz, dxdz,
            ...                            dxdy/2., dxdy/2., dxdz, numerix.sqrt(dx**2 + dy**2) * dz))
            >>> numerix.allclose(faceAreas, mesh._getFaceAreas(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = numerix.take(vertices, MA.filled(faces, 0), axis=1)
            >>> faceCenters = faceCoords[...,0,:] + faceCoords[...,1,:] + faceCoords[...,2,:] + faceCoords[...,3,:]
            >>> numVex = numerix.array((4., 4., 4., 4., 4., 4., 3., 3., 4., 4.))
            >>> faceCenters /= numVex
            >>> numerix.allclose(faceCenters, mesh.getFaceCenters(), atol = 1e-10, rtol = 1e-10)
            1

            >>> faceNormals = numerix.array((( 0., 0., -1., 1.,  0., 0.,  0., 0.,  0., dy / numerix.sqrt(dy**2 + dx**2)),
            ...                              ( 0., 0.,  0., 0., -1., 1.,  0., 0., -1., dx / numerix.sqrt(dy**2 + dx**2)),
            ...                              (-1., 1.,  0., 0.,  0., 0., -1., 1.,  0., 0.)))
            >>> numerix.allclose(faceNormals, mesh._getFaceNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = MA.masked_values(((1, -1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, -2)), -2)
            >>> numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            1
                                             
            >>> cellVolumes = numerix.array((dx*dy*dz, dx*dy*dz / 2.))
            >>> numerix.allclose(cellVolumes, mesh.getCellVolumes(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellCenters = numerix.array(((dx/2., dx+dx/3.),
            ...                              (dy/2.,    dy/3.),
            ...                              (dz/2.,    dz/2.)))
            >>> numerix.allclose(cellCenters, mesh.getCellCenters(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> d1 = numerix.sqrt((dx / 3.)**2 + (dy / 6.)**2)
            >>> d2 = numerix.sqrt((dx / 6.)**2 + (dy / 3.)**2)
            >>> d3 = numerix.sqrt((dx / 6.)**2 + (dy / 6.)**2)
            >>> d4 = numerix.sqrt((5 * dx / 6.)**2 + (dy / 6.)**2)
            >>> faceToCellDistances = MA.masked_values(((dz / 2., dz / 2., dx / 2., dx / 2., dy / 2., dy / 2., dz / 2., dz / 2., d2, d3),
            ...                                         (     -1,      -1,      -1,      d1,      -1,      -1,      -1,      -1, -1, -1)), -1)
            >>> numerix.allclose(faceToCellDistances, mesh._getFaceToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> cellDistances = numerix.array((dz / 2., dz / 2., dx / 2.,
            ...                                d4,
            ...                                dy / 2., dy / 2., dz / 2., dz / 2.,
            ...                                d2,
            ...                                d3))
            >>> numerix.allclose(cellDistances, mesh._getCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> numerix.allclose(faceToCellDistanceRatios, mesh._getFaceToCellDistanceRatio(), atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas
            >>> numerix.allclose(areaProjections, mesh._getAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> v1 = numerix.take(vertices, numerix.array(faces[0]), axis=1)
            >>> tmp = faceCenters - v1
            >>> tangents1 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> numerix.allclose(tangents1, mesh._getFaceTangents1(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tmp = numerix.cross(tangents1, faceNormals, axis=0)
            >>> tangents2 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> numerix.allclose(tangents2, mesh._getFaceTangents2(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((-1,  0),
            ...                                   (-1, -1),
            ...                                   (-1, -1), 
            ...                                   ( 1, -1), 
            ...                                   (-1, -1),
            ...                                   (-1, -1)), -1)
            >>> numerix.allequal(cellToCellIDs, mesh._getCellToCellIDs())
            1

            >>> cellToCellDistances = MA.masked_values(((dz / 2., d4),
            ...                                         (dz / 2., -1),
            ...                                         (dx / 2., -1),
            ...                                         (     d4, -1),
            ...                                         (dy / 2., -1),
            ...                                         (dy / 2., -1)), -1)
            >>> numerix.allclose(cellToCellDistances, mesh._getCellToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1

            >>> interiorCellIDs = numerix.array(())
            >>> numerix.allequal(interiorCellIDs, mesh._getInteriorCellIDs())
            1

            >>> exteriorCellIDs = numerix.array((0, 1))
            >>> numerix.allequal(exteriorCellIDs, mesh._getExteriorCellIDs())
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
            >>> numerix.allclose(cellNormals, mesh._getCellNormals(), atol = 1e-10, rtol = 1e-10)
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
            >>> numerix.allclose(cellAreaProjections, mesh._getCellAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellVertexIDs = MA.masked_values(((7, 6, 5, 4, 3, 2, 1, 0), (9, 8, 6, 5, 2, 1, -1000, -1000)), -1000)
            >>> cellVertexIDs = MA.masked_values(((7, 9),
            ...                                   (6, 8),
            ...                                   (5, 6),
            ...                                   (4, 5),
            ...                                   (3, 2),
            ...                                   (2, 1),
            ...                                   (1, -1000),
            ...                                   (0, -1000)), -1000)
            >>> numerix.allclose(cellVertexIDs, mesh._getCellVertexIDs())
            1


            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')
            >>> unpickledMesh = dump.read(filename, f)

            >>> numerix.allequal(mesh.getCellCenters(), unpickledMesh.getCellCenters())
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
            >>> numerix.allclose(bigMesh.getCellVolumes(), volumes)
            1
            
        """

                      
### test test test    

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
