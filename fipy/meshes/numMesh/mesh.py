#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
 #                                    created: 11/10/03 {2:44:42 PM} 
 #                                last update: 5/13/04 {1:37:41 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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
"""



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

    """calc geometry methods"""

    def calcFaceAreas(self):
        faceVertexIDs = MA.filled(self.faceVertexIDs, -1)
        substitute = Numeric.reshape(Numeric.repeat(faceVertexIDs[:,0],len(faceVertexIDs[0])), Numeric.shape(faceVertexIDs))
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
        self.cellToCellDistances = MAtake(self.getCellDistances(), self.getCellFaceIDs())

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
        
                      
    
