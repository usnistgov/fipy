#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
 #                                    created: 11/10/03 {2:44:42 PM} 
 #                                last update: 5/6/04 {11:28:45 AM} 
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

"""Generic mesh class defining implementation-agnostic behavior

    Meshes contain cells, faces, and vertices.
"""

## import sets

import Numeric

## from fipy.tools.dimensions.physicalField import PhysicalField



class Mesh:
    def __init__(self):
	self.calcTopology()
	self.calcGeometry()
    
    """topology methods"""
    
    def calcTopology(self):
	self.calcInteriorAndExteriorFaceIDs()
	self.calcInteriorAndExteriorCellIDs()
	self.calcCellToFaceOrientations()
	self.calcAdjacentCellIDs()
	self.calcCellToCellIDs()
       
    """get Topology methods"""

    def getCellFaceIDs(self):
        return self.cellFaceIDs

    def getExteriorFaceIDs(self):
	return self.exteriorFaceIDs
	
    def getExteriorFaces(self):
	pass

    def getInteriorFaceIDs(self):
	return self.interiorFaceIDs
	    
    def getInteriorFaces(self):
	pass
	
    def getExteriorCellIDs(self):
	return self.exteriorCellIDs

    def getInteriorCellIDs(self):
	return self.interiorCellIDs

    def getCellFaceOrientations(self):
        return self.cellToFaceOrientations

    def getNumberOfCells(self):
	return self.numberOfCells
	
    def getAdjacentCellIDs(self):
        return self.adjacentCellIDs

    def getDim(self):
        return self.dim

    def getFaces(self):
        pass

    def getCellsByID(self, ids = None):
	pass
	    
    def getCells(self, ids = None, filter = None, **args):
	"""Return Cells of Mesh."""
	cells = self.getCellsByID(ids)
	
        if filter is None:
            return cells
        else:        
	    return [cell for cell in cells if filter(cell, **args)]
	    
    def getMaxFacesPerCell(self):
	pass

    def getNumberOfFaces(self):
	return self.numberOfFaces

    def getCellToCellIDs(self):
        return self.cellToCellIDs
	
    """geometry methods"""
    
    def calcGeometry(self):
	self.calcFaceAreas()
	self.calcFaceNormals()
	self.calcOrientedFaceNormals()
	self.calcCellVolumes()
	self.calcCellCenters()
	self.calcFaceToCellDistances()
	self.calcCellDistances()        
	self.calcAreaProjections()
	self.calcOrientedAreaProjections()
	self.calcFaceTangents()
	self.calcFaceToCellDistanceRatio()
	self.calcFaceAspectRatios()
	self.calcCellToCellDistances()
##        self.setScale()
       
    """calc geometry methods"""
    
    def calcFaceAspectRatios(self):
	self.faceAspectRatios = self.getFaceAreas() / self.getCellDistances()

    """get geometry methods"""
        
    def getFaceAreas(self):
        return self.faceAreas

##     def getFaceCenters(self):
##         return self.faceCenters

    def getFaceNormals(self):
        return self.faceNormals
	
    def getCellVolumes(self):
	return self.cellVolumes

    def getCellCenters(self):
	return self.cellCenters

    def getFaceToCellDistances(self):
        return self.faceToCellDistances

    def getCellDistances(self):
        return self.cellDistances

    def getFaceToCellDistanceRatio(self):
        return self.faceToCellDistanceRatio

    def getOrientedAreaProjections(self):
        return self.orientedAreaProjections

    def getAreaProjections(self):
        return self.areaProjections

    def getOrientedFaceNormals(self):
        return self.orientedFaceNormals

    def getFaceTangents1(self):
        return self.faceTangents1

    def getFaceTangents2(self):
        return self.faceTangents2
	
    def getFaceAspectRatios(self):
	return self.faceAspectRatios
    
    def getCellToCellDistances(self):
	return self.cellToCellDistances
	    
    """scaling"""

    def setScale(self, value = 1.):
        self.scale = 1.
    
##    def setScale(self):
##	self.scale = PhysicalField(value = 1.)
##        self.faceAreas = self.faceAreas * self.scale**2
##        self.faceCenters = self.faceCenters * self.scale
##        self.cellVolumes = self.cellVolumes * self.scale**3
##        self.cellCenters = self.cellCenters * self.scale
##        self.faceToCellDistances = self.faceToCellDistances * self.scale
##        self.cellDistances = self.cellDistances * self.scale
##        self.areaProjections = self.areaProjections * self.scale**2

    """point to cell distances"""
    
    def getPointToCellDistances(self, point):
	tmp = self.getCellCenters() - Numeric.array(point)
	
	import fipy.tools.array
	return fipy.tools.array.sqrtDot(tmp, tmp)

    def getNearestCell(self, point):
        d = self.getPointToCellDistances(point)
        i = Numeric.argsort(d)
	return self.getCellsByID([i[0]])[0]
    

## pickling

##    self.__getinitargs__(self):
##        return (self.vertexCoords, self.faceVertexIDs, self.cellFaceIDs)
    

## ##     def __getstate__(self):
## ##         dict = {
## ##             'vertexCoords' : self.vertexCoords,            
## ##             'faceVertexIDs' : self.faceVertexIDs,
## ##             'cellFaceIDs' : self.cellFaceIDs }
## ##         return dict
## ## 
## ##     def __setstate__(self, dict):
## ##         self.__init__(dict['vertexCoords'], dict['faceVertexIDs'], dict['cellFaceIDs'])
        
                      
    
