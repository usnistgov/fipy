#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "mesh.py"
 #                                    created: 11/10/03 {2:44:42 PM} 
 #                                last update: 12/10/03 {11:59:14 AM} 
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

"""Generic mesh class

    Meshes contain cells, faces, and vertices.
"""

import Numeric

class Mesh:
    def __init__(self, cells, faces, interiorFaces, vertices):
        self.cells = cells
        self.faces = faces
        self.vertices = vertices
        self.interiorFaces = interiorFaces
        self.dim = len(self.vertices[0].getCoordinates())
	
	self.calcCellFaceIDs(cells)
	self.calcCellCenters(cells)
	self.calcFaceTangents(faces)
	
	self.calcAdjacentCellIDs()

    def getCells(self,func = None):
	"""Return Cells of Mesh.
	"""
        if func == None:
            return self.cells
        else:        
            returnCells = ()
            for cell in self.cells:
                if func(cell.getCenter()):
                    returnCells += (cell,)
            return returnCells

    def getFaces(self):
	"""Return Faces of Mesh.
	"""
        return self.faces

    def getInteriorFaces(self):
	"""Return Faces of Mesh that are not on a Mesh boundary.
	"""
        return self.interiorFaces

    def getExteriorFaces(self):
        """Return all exterior faces
        """
        exteriorFaces = ()
        for face in self.faces:
            if face not in self.getInteriorFaces():
                exteriorFaces += (face,)
        return exteriorFaces
    
    def makeGridData(self,array):
	"""Return array data mapped onto cell geometry of Mesh.
	"""
        pass

    def getDim(self):
        return self.dim 

##    def removeCell(self,cell):
##        """Remove cell from Mesh.
##        """
##        rCellId = cell.getId()
##        for face in cell.getFaces():
##            face.removeBoundingCell(cell)
##            if face.cells == ():
##                self.removeFace(face)
##        cell = self.cells[-1]
##        cell.setId(rCellId)
##        self.cells = self.cells[:-1]

##    def removeFace(self,face):
##        """Remove face from Mesh.
##        """
##        rFaceId = face.getId()
##        face = self.faces[-1]
##        print len(self.faces)
##        face.setId(rFaceId)
##        self.faces = self.faces[:-1]
##        print len(self.faces)
##        raw_input()

    def getPhysicalShape(self):
	"""Return physical dimensions of Mesh.
	"""
        pass

    def getFaceAreas(self):
	pass

    def getAdjacentCellIDs(self):
	return (self.cellId1, self.cellId2)
	
    def calcAdjacentCellIDs(self):
	N = len(self.faces)
	self.cellId1 = Numeric.zeros(N)
	self.cellId2 = Numeric.zeros(N)
	for i in range(N):
	    self.cellId1[i] = self.faces[i].getCellId(0)
	    self.cellId2[i] = self.faces[i].getCellId(1)

    def getFaceOrientations(self):
        return self.faceOrientations
    
    def calcFaceOrientations(self,faces):
        N = len(faces)
        orientations = Numeric.zeros((N),'d')
        for i in range(N):
            orientations[i] = faces[i].getOrientation()
        self.faceOrientations = orientations

    def getMaxFacesPerCell(self):
        pass
    
    def getCellFaceOrientations(self):
        return self.cellFaceOrientations

    def calcCellFaceOrientations(self,cells):
        N = len(cells)
        M = self.getMaxFacesPerCell()
        self.cellFaceOrientations = Numeric.zeros((N,M,1))
        for i in range(N):
            orientations = cells[i].getFaceOrientations()
            orientations = Numeric.reshape(orientations,(len(cells[i].getFaces()),))
            for j in range(len(orientations)):
                self.cellFaceOrientations[i,j,0] = orientations[j]

    def getCellFaceIDs(self):
        return self.cellFaceIDs

    def calcCellFaceIDs(self,cells):
	for cell in cells:
	    cell.calcFaceIDs()
        self.cellFaceIDs = ()
        self.cellFaceIDIndices = ()
        for i in range(len(cells)):
            cell = cells[i]
            ids = cell.getFaceIDs()
            self.cellFaceIDs += ids
        
    def getFaceAreas(self):
	return self.faceAreas
	
    def calcFaceAreas(self,faces):
	N = len(faces)
	self.faceAreas = Numeric.zeros((N),'d')
	for i in range(N):
	    self.faceAreas[i] = faces[i].getArea()
	    
    def getCellVolumes(self):
	return self.cellVolumes
	
    def calcCellVolumes(self,cells):
	N = len(cells)
	self.cellVolumes = Numeric.zeros((N),'d')
	for i in range(N):
	    self.cellVolumes[i] = cells[i].getVolume()	    
	    
    def getCellCenters(self):
	return self.cellCenters
	
    def calcCellCenters(self,cells):
	N = len(cells)
	self.cellCenters = Numeric.zeros((N,self.dim),'d')
	for i in range(N):
	    self.cellCenters[i] = cells[i].getCenter()	    
	
    def getCellDistances(self):
	return self.cellDistances
	
    def calcCellDistances(self,faces):
	N = len(faces)
	self.cellDistances = Numeric.zeros((N),'d')
	for i in range(N):
	    self.cellDistances[i] = faces[i].getCellDistance()
	
    def getFaceToCellDistances(self):
	return self.faceToCellDistances
	
    def calcFaceToCellDistances(self,faces):
	N = len(faces)
	self.faceToCellDistances = Numeric.zeros((N),'d')
	for i in range(N):
	    self.faceToCellDistances[i] = faces[i].getFaceToCellDistance()

    def getFaceNormals(self):
	return self.faceNormals

    def calcFaceNormals(self, faces):
	N = len(faces)
	dim = len(faces[0].getCenter())
	self.faceNormals = Numeric.zeros((N,dim),'d')
	for i in range(N):
	    self.faceNormals[i] = faces[i].calcNormal()
	    
    def getOrientedFaceNormals(self):
	return self.orientedFaceNormals
	
    def calcOrientedFaceNormals(self):
	self.orientedFaceNormals = self.getFaceNormals().copy()
	self.orientedFaceNormals *= Numeric.reshape(self.getFaceOrientations(),(len(self.orientedFaceNormals),1))
	    
    def getFaceAreas(self):
	return self.faceAreas
	
    def calcFaceAreas(self, faces):
	N = len(faces)
	self.faceAreas = Numeric.zeros(N,'d')
	for i in range(N):
	    self.faceAreas[i] = faces[i].getArea()
	    
    def getAreaProjections(self):
	return self.areaProjections
	
    def calcAreaProjections(self):
	N = len(self.faceNormals)
	self.areaProjections = self.faceNormals * Numeric.reshape(self.faceAreas,(N,1))

    def getOrientedAreaProjections(self):
	return self.orientedAreaProjections
	
    def calcOrientedAreaProjections(self):
	N = len(self.faceNormals)
	self.orientedAreaProjections = self.areaProjections * Numeric.reshape(self.getFaceOrientations(),(N,1))
	
    def calcFaceTangents(self, faces):
	N = len(faces)
	dim = len(faces[0].getCenter())
	self.faceTangents1 = Numeric.zeros((N,dim),'d')
	self.faceTangents2 = Numeric.zeros((N,dim),'d')
	for i in range(N):
	    self.faceTangents1[i] = faces[i].calcTangent1()
	    self.faceTangents2[i] = faces[i].calcTangent2()

    def getFaceTangents1(self):
	return self.faceTangents1

    def getFaceTangents2(self):
	return self.faceTangents2
	
    def getFaceToCellDistanceRatio(self):
	return self.faceToCellDistanceRatio
	
    def calcFaceToCellDistanceRatio(self):
	dAP = self.getCellDistances()
	dFP = self.getFaceToCellDistances()
	self.faceToCellDistanceRatio = dFP / dAP

    def refreshFaces(self,faces):
	self.calcFaceOrientations(faces)
	self.calcFaceAreas(faces)
	self.calcCellDistances(faces)
	self.calcFaceToCellDistances(faces)
	self.calcFaceToCellDistanceRatio()
	self.calcFaceNormals(faces)
	self.calcOrientedFaceNormals()
	self.calcAreaProjections()
	self.calcOrientedAreaProjections()

        
