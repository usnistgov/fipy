#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
 #                                    created: 11/10/03 {2:44:42 PM} 
 #                                last update: 9/3/04 {10:39:00 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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

"""Generic mesh class

    Meshes contain cells, faces, and vertices.
"""

import sets

import Numeric

import fipy.meshes.common.mesh

from fipy.tools.dimensions.physicalField import PhysicalField

class Mesh(fipy.meshes.common.mesh.Mesh):
    def __init__(self, cells, faces, interiorFaces, vertices):
        self.cells = cells
        self.faces = faces
        self.vertices = vertices
        self.interiorFaces = interiorFaces
        self.dim = len(self.vertices[0].getCoordinates())
	if not self.__dict__.has_key("scale"):
	    self.scale = 1
	
	fipy.meshes.common.mesh.Mesh.__init__(self)

    """topology methods"""
    
    def calcTopology(self):
	self.dim = len(self.vertices[0].getCoordinates())
	self.numberOfFaces = len(self.faces)
	self.numberOfCells = len(self.cells)
	self.calcCellFaceIDs()
	
	fipy.meshes.common.mesh.Mesh.calcTopology(self)

	self.calcFaceOrientations()

    """calc topology methods"""
	
    def calcCellFaceIDs(self):
	cells = self.getCells()
	for cell in cells:
	    cell.calcFaceIDs()
	self.cellFaceIDs = ()
	self.cellFaceIDIndices = ()
	for i in range(len(cells)):
	    cell = cells[i]
	    ids = cell.getFaceIDs()
	    self.cellFaceIDs += ids
	self.cellFaceIDs = Numeric.array(self.cellFaceIDs)
	self.cellFaceIDs = Numeric.reshape(self.cellFaceIDs, (len(cells), self.getMaxFacesPerCell()))
	
    def calcInteriorAndExteriorFaceIDs(self):
	self.exteriorFaceIDs = Numeric.arange(len(self.interiorFaces),len(self.faces))
	self.interiorFaceIDs = Numeric.arange(len(self.interiorFaces))
	
    def calcExteriorCellIDs(self):
	# we want a list of unique ids
	self.exteriorCellIDs = list(sets.Set([face.getCellID() for face in self.getExteriorFaces()]))
	
    def calcCellToFaceOrientations(self):
	cells = self.getCells()
	N = len(cells)
	M = self.getMaxFacesPerCell()
	self.cellToFaceOrientations = Numeric.zeros((N,M,1))
	for i in range(N):
	    orientations = cells[i].getFaceOrientations()
	    orientations = Numeric.reshape(orientations,(len(cells[i].getFaces()),))
	    for j in range(len(orientations)):
		self.cellToFaceOrientations[i,j,0] = orientations[j]
		
    def calcAdjacentCellIDs(self):
	self.adjacentCellIDs = (Numeric.array([face.getCellID(0) for face in self.faces]),
				Numeric.array([face.getCellID(1) for face in self.faces]))

    def calcCellToCellIDs(self):
	pass
## 	self.cellToCellIDs = MAtake(self.faceCellIDs, self.cellFaceIDs)
## 	self.cellToCellIDs = MA.where(self.cellToFaceOrientations == 1, self.cellToCellIDs[:,:,1], self.cellToCellIDs[:,:,0])

    def calcFaceOrientations(self):
	faces = self.getFaces()
	N = len(faces)
	orientations = Numeric.zeros((N),'d')
	for i in range(N):
	    orientations[i] = faces[i].getOrientation()
	self.faceOrientations = orientations

    """get topology methods"""
	
    def getCellsByID(self, ids = None):
	if ids is None:
	    ids = range(self.numberOfCells) 
	return [self.cells[id] for id in ids]
	
    def getFaces(self):
        return self.faces

    def getInteriorFaces(self):
	"""Return Faces of Mesh that are not on a Mesh boundary."""
        return self.interiorFaces

    def getExteriorFaces(self):
        """Return all Faces of Mesh that are on a Mesh boundary"""
        return self.faces[len(self.getInteriorFaces()):]
    
    def makeGridData(self,array):
	"""Return array data mapped onto cell geometry of Mesh."""
        pass

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

    def getFaceOrientations(self):
        return self.faceOrientations
    
    """geometry methods"""
    
    def getPhysicalShape(self):
	"""Return physical dimensions of Mesh.
	"""
	pass

    def getScale(self):
	return self.scale

    def setScale(self, scale):
	self.scale = PhysicalField(value = scale)
	
    """calc geometry methods"""
    
    def calcFaceAreas(self):
	faces = self.getFaces()
	N = len(faces)
	self.faceAreas = Numeric.zeros((N),'d')
	# get the units right
	self.faceAreas = self.faceAreas * faces[0].getArea()
	for i in range(N):
	    self.faceAreas[i] = faces[i].getArea()
	    
    def calcCellVolumes(self):
	cells = self.getCells()
	N = len(cells)
	self.cellVolumes = Numeric.zeros((N),'d')
	# get the units right
	self.cellVolumes = self.cellVolumes * cells[0].getVolume()	    
	for i in range(N):
	    self.cellVolumes[i] = cells[i].getVolume()	    
	    
    def calcCellCenters(self):
	cells = self.getCells()
	N = len(cells)
	self.cellCenters = Numeric.zeros((N,self.dim),'d')
	# get the units right
	self.cellCenters = self.cellCenters * cells[0].getCenter()
	for i in range(N):
	    self.cellCenters[i] = cells[i].getCenter()	    
	
    def calcCellDistances(self):
	faces = self.getFaces()
	N = len(faces)
	self.cellDistances = Numeric.zeros((N),'d')
	# get the units right
	self.cellDistances = self.cellDistances * faces[0].getCellDistance()
	for i in range(N):
	    self.cellDistances[i] = faces[i].getCellDistance()
	
    def calcFaceToCellDistances(self):
	faces = self.getFaces()
	N = len(faces)
	self.faceToCellDistances = Numeric.zeros((N),'d')
	# get the units right
	self.faceToCellDistances = self.faceToCellDistances * faces[0].getFaceToCellDistance()
	for i in range(N):
	    self.faceToCellDistances[i] = faces[i].getFaceToCellDistance()

    def calcFaceNormals(self):
	faces = self.getFaces()
	N = len(faces)
	dim = len(faces[0].getCenter())
	self.faceNormals = Numeric.zeros((N,dim),'d')
	# get the units right
	self.faceNormals = self.faceNormals * faces[0].calcNormal()
	for i in range(N):
	    self.faceNormals[i] = faces[i].calcNormal()
	    
    def calcOrientedFaceNormals(self):
        self.orientedFaceNormals = self.getFaceNormals() * self.getFaceOrientations()[:,Numeric.NewAxis]
	    
    def calcAreaProjections(self):
        self.areaProjections = self.getFaceNormals() * self.getFaceAreas()[:,Numeric.NewAxis] 
	
    def calcOrientedAreaProjections(self):
	self.orientedAreaProjections = self.getAreaProjections() * self.getFaceOrientations()[:,Numeric.NewAxis]
	
    def calcFaceTangents(self):
	faces = self.getFaces()
	N = len(faces)
	dim = len(faces[0].getCenter())
	self.faceTangents1 = Numeric.zeros((N,dim),'d')
	self.faceTangents2 = Numeric.zeros((N,dim),'d')
	# get the units right
	self.faceTangents1 = self.faceTangents1 * faces[0].calcTangent1()
	self.faceTangents2 = self.faceTangents2 * faces[0].calcTangent2()
	for i in range(N):
	    self.faceTangents1[i] = faces[i].calcTangent1()
	    self.faceTangents2[i] = faces[i].calcTangent2()

    def calcFaceToCellDistanceRatio(self):
	dAP = self.getCellDistances()
	dFP = self.getFaceToCellDistances()
	self.faceToCellDistanceRatio = dFP / dAP

    def calcCellToCellDistances(self):
	self.cellToCellDistances = fipy.tools.array.take(self.getCellDistances(), self.getCellFaceIDs())


