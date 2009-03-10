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

from fipy.tools import numerix

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
    
    def _calcTopology(self):
	self.dim = len(self.vertices[0].getCoordinates())
	self.numberOfFaces = len(self.faces)
	self.numberOfCells = len(self.cells)
	self._calcCellFaceIDs()
	
	fipy.meshes.common.mesh.Mesh._calcTopology(self)

	self._calcFaceOrientations()

    """calc topology methods"""
	
    def _calcCellFaceIDs(self):
	cells = self.getCells()
	for cell in cells:
	    cell._calcCenter()
	self.cellFaceIDs = ()
	self.cellFaceIDIndices = ()
	for i in range(len(cells)):
	    cell = cells[i]
	    ids = cell.getFaceIDs()
	    self.cellFaceIDs += ids
	self.cellFaceIDs = numerix.array(self.cellFaceIDs)
	self.cellFaceIDs = numerix.reshape(self.cellFaceIDs, (len(cells), self._getMaxFacesPerCell()))
	
    def _calcInteriorAndExteriorFaceIDs(self):
	self.exteriorFaceIDs = numerix.arange(len(self.interiorFaces),len(self.faces))
	self.interiorFaceIDs = numerix.arange(len(self.interiorFaces))
	
    def _calcExteriorCellIDs(self):
	# we want a list of unique ids
	self.exteriorCellIDs = list(sets.Set([face.getCellID() for face in self.getExteriorFaces()]))
	
    def _calcCellToFaceOrientations(self):
	cells = self.getCells()
	N = len(cells)
	M = self._getMaxFacesPerCell()
	self.cellToFaceOrientations = numerix.zeros((N,M,1))
	for i in range(N):
	    orientations = cells[i].getFaceOrientations()
	    orientations = numerix.reshape(orientations,(len(cells[i].getFaces()),))
	    for j in range(len(orientations)):
		self.cellToFaceOrientations[i,j,0] = orientations[j]
		
    def _calcAdjacentCellIDs(self):
	self.adjacentCellIDs = (numerix.array([face.getCellID(0) for face in self.faces]),
				numerix.array([face.getCellID(1) for face in self.faces]))

    def _calcCellToCellIDs(self):
	pass

    def _calcFaceOrientations(self):
	faces = self.getFaces()
	N = len(faces)
	orientations = numerix.zeros((N),'d')
	for i in range(N):
	    orientations[i] = faces[i]._getOrientation()
	self.faceOrientations = orientations

    """get topology methods"""
	
    def _getCellsByID(self, ids = None):
	if ids is None:
	    ids = range(self.numberOfCells) 
	return [self.cells[id] for id in ids]
	
    def _getFaces(self):
        return self.faces

    def _getInteriorFaces(self):
        """
        Return only the faces that have two neighboring cells.
        """
        return self.interiorFaces

    def getExteriorFaces(self):
        """
        Return only the faces that have one neighboring cell.
        """
        return self.faces[len(self._getInteriorFaces()):]
    
##    def removeCell(self,cell):
##        """Remove cell from Mesh.
##        """
##        rCellId = cell.getId()
##        for face in cell.getFaces():
##            face._removeBoundingCell(cell)
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
    
    def _calcFaceAreas(self):
	faces = self.getFaces()
	N = len(faces)
	self.faceAreas = numerix.zeros((N),'d')
	# get the units right
	self.faceAreas = self.faceAreas * faces[0].getArea()
	for i in range(N):
	    self.faceAreas[i] = faces[i].getArea()
	    
    def _calcCellVolumes(self):
	cells = self.getCells()
	N = len(cells)
	self.cellVolumes = numerix.zeros((N),'d')
	# get the units right
	self.cellVolumes = self.cellVolumes * cells[0].getVolume()	    
	for i in range(N):
	    self.cellVolumes[i] = cells[i].getVolume()	    
	    
    def _calcCellCenters(self):
	cells = self.getCells()
	N = len(cells)
	self.cellCenters = numerix.zeros((N,self.dim),'d')
	# get the units right
	self.cellCenters = self.cellCenters * cells[0].getCenter()
	for i in range(N):
	    self.cellCenters[i] = cells[i].getCenter()	    
	
    def _calcCellDistances(self):
	faces = self.getFaces()
	N = len(faces)
	self.cellDistances = numerix.zeros((N),'d')
	# get the units right
	self.cellDistances = self.cellDistances * faces[0].getCellDistance()
	for i in range(N):
	    self.cellDistances[i] = faces[i].getCellDistance()
	
    def _calcFaceToCellDistances(self):
	faces = self.getFaces()
	N = len(faces)
	self.faceToCellDistances = numerix.zeros((N),'d')
	# get the units right
	self.faceToCellDistances = self.faceToCellDistances * faces[0]._getFaceToCellDistance()
	for i in range(N):
	    self.faceToCellDistances[i] = faces[i]._getFaceToCellDistance()

    def _calcFaceNormals(self):
	faces = self.getFaces()
	N = len(faces)
	dim = len(faces[0].getCenter())
	self.faceNormals = numerix.zeros((N,dim),'d')
	# get the units right
	self.faceNormals = self.faceNormals * faces[0]._calcNormal()
	for i in range(N):
	    self.faceNormals[i] = faces[i]._calcNormal()
	    
    def _calcOrientedFaceNormals(self):
        self.orientedFaceNormals = self._getFaceNormals() * self.getFaceOrientations()[:,numerix.NewAxis]
	    
    def _calcAreaProjections(self):
        self.areaProjections = self._getFaceNormals() * self._getFaceAreas()[:,numerix.NewAxis] 
	
    def _calcOrientedAreaProjections(self):
	self.orientedAreaProjections = self._getAreaProjections() * self.getFaceOrientations()[:,numerix.NewAxis]
	
    def _calcFaceTangents(self):
	faces = self.getFaces()
	N = len(faces)
	dim = len(faces[0].getCenter())
	self.faceTangents1 = numerix.zeros((N,dim),'d')
	self.faceTangents2 = numerix.zeros((N,dim),'d')
	# get the units right
	self.faceTangents1 = self.faceTangents1 * faces[0]._calcTangent1()
	self.faceTangents2 = self.faceTangents2 * faces[0]._calcTangent2()
	for i in range(N):
	    self.faceTangents1[i] = faces[i]._calcTangent1()
	    self.faceTangents2[i] = faces[i]._calcTangent2()

    def _calcFaceToCellDistanceRatio(self):
	dAP = self._getCellDistances()
	dFP = self._getFaceToCellDistances()
	self.faceToCellDistanceRatio = dFP / dAP

    def _calcCellToCellDistances(self):
        from fipy.tools import numerix
	self.cellToCellDistances = numerix.take(self._getCellDistances(), self._getCellFaceIDs())


