#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "face.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 12/1/03 {3:55:22 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
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
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

"""Face within a Mesh
"""

import tools
import Numeric

class Face:
    """Face within a Mesh

    Faces are bounded by Vertices. Faces separate Cells.
    """
    
    def __init__(self, vertices, id):
	"""Face is initialized by Mesh
	
	Arguments:
	    
	    'vertices' -- the 'Vertex' points that bound the 'Face'
	    
	    'id' -- a unique identifier
	"""
        self.vertices = vertices
        self.cells = ()
        self.cellsId = ()
	self.id = id
	self.center = self.calcCenter()
        self.area = self.calcArea()
	self.setNormal()
            
    def addBoundingCell(self, cell, orientation):
	"""Add cell to the list of Cells which lie on either side of this Face.
	"""
        self.cells += (cell,)
        self.cellsId += (cell.getId(),)
	self.orientation = -orientation
	self.setCellDistance()
	self.setFaceToCellDistances()

        
    def getCells(self):
	"""Return the Cells which lie on either side of this Face.
	"""
        return self.cells
		
    def getId(self):
	"""Return the id of this Face.
	"""
	return self.id

    def getCellId(self, index = 0):
	"""Return the id of the specified Cell on one side of this Face.
	"""
        if index == 1 and len(self.cellsId) == 1:
            index = 0
        return self.cellsId[index]

##    def getCellId(self,index=0):
##        return self.cells[index].getId()
	
    def getCenter(self):
	"""Return the coordinates of the Face center.
	"""
	return self.center
    
    def calcCenter(self):
	"""Calculate the coordinates of the Face center.
	
	Cell center is the average of the bounding Vertex centers.
	"""
	ctr = self.vertices[0].getCoordinates().copy()
	for vertex in self.vertices[1:]:
	    ctr += vertex.getCoordinates()
	return ctr / float(len(self.vertices))

    def getArea(self):
	"""Return the area of the Face.
	"""
        return self.area
    
    def calcArea(self):
	"""Calculate the area of the Face.
	
	Area is the signed sum of the area of the triangles bounded by
	each polygon edge and the origin.
	"""
	a=0.
	p1 = self.vertices[0].getCoordinates().copy()
	for vertex in self.vertices[2:-1]:
	    p2=vertex.getCoordinates().copy()
	    a += tools.crossProd(p1,p2)
	    p1 = p2
	return abs(a/2.)
        
    def getNormal(self):
	"""Return the unit normal vector
	"""
	return self.normal
		
    def setNormal(self):
	"""Cache the unit normal vector.
	
	Called by Cell initializer after bounding Cells have been added to
	Faces.
	"""
        self.normal = self.calcNormal()
	
    def calcNormal(self):
	"""Calculate the unit normal vector, accounting for whether the Face
	points toward or away from the specified Cell.	
	
	Unit normal vector is calculated from cross-product of two
	tangent vectors.
    
	**Doesn't work if t1 and t2 are colinear!**
	"""
	t1 = self.vertices[1].getCoordinates() - self.vertices[0].getCoordinates()
	t2 = self.vertices[2].getCoordinates() - self.vertices[1].getCoordinates()
	norm = tools.crossProd(t1,t2)
	norm /= tools.sqrtDot(norm,norm)
## we calculate the orientation after we know the normal
##	norm *= self.orientation
	
	return norm
	
    def calcTangent1(self):
	norm = self.normal
	mag = Numeric.sqrt(norm[0]**2 + norm[1]**2)
	tan1 = Numeric.array((-norm[1],norm[0],0))
	return tan1/mag

    def calcTangent2(self):
	norm = self.normal
	mag = Numeric.sqrt(norm[0]**2 + norm[1]**2)
	tan2 = Numeric.array(norm[0] * norm[2], norm[1] * norm[2], -mag**2)
	return tan2/mag
	    
    def getCellDistance(self):
	"""Return the distance between adjacent cell centers.
	"""
        return self.cellDistance

    def setCellDistance(self):
	"""Assign the cached distance between adjacent cell centers.
	"""
        self.cellDistance = self.calcCellDistance()

    def calcCellDistance(self):
	"""Calculate the distance between adjacent Cell centers.
	
	If the Face is on a boundary and has only one bordering Cell,
	the distance is from the Cell center to the Face center.
	"""
        if(len(self.cells)==2):
            vec=self.cells[1].getCenter()-self.cells[0].getCenter()
        else:
            vec=self.center-self.cells[0].getCenter()        
        return tools.sqrtDot(vec,vec)

    def calcFaceToCellDistance(self, cell):
        vec=self.center-cell.getCenter()        
        return tools.sqrtDot(vec,vec)

    def setFaceToCellDistances(self):
        faceToCellDistances = ()
        for cell in self.cells:
            faceToCellDistances += (self.calcFaceToCellDistance(cell),)
        self.faceToCellDistances = faceToCellDistances

    def getFaceToCellDistance(self, cell = 'None'):
        if cell == self.cells[0] or cell == 'None':
            return self.faceToCellDistances[0]
        elif cell == self.cells[1]:
            return self.faceToCellDistances[1]
        else:
            return self.calcFaceToCellDistance(cell)

#     def __repr__(self):
# 	"""Textual representation of Face.
# 	"""
# 	return "<id = " + str(self.id) + ">"
# 	+ ", area = " + str(self.area()) + ", normal = " + str(self.normal()) + ", vertices = " + str(self.vertices) + ", centers = " + str(self.center) + ">\n"
    def __repr__(self):
	"""Textual representation of Face.
	"""
	return "<id = " + str(self.id) + ", area = " + str(self.area) + ", normal = " + str(self.normal) + ", vertices = " + str(self.vertices) + ", center = " + str(self.center) + ">\n"

    def removeBoundingCell(self,cell):
	"""Remove cell from the list of bounding Cells.
	
	Called by the Mesh when a Cell is removed.
	"""
        if cell in self.cells:
            if cell == self.cells[0]:
                self.cells = self.cells[:-1]
            else:
                self.cells = self.cells[1:]
        else:
            print "error in removeBoundingCell:"
            print "cell not in face.cells"
            
        
    def setId(self,id):
	"""Set the id of the Face.
	"""
        self.id = id
            
            
    def getOrientation(self):
        return self.orientation
