#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "face.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 11/21/03 {5:23:12 PM} 
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

    Faces are bounded by Vertices. Faces separate Cells.
"""

import tools
import Numeric

class Face:
    def __init__(self, vertices, id):
	"""Face is initialized by Mesh with its bounding vertices and a unique
	id.
	
	"""
        self.vertices = vertices
        self.cells = ()
        self.cellsId = ()
	self.id = id
	self.center = self.calcCenter()
        self.area = self.calcArea()
            
    def addBoundingCell(self, cell):
	"""Add cell to the list of Cells which lie on either side of this Face.
	
	"""
        self.cells += (cell,)
        self.cellsId += (cell.getId(), )
        
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
	    the origin and each polygon edge.  Properly calculates the area
	    of concave and convex polygons.
	"""
	a=0.
	p1 = self.vertices[0].getCoordinates().copy()
	for vertex in self.vertices[2:-1]:
	    p2=vertex.getCoordinates().copy()
	    a += tools.crossProd(p1,p2)
	    p1 = p2
	return abs(a/2.)
        
    def orientNormal(self, norm, cell):
	"""Determine if normal points into or out of the cell in question.
	
	    *Maybe the cell should keep track of this, rather than the face?*
	"""
        if len(self.cells) == 0:
            return "abnormal"
	elif len(self.cells) == 1:
	    """Boundary faces only have one cell
	    
	    center-to-center vector is from face center to cell center
	    """
	    cc = self.cells[0].getCenter() - self.center
	else:
	    cc = self.cells[0].getCenter() - self.cells[1].getCenter()
            
	if cell == 'None' or cell == self.cells[0]:
	    cc *= -1
            
	if Numeric.dot(cc,norm) < 0:
	    norm *= -1
	    
	return norm

    def getNormal(self, cell = 'None'):
	"""Return the unit normal vector, accounting for whether the Face
	points toward or away from the specified Cell.	
	"""
        if cell == self.cells[0] or cell == 'None':
            return self.normal
        else:
            return -self.normal

    def setNormal(self):
	"""Cache the unit normal vector.
	
	Called by Cell initializer after bounding Cells have been added to
	Faces.
	"""
        self.normal = self.calcNormal(self.cells[0])
	
    def calcNormal(self, cell = 'None'):
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
	
	return self.orientNormal(norm, cell)

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

    def __repr__(self):
	"""Textual representation of Face.
	"""
	return "<id = " + str(self.id) + ", area = " + str(self.area()) + ", normal = " + str(self.normal()) + ", vertices = " + str(self.vertices) + ", centers = " + str(self.center) + ">\n"

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
            
            
