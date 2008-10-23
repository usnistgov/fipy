#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "face.py"
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
 #  
 # ###################################################################
 ##

"""`Face` within a `Mesh`
"""
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.tools.dimensions.physicalField import PhysicalField

class Face:
    """`Face` within a `Mesh`

    `Face` objects are bounded by `Vertex` objects. 
    `Face` objects separate `Cell` objects.
    """
    
    def __init__(self, vertices, id):
	"""`Face` is initialized by `Mesh`
	
	:Parameters:
	    
	  - `vertices`: the `Vertex` points that bound the `Face`
	  - `id`: a unique identifier
	"""
        self.vertices = vertices
##        self.cells = ()
        self.cellsID = ()
        self.cellCenters = ()
	self.id = id
	self.center = self._calcCenter()
        self.area = self._calcArea()
	self._setNormal()
	self.orientation = 1
            
##    def addBoundingCell(self, cell, orientation):
##	"""Add cell to the list of Cells which lie on either side of this Face.
##	"""
##        self.cells += (cell,)
##        self.cellsId += (cell.getId(),)
##	if len(self.cells) == 1:
##	    self.orientation = orientation
##	else:
##	    self.orientation = -orientation	    
##	self._setCellDistance()
##	self._setFaceToCellDistances()

    def addBoundingCell(self, cell, orientation):
	"""Add `cell` to the list of `Cell` objects which lie 
	on either side of this Face.
	"""
##        self.cells += (cell,)
        self.cellsID += (cell.getID(),)
	if len(self.cellsID) == 1:
	    self.orientation = orientation
	else:
	    self.orientation = -orientation
        self.cellCenters += (cell.getCenter(),)
	self._setCellDistance()
	self._setFaceToCellDistances()
        
    def getCells(self):
	"""Return the `Cell` objects which lie on either side of this `Face`.
	"""
        return self.cells
		
    def getID(self):
	return self.id

    def getCellID(self, index = 0):
	"""Return the `id` of the specified `Cell` on one side of this `Face`.
	"""
        if index == 1 and len(self.cellsID) == 1:
            index = 0
        return self.cellsID[index]

##    def getCellId(self,index=0):
##        return self.cells[index].getId()
	
    def getCenter(self):
	"""Return the coordinates of the `Face` center.
	"""
	return self.center
    
    def _calcCenter(self):
	"""Calculate the coordinates of the `Face` center.
	
	Cell center is the average of the bounding Vertex centers.
	"""
	ctr = self.vertices[0].getCoordinates().copy()
	for vertex in self.vertices[1:]:
	    ctr += vertex.getCoordinates()
	return ctr / float(len(self.vertices))

    def getArea(self):
	"""Return the area of the `Face`.
	"""
        return self.area
    
    def _calcArea(self):
	"""Calculate the area of the `Face`.
	
	Area is the signed sum of the area of the triangles bounded by
	each polygon edge and the origin.
	"""
	a=0.
	p1 = self.vertices[0].getCoordinates().copy()
	for vertex in self.vertices[2:-1]:
	    p2=vertex.getCoordinates().copy()
	    a += numerix.cross(p1,p2)
	    p1 = p2
	return abs(a/2.)
        
    def getNormal(self):
	"""Return the unit normal vector
	"""
	return self.normal
		
    def _setNormal(self):
	"""Cache the unit normal vector.
	
	Called by Cell initializer after bounding Cells have been added to
	Faces.
	"""
        self.normal = self._calcNormal()
	
    def _calcNormal(self):
	"""Calculate the unit normal vector.	
	
	Unit normal vector is calculated from cross-product of two
	tangent vectors.
    
	.. warning::
	   
	   Doesn't work if t1 and t2 are colinear!
	"""
	t1 = self.vertices[1].getCoordinates() - self.vertices[0].getCoordinates()
	t2 = self.vertices[2].getCoordinates() - self.vertices[1].getCoordinates()
	norm = numerix.crossProd(t1,t2)
	norm /= numerix.sqrtDot(norm,norm)
## we calculate the orientation after we know the normal
##	norm *= self.orientation
	
	return norm
	
    def _calcTangent1(self):
	norm = self.normal
	mag = numerix.sqrt(norm[0]**2 + norm[1]**2)
# 	tan1 = numerix.array((-norm[1],norm[0],0))
	tan1 = PhysicalField(value = (-norm[1],norm[0],0))
	return tan1/mag

    def _calcTangent2(self):
	norm = self.normal
	mag = numerix.sqrt(norm[0]**2 + norm[1]**2)
# 	tan2 = numerix.array(norm[0] * norm[2], norm[1] * norm[2], -mag**2)
	tan2 = PhysicalField(value = (norm[0] * norm[2], norm[1] * norm[2], -mag**2))
	return tan2/mag
	    
    def getCellDistance(self):
	"""Return the distance between adjacent `Cell` centers.
	"""
        return self.cellDistance

##    def _setCellDistance(self):
##	"""Assign the cached distance between adjacent cell centers.
##	"""
##        self.cellDistance = self._calcCellDistance()

    
    def _setCellDistance(self):
	"""Assign the cached distance between adjacent `Cell` centers.
	"""
        self.cellDistance = self._calcCellDistance()


##    def _calcCellDistance(self):
##	"""Calculate the distance between adjacent Cell centers.
	
##	If the Face is on a boundary and has only one bordering Cell,
##	the distance is from the Cell center to the Face center.
##	"""
##        if(len(self.cells)==2):
##            vec=self.cells[1].getCenter()-self.cells[0].getCenter()
##        else:
##            vec=self.center-self.cells[0].getCenter()        
##        return numerix.sqrtDot(vec,vec)

    def _calcCellDistance(self):
	"""Calculate the distance between adjacent `Cell` centers.
	
	If the Face is on a boundary and has only one bordering Cell,
	the distance is from the Cell center to the Face center.
	"""
        if(len(self.cellsID)==2):
            vec=self.cellCenters[0]-self.cellCenters[1]
        else:
            vec=self.center-self.cellCenters[0]        
        return numerix.sqrtDot(vec,vec)

##    def _calcFaceToCellDistance(self, cell):
##        vec=self.center-cell.getCenter()        
##        return numerix.sqrtDot(vec,vec)

    def _calcFaceToCellDistance(self, id):
        vec=self.center - self.cellCenters[id]
        return numerix.sqrtDot(vec,vec)

##    def _setFaceToCellDistances(self):
##        faceToCellDistances = ()
##        for cell in self.cells:
##            faceToCellDistances += (self._calcFaceToCellDistance(cell),)
##        self.faceToCellDistances = faceToCellDistances

    def _setFaceToCellDistances(self):
        faceToCellDistances = ()
        for id in range(len(self.cellsID)):
            faceToCellDistances += (self._calcFaceToCellDistance(id),)
        self.faceToCellDistances = faceToCellDistances
        

##    def _getFaceToCellDistance(self, cell = None):
##        if cell == self.cells[0] or cell == None:
##            return self.faceToCellDistances[0]
##        elif cell == self.cells[1]:
##            return self.faceToCellDistances[1]
##        else:
##            return self._calcFaceToCellDistance(cell)

    def _getFaceToCellDistance(self, cellID = None):
        if cellID == self.cellsID[0] or cellID == None:
            return self.faceToCellDistances[0]
        elif cellID == self.cellsID[1]:
            return self.faceToCellDistances[1]
        else:
            raise Exception


#     def __repr__(self):
# 	"""Textual representation of Face.
# 	"""
# 	return "<id = " + str(self.id) + ">"
# 	+ ", area = " + str(self.area()) + ", normal = " + str(self.normal()) + ", vertices = " + str(self.vertices) + ", centers = " + str(self.center) + ">\n"
    def __repr__(self):
	"""Textual representation of `Face`.
	"""
	return "<id = " + str(self.id) + ", area = " + str(self.area) + ", normal = " + str(self.normal) + ", vertices = " + str(self.vertices) + ", center = " + str(self.center) + ">\n"

    def _removeBoundingCell(self,cell):
	"""Remove `cell` from the list of bounding `Cell` objects.
	
	Called by the Mesh when a Cell is removed.
	"""
        if cell in self.cells:
            if cell == self.cells[0]:
                self.cells = self.cells[:-1]
            else:
                self.cells = self.cells[1:]
        else:
            print "error in _removeBoundingCell:"
            print "cell not in face.cells"
            
        
    def _setID(self,id):
	"""Set the `id` of the `Face`.
	"""
        self.id = id
            
    def _getOrientation(self):
        return self.orientation
