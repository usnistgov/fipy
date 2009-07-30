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

__docformat__ = 'restructuredtext'

class Face:
    """`Face` within a `Mesh`

    `Face` objects are bounded by `Vertex` objects. 
    `Face` objects separate `Cell` objects.
    """
    
    def __init__(self, mesh, id):
	"""
	`Face` is initialized by `Mesh`
	
	:Parameters:
	    
	  - `mesh`: the `Mesh` that contains this `Face`
	  - `id`:   a unique identifier
	"""
        import warnings
        warnings.warn("Face masks (e.g., `(x < 0.5) & (y > 1.3)`) or ID numbers should be used instead of 'Face'", DeprecationWarning, stacklevel=2)
        self.mesh = mesh
	self.id = id
    
    def getMesh(self):
	return self.mesh
		
    def getID(self):
	return self.id

    def getCellID(self, index = 0):
	"""Return the `id` of the specified `Cell` on one side of this `Face`.
	"""
	return self.mesh.getFaceCellIDs()[index, self.id]

    def getCenter(self):
 	"""Return the coordinates of the `Face` center.
 	"""
 	return self.mesh.getFaceCenters()[...,self.id]

     
    def getArea(self):
	return self.mesh._getFaceAreas()[...,self.id]
##     
##     def getNormal(self):
## 	"""Return the unit normal vector
## 	"""
## 	return self.normal
## 		
##     def getCellDistance(self):
## 	"""Return the distance between adjacent cell centers.
## 	"""
##         return self.cellDistance
## 
##     def _getFaceToCellDistance(self, cellId = None):
##         if cellId == self.cellsId[0] or cellId == None:
##             return self.faceToCellDistances[0]
##         elif cellId == self.cellsId[1]:
##             return self.faceToCellDistances[1]
##         else:
##             raise Exception
## 
##     def __repr__(self):
## 	"""Textual representation of Face.
## 	"""
## 	return "<id = " + str(self.id) + ", area = " + str(self.area) + ", normal = " + str(self.normal) + ", vertices = " + str(self.vertices) + ", center = " + str(self.center) + ">\n"
##             
##     def _getOrientation(self):
##         return self.orientation
