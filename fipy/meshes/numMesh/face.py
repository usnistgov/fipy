#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "face.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 3/5/04 {1:45:57 PM} 
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

class Face:
    """Face within a Mesh

    Faces are bounded by Vertices. Faces separate Cells.
    """
    
    def __init__(self, mesh, id):
	"""Face is initialized by Mesh
	
	Arguments:
	    
	    'mesh' -- the 'Mesh' that contains this 'Face'
	    
	    'id' -- a unique identifier
	"""
        self.mesh = mesh
	self.id = id
	
##     def getCells(self):
## 	"""Return the Cells which lie on either side of this Face.
## 	"""
##         return self.cells
		
    def getID(self):
	""" If you can't figure out what the method does, then Daniel thinks you're an idiot.
	"""
	return self.id

    def getCellID(self, index = 0):
	"""Return the id of the specified Cell on one side of this Face.
	"""
	return self.mesh.getFaceCellIDs()[self.id,index]

    def getCenter(self):
 	"""Return the coordinates of the Face center.
 	"""
 	return self.mesh.getFaceCenters()[self.id]

     
    def getArea(self):
	return self.mesh.getFaceAreas()[self.id]
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
##     def getFaceToCellDistance(self, cellId = None):
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
##     def getOrientation(self):
##         return self.orientation
