#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "cell.py"
 #                                    created: 11/10/03 {3:23:11 PM} 
 #                                last update: 11/24/03 {10:25:57 AM} 
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

"""Cell within a mesh
"""
 
class Cell:
    """Cell within a mesh

	Cells are bounded by Faces.
    """
    
    def __init__(self, faces, id):
	"""Cell is initialized by Mesh with 
	
	Arguments:
	    
	    'faces' -- 'list' or 'tuple' of bounding faces that define the cell
	    
	    'id' -- unique identifier
	"""
        self.faces = faces
        self.id = id
        for face in self.faces:
            face.addBoundingCell(self)
        self.center = self.calcCenter()
## can not calculate face normals until this point
## the face needs to know its cells in order to
## calculate the normal orientation. Initially the
## normal goes from cell 1 to cell 2.
        for face in self.faces:
            face.setNormal()
            face.setCellDistance()
        self.volume = self.calcVolume()
        

    def getId(self):
	"""Return the id of this Cell.
	"""
        return self.id
	
    def calcVolume(self):
	"""Calculate the volume of the Cell.
	
	Sums the projected volumes of the faces.  
	
	Projected volume of a face is a right-prism bounded by its center
	and the y-z plane, whose cross-section is the projection of the face
	on the y-z plane.
	"""
	vol = 0.
	for face in self.faces:
	    vol += face.getCenter()[0] * face.getArea() * face.getNormal(self)[0]
	return vol

    def getVolume(self):
	"""Return the volume of the Cell.
	"""
        return self.volume

    def getCenter(self):
	"""Return the coordinates of the Cell center.
	"""
        return self.center

    def calcCenter(self):
	"""Calculate the coordinates of the Cell center.
	
	Cell center is the average of the bounding Face centers.
	
	**Is this right? Seems to give too much weight to small faces.**
	"""
        ctr = self.faces[0].getCenter().copy()
        for face in self.faces[1:]:
            ctr += face.getCenter()
        return ctr/float(len(self.faces))
            
    def __repr__(self):
	"""Textual representation of Cell.
	"""
	rep = "<id = " + str(self.id) + ", volume = " + str(self.volume()) + ", center = " + str(self.getCenter()) + ", faces = \n" 
	
	for face in self.faces:
	    rep += "id = " + str(face.getId()) + ", normal = " + str(face.getNormal(self)) + "\n"
	
	rep += ">\n"
	
	return rep

    def getFaces(self):
	"""Return the faces bounding the Cell.
	"""
        return self.faces

    def setId(self,id):
	"""Set the id of the Cell.
	"""
        self.id = id
    
