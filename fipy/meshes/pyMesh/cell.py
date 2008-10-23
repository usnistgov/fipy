#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cell.py"
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

"""Cell within a mesh
"""
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

class Cell:
    """Cell within a mesh

    `Cell` objects are bounded by `Face` objects.
    """
    
    def __init__(self, faces, faceOrientations, id):
	"""`Cell` is initialized by `Mesh` 
	
	:Parameters:
	    
	  - `faces`: `list` or `tuple` of bounding faces that define the cell
	
	  - `faceOrientations`: `list`, `tuple`, or `numerix.array` of
	    orientations (+/-1) to indicate whether a face points into this
	    face or out of it.  Can be calculated, but the mesh typically
	    knows this information already.
	
	  - `id`: unique identifier
	"""
        self.faces = faces
	self.faceOrientations = numerix.array(faceOrientations)
	self.faceOrientations = numerix.reshape(faceOrientations,(len(faces),1))
        self.id = id
	self.center = self._calcCenter()
	self.volume = self._calcVolume()
	for i in range(len(self.faces)):
	    self.faces[i].addBoundingCell(self,faceOrientations[i])

    def getID(self):
	"""Return the id of this Cell.
	"""
        return self.id
	
    def getFaceOrientations(self):
	return self.faceOrientations
	
    def _calcVolume(self):
	"""Calculate the volume of the Cell.
	
	Sums the projected volumes of the faces.  
	
	Projected volume of a face is a right-prism bounded by its center
	and the y-z plane, whose cross-section is the projection of the face
	on the y-z plane.
	"""
	vol = 0.
	for i in range(len(self.faces)):
	    face = self.faces[i]
	    vol += face.getCenter()[0] * face.getArea() * face.getNormal()[0] * self.faceOrientations[i][0]
	return vol

    def getVolume(self):
	"""Return the volume of the `Cell`.
	"""
        return self.volume

    def getCenter(self):
	"""Return the coordinates of the `Cell` center.
	"""
        return self.center

    def _calcCenter(self):
	"""Calculate the coordinates of the `Cell` center.
	
	Cell center is the average of the bounding Face centers.
	
	.. attention::
	    
	   Is this right? Seems to give too much weight to small faces.
	"""
        ctr = self.faces[0].getCenter().copy()
        for face in self.faces[1:]:
            ctr += face.getCenter()
        return ctr/float(len(self.faces))
            
    def __repr__(self):
	"""Textual representation of `Cell`.
	"""
	rep = "<id = " + str(self.id) + ", volume = " + str(self.getVolume()) + ", center = " + str(self.getCenter()) + ", faces = \n" 
	
	for i in range(len(self.faces)):
	    face = self.faces[i]
	    rep += "id = " + str(face.getID()) + ", normal = " + str(face.getNormal() * self.faceOrientations[i]) + "\n"
	
	rep += ">\n"
	
	return rep

    def getFaces(self):
	"""Return the faces bounding the `Cell`.
	"""
        return self.faces

    def _setID(self,id):
	"""Set the `id` of the `Cell`.
	"""
        self.id = id
    
    def getFaceIDs(self):
	return self.faceIDs
	
#    def _calcFaceIDs(self):
#	self.faceIDs = numerix.zeros(len(self.faces))
#	for i in range(len(self.faces)):
#	    self.faceIDs[i] = self.faces[i].getId()
            
    def _calcFaceIDs(self):
	self.faceIDs = ()
	for i in range(len(self.faces)):
	    self.faceIDs += (self.faces[i].getID(),)

    def getBoundingCells(self):
        boundingCells = ()
        for face in self.faces:
            faceCells = face.getCells()
            if len(faceCells) == 2:
                if faceCells[0].getID() == self.id:
                    boundingCells += (faceCells[1],)
                else:
                    boundingCells += (faceCells[0],)

        return boundingCells
                
                
