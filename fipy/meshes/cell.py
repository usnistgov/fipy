"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "cell.py"
                                   created: 11/10/03 {3:23:11 PM} 
                               last update: 11/18/03 {12:13:00 PM} 
 Author: Jonathan Guyer
 E-mail: guyer@nist.gov
 Author: Daniel Wheeler
 E-mail: daniel.wheeler@nist.gov
   mail: NIST
    www: http://ctcms.nist.gov
 
========================================================================
This software was developed at the National Institute of Standards
and Technology by employees of the Federal Government in the course
of their official duties.  Pursuant to title 17 Section 105 of the
United States Code this software is not subject to copyright
protection and is in the public domain.  PFM is an experimental
system.  NIST assumes no responsibility whatsoever for its use by
other parties, and makes no guarantees, expressed or implied, about
its quality, reliability, or any other characteristic.  We would
appreciate acknowledgement if the software is used.

This software can be redistributed and/or modified freely
provided that any derivative works bear some notice that they are
derived from it, and any modified versions bear some notice that
they have been modified.
========================================================================
 
 Description: 

 History

 modified   by  rev reason
 ---------- --- --- -----------
 2003-11-10 JEG 1.0 original
###################################################################
"""

class Cell:
    def __init__(self, faces, id):
        self.faces = faces
        self.id = id
        for face in self.faces:
            face.addBoundingCell(self)
        self.center = self.calcCenter()
## can not calculate face normals untill this point
## the face needs to know it's cells in order to
## calculate the normal orientation. Initially the
## normal goes from cell 1 to cell 2.
        for face in self.faces:
            face.setNormal()

    def getId(self):
        return self.id
	
    def volume(self):
	vol = 0.
	for face in self.faces:
	    vol += face.getCenter()[0] * face.area() * face.getNormal(self)[0]
	return vol

    def getCenter(self):
        return self.center

    def calcCenter(self):
        ctr = self.faces[0].getCenter().copy()
        for face in self.faces[1:]:
            ctr += face.getCenter()
        return ctr/float(len(self.faces))
            
    def __repr__(self):
	rep = "<id = " + str(self.id) + ", volume = " + str(self.volume()) + ", center = " + str(self.getCenter()) + ", faces = \n" 
	
	for face in self.faces:
	    rep += "id = " + str(face.getId()) + ", normal = " + str(face.getNormal(self)) + "\n"
	
	rep += ">\n"
	
	return rep

    def getFaces(self):
        return self.faces

    def setId(self,id):
        self.id = id
    
