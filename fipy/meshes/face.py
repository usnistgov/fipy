"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "face.py"
                                   created: 11/10/03 {3:23:47 PM}
                               last update: 11/18/03 {12:05:49 PM} 
 Author: Jonathan Guyer
 E-mail: guyer@nist.gov
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

import tools
import Numeric

class Face:
    def __init__(self, vertices, id):
        self.vertices = vertices
        self.cells = ()
	self.id = id
	
    def addBoundingCell(self, cell):
        self.cells += (cell,)

    def getCells(self):
        return self.cells
		
    def getId(self):
	return self.id

    def getCellId(self,index=0):
        return self.cells[index].getId()
	
    def center(self):
	ctr = self.vertices[0].getCoordinates().copy()
	for vertex in self.vertices[1:]:
	    ctr += vertex.getCoordinates()
	return ctr / float(len(self.vertices))
	
    def area(self):
	a=0.
	p1 = self.vertices[0].getCoordinates()
	for vertex in self.vertices[2:-1]:
	    p2=vertex.getCoordinates()
	    a += tools.crossProd(p1,p2)
	    p1 = p2
	return abs(a/2.)
        
    def orientNormal(self, norm, cell):
	"""
	Determine if normal points into or out of the cell in question
	"""
	if len(self.cells) == 1:
	    """Boundary faces only have one cell
	    
	    center-to-center vector is from face center to cell center
	    """
	    cc = self.cells[0].center() - self.center()
	else:
            print len(self.cells)
	    cc = self.cells[0].center() - self.cells[1].center()
	if cell == None or cell == self.cells[0]:
	    cc *= -1
	if Numeric.dot(cc,norm) < 0:
	    norm *= -1
	    
	return norm
	
	
    def normal(self, cell = None):	
	t1 = self.vertices[1].getCoordinates() - self.vertices[0].getCoordinates()
	t2 = self.vertices[2].getCoordinates() - self.vertices[1].getCoordinates()
	norm = tools.crossProd(t1,t2)
	norm /= tools.sqrtDot(norm,norm)
	
	return self.orientNormal(norm, cell)

    def cellDistance(self):
        if(len(self.cells)==2):
            vec=self.cells[1].center()-self.cells[0].center()
        else:
            vec=self.center()-self.cells[0].center()        
        return tools.sqrtDot(vec,vec)

    def __repr__(self):
	return "<id = " + str(self.id) + ", area = " + str(self.area()) + ", normal = " + str(self.normal(self.cells[0])) + ", vertices = " + str(self.vertices) + ">\n"


    def removeBoundingCell(self,cell):
        if cell in self.cells:
            if cell == self.cells[0]:
                self.cells = self.cells[:-1]
            else:
                self.cells = self.cells[1:]
        else:
            print "error in removeBoundingCell:"
            print "cell not in face.cells"
            
        
    def setId(self,id):
        self.id = id
            
            
