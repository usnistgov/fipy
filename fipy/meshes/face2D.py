"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "face2D.py"
                                   created: 11/10/03 {3:23:47 PM}
                               last update: 11/20/03 {10:37:12 AM} 
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
from face import Face
import Numeric

class Face2D(Face):
    def __init__(self, vertices, id):
        Face.__init__(self,vertices,id)
	
    def calcArea(self):
        tangent=self.vertices[0].getCoordinates()-self.vertices[1].getCoordinates()
        return tools.sqrtDot(tangent,tangent)
	
    def calcNormal(self, cell = 'None'):
	tangent = self.vertices[1].getCoordinates() - self.vertices[0].getCoordinates()
	norm = Numeric.array([-tangent[1],tangent[0]])
	norm /= tools.sqrtDot(norm,norm)
	    
	return self.orientNormal(norm, cell)
