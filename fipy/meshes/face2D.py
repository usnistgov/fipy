#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "face2D.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 12/1/03 {3:54:52 PM} 
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

"""1D (edge) Face in a 2D Mesh
"""

import tools
from face import Face
import Numeric

class Face2D(Face):
    """1D (edge) Face in a 2D Mesh

	Face2D is bounded by two Vertices.
    """
    
    def calcArea(self):
	"""Area is length of vector between vertices.
	"""
        tangent=self.vertices[0].getCoordinates()-self.vertices[1].getCoordinates()
        return tools.sqrtDot(tangent,tangent)
	
    def calcNormal(self):
	"""Normal is perpendicular to vector between vertices.
	"""
	tangent = self.vertices[1].getCoordinates() - self.vertices[0].getCoordinates()
	norm = Numeric.array([-tangent[1],tangent[0]])
	norm /= tools.sqrtDot(norm,norm)
## we calculate the orientation after we know the normal
##	norm *= self.orientation

	return norm

    def calcTangent1(self):
	norm = self.normal
	mag = Numeric.sqrt(norm[0]**2 + norm[1]**2)
	tan1 = Numeric.array((-norm[1],norm[0]))
	return tan1/mag
	    
    def calcTangent2(self):
	return Numeric.array((0.,0.))
