## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "finiteDifferenceAspectRatioGradient.py"
 #                                    created: 4/30/04 {2:01:50 PM} 
 #                                last update: 5/7/04 {9:54:51 AM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  finiteDifferenceAspectRatioGradient.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2004-04-30 JEG 1.0 original
 # ###################################################################
 ##

import fipy.tools.array as array
from fipy.variables.faceVariable import FaceVariable

class FaceDifferenceVariable(FaceVariable):
    """
    Calculates the difference between the cell values adjacent to the face,
    weighted by the aspect ratio of the face.
    """
    
    def __init__(self, var):
	FaceVariable.__init__(self, var.getMesh())
	self.var = self.requires(var)
	
    def calcValue(self):
	id1, id2 = self.mesh.getAdjacentCellIDs()
	self.value = self.mesh.getFaceAspectRatios()[:] * (array.take(self.var,id2) - array.take(self.var,id1))

