## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "vectorFaceDifferenceVariable.py"
 #                                    created: 3/17/05 {10:26:32 AM} 
 #                                last update: 12/22/05 {10:53:45 AM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
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
 #  2005-03-17 JEG 1.0 original
 # ###################################################################
 ##

from fipy.tools import numerix
from fipy.variables.vectorFaceVariable import VectorFaceVariable

class _VectorFaceDifferenceVariable(VectorFaceVariable):
    """
    Calculates the difference between the cell vector values adjacent to the face,
    weighted by the aspect ratio of the face.
    """
    
    def __init__(self, var):
	VectorFaceVariable.__init__(self, var.getMesh())
	self.var = self._requires(var)
	
    def _calcValue(self):
	id1, id2 = self.mesh._getAdjacentCellIDs()
	return self.mesh._getFaceAspectRatios()[:,numerix.NewAxis] * (numerix.take(self.var,id2) - numerix.take(self.var,id1))

