## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "vectorFaceDifferenceVariable.py"
 #                                    created: 3/17/05 {10:26:32 AM} 
 #                                last update: 4/2/05 {7:30:55 PM} 
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

from fipy.tools import array
from fipy.variables.vectorFaceVariable import VectorFaceVariable
import Numeric

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
        print self.mesh._getFaceAspectRatios()[:,Numeric.NewAxis].shape, (array.take(self.var,id2) - array.take(self.var,id1)).shape
	self.value = self.mesh._getFaceAspectRatios()[:,Numeric.NewAxis] * (array.take(self.var,id2) - array.take(self.var,id1))

