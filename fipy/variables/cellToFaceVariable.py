## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "cellToFaceVariable.py"
 #                                    created: 12/18/03 {2:23:41 PM} 
 #                                last update: 12/18/03 {3:23:02 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

from faceVariable import FaceVariable
import Numeric

class CellToFaceVariable(FaceVariable):
    def __init__(self, var):
	FaceVariable.__init__(self, var.getMesh())
	self.var = self.requires(var)

    def calcValue(self):
	alpha = self.mesh.getFaceToCellDistanceRatio()
	id1, id2 = self.mesh.getAdjacentCellIDs()
	cell1 = Numeric.take(self.var[:], id1)
	cell2 = Numeric.take(self.var[:], id2)
	self.value = (cell1 - cell2) * alpha + cell2
	self.value = Numeric.reshape(self.value, (len(self.value),1))
