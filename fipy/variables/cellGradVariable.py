"""
## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "cellGradVariable.py"
 #                                    created: 12/18/03 {2:28:00 PM} 
 #                                last update: 01/06/04 { 5:38:23 PM}
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
"""

from vectorCellVariable import VectorCellVariable
import Numeric

class CellGradVariable(VectorCellVariable):
    def __init__(self, var):
	VectorCellVariable.__init__(self, var.getMesh())
	self.var = self.requires(var)
	self.faceGradientContributions = self.mesh.getAreaProjections() * self.var.getFaceValue().transpose()
	
    def calcValue(self):
	N = len(self.var[:])
	M = self.mesh.getMaxFacesPerCell()
	
	ids = self.mesh.getCellFaceIDs()

	contributions = Numeric.take(self.faceGradientContributions[:], ids)
	contributions = Numeric.reshape(contributions,(N,M,self.mesh.getDim()))

	orientations = self.mesh.getCellFaceOrientations()

	grad = Numeric.sum(orientations*contributions,1)

	volumes = self.mesh.getCellVolumes()

	grad = grad/volumes[:,Numeric.NewAxis]

	self.value = grad
