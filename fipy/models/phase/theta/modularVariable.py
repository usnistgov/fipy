"""
## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "modularVariable.py"
 #                                    created: 12/8/03 {5:47:27 PM} 
 #                                last update: 12/19/03 {3:06:22 PM} 
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

from variables.cellVariable import CellVariable
import Numeric

class ModularVariable(CellVariable):
    def __init__(self, mesh, name = '', value=0., scaling = None, unit = None, hasOld = 1):
	CellVariable.__init__(self, mesh = mesh, name = name, value = value, scaling = scaling, unit = unit, hasOld = hasOld)

        self.mod = lambda array, pi=Numeric.pi: (array + 3. * pi) % (2 * pi) - pi

    def getGrad(self):
	if self.grad is None:
	    from variables.cellGradVariable import CellGradVariable
	    gridSpacing = self.mesh.getMeshSpacing()
	    self.grad = self.mod(CellGradVariable(self) * gridSpacing) / gridSpacing 
	
	return self.grad

    def getFaceValue(self):
	if self.faceValue is None:
	    from variables.cellToFaceVariable import CellToFaceVariable
	    self.faceValue = CellToFaceVariable(self, self.mod)

	return self.faceValue

    def getFaceGrad(self):
	if self.faceGrad is None:
	    from variables.faceGradVariable import FaceGradVariable
	    self.faceGrad = FaceGradVariable(self, self.mod)

	return self.faceGrad
	
    def updateOld(self):
        if self.old != None:
            self.value = self.mod(self.value)
	    self.old.setValue(self.value)

    
	


