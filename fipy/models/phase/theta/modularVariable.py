#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "modularVariable.py"
 #                                    created: 12/8/03 {5:47:27 PM} 
 #                                last update: 9/3/04 {10:35:53 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
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

 
import Numeric

from fipy.variables.cellVariable import CellVariable
from fipy.variables.cellGradVariable import CellGradVariable
from fipy.models.phase.theta.modCellGradVariable import ModCellGradVariable
from modCellToFaceVariable import ModCellToFaceVariable
from fipy.models.phase.theta.modPhysicalField import ModPhysicalField

class ModularVariable(CellVariable):
    def __init__(self, mesh, name = '', value=0., unit = None, hasOld = 0):
	CellVariable.__init__(self, mesh = mesh, name = name, value = value, unit = unit, hasOld = hasOld)
        self.arithmeticFaceValue = None
        self.grad = None
        
    _modIn = """
    # define pi 3.141592653589793
    # define mod(x) (fmod(x + 3. * pi, 2. * pi) - pi)
    """
	
    def _setValue(self, value, unit = None, array = None):
	self.value = ModPhysicalField(value = value, unit = unit, array = array)
	
    def updateOld(self):
	self.setValue(self.value.mod(self()))
        if self.old != None:
	    self.old.setValue(self())

    def getGrad(self):
	if self.grad is None:
##	    gridSpacing = self.mesh.getMeshSpacing()
##            self.grad = self.value.mod(CellGradVariable(self) * gridSpacing) / gridSpacing
            self.grad = ModCellGradVariable(self, self._modIn, self.value.mod)

	return self.grad

    def getArithmeticFaceValue(self):
	if self.arithmeticFaceValue is None:
	    from modCellToFaceVariable import ModCellToFaceVariable
	    self.arithmeticFaceValue = ModCellToFaceVariable(self, self._modIn)

 	return self.arithmeticFaceValue

    def getFaceGrad(self):
	if self.faceGrad is None:
	    from modFaceGradVariable import ModFaceGradVariable
	    self.faceGrad = ModFaceGradVariable(self, self._modIn)

	return self.faceGrad
