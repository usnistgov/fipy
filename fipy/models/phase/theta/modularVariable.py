#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "modularVariable.py"
 #                                    created: 12/8/03 {5:47:27 PM} 
 #                                last update: 1/16/04 {9:20:28 PM} 
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

from Numeric import pi

from fivol.variables.cellVariable import CellVariable
from fivol.variables.cellGradVariable import CellGradVariable
from modCellToFaceVariable import ModCellToFaceVariable
from fivol.examples.phase.theta.modPhysicalField import ModPhysicalField


class ModularVariable(CellVariable):
    def __init__(self, mesh, name = '', value=0., unit = None, hasOld = 1):
	CellVariable.__init__(self, mesh = mesh, name = name, value = value, unit = unit, hasOld = hasOld)
        self.faceValue = None
        self.grad = None
        self.mod = """
        # define pi 3.141592653589793
        # define mod(x) (fmod(x + 3. * pi, 2. * pi) - pi)
        """
        
    def getPhysicalFieldClass(self):
	return ModPhysicalField

    def updateOld(self):
        self.value.value = self.value.mod(self.value)
        if self.old != None:
	    self.old.setValue(self.value.value)

    def getGrad(self):
	if self.grad is None:
	    gridSpacing = self.mesh.getMeshSpacing()
            self.grad = self.value.mod(CellGradVariable(self) * gridSpacing) / gridSpacing 

	return self.grad

    def getFaceValue(self):
	if self.faceValue is None:
	    from modCellToFaceVariable import ModCellToFaceVariable
	    self.faceValue = ModCellToFaceVariable(self, self.mod)

 	return self.faceValue

    def getFaceGrad(self):
	if self.faceGrad is None:
	    from modFaceGradVariable import ModFaceGradVariable
	    self.faceGrad = ModFaceGradVariable(self, self.mod)

	return self.faceGrad
