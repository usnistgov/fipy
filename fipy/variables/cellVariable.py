#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "cellVariable.py"
 #                                    created: 12/9/03 {2:03:28 PM} 
 #                                last update: 1/16/04 {11:02:07 AM} 
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

from variable import Variable
import Numeric

class CellVariable(Variable):
    def __init__(self, mesh, name = '', value=0., unit = None, hasOld = 1):
	array = Numeric.zeros([len(mesh.getCells())],'d')
# 	array[:] = value
	
	Variable.__init__(self, mesh = mesh, name = name, value = value, unit = unit, array = array)

	if hasOld:
	    self.old = self.copy()
	else:
            self.old = None
	    
	self.faceValue = self.grad = self.faceGrad = None
	
    def getVariableClass(self):
	return CellVariable

    def copy(self):
	    
	return CellVariable(
	    mesh = self.mesh, 
	    name = self.name + "_old", 
	    value = self.getValue(),
	    hasOld = 0)

    def getGridArray(self):
	return self.mesh.makeGridData(self.value)
	
    def setValue(self,value,cells = ()):
	if cells == ():
	    self[:] = value
# 	    if type(value) == type(Numeric.array((1.))):
# 		self[:] = value[:]
# 	    elif type(value) in [type(1.),type(1)]:
# 		self[:] = value
# 	    else:
# 		raise TypeError, str(value) + " is not numeric or a Numeric.array"
	else:
	    for cell in cells:
		self[cell.getId()] = value
	
    def getGrad(self):
	if self.grad is None:
	    from cellGradVariable import CellGradVariable
	    self.grad = CellGradVariable(self)
	
	return self.grad

    def getFaceValue(self):
	if self.faceValue is None:
	    from cellToFaceVariable import CellToFaceVariable
	    self.faceValue = CellToFaceVariable(self)

	return self.faceValue

    def getFaceGrad(self):
	if self.faceGrad is None:
	    from faceGradVariable import FaceGradVariable
	    self.faceGrad = FaceGradVariable(self)

	return self.faceGrad

    def getOld(self):
	if self.old is None:
	    return self
	else:
	    return self.old

    def updateOld(self):
	if self.old is not None:
	    self.old.setValue(self.value)
	    
