"""
## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "CellVariable.py"
 #                                    created: 12/9/03 {2:03:28 PM} 
 #                                last update: 12/19/03 {4:04:46 PM} 
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
from variable import Variable
import Numeric

class CellVariable(Variable):
    def __init__(self, mesh, name = '', value=0., scaling = None, unit = None, viewer = None, hasOld = 1):
	array = Numeric.zeros([len(mesh.getCells())],'d')
# 	array[:] = value
	
	if viewer is not None:
	    self.viewer = viewer(var = self)
	else:
	    self.viewer = None
	    
	Variable.__init__(self, mesh, name = name, value = value, array = array, scaling = scaling, unit = unit)

	if hasOld:
	    self.old = self.copy()
	else:
            self.old = None
	    
	self.faceValue = self.grad = self.faceGrad = None
	
    def getVariableClass(self):
	return CellVariable

    def copy(self):
	if self.viewer is None:
	    viewer = None
	else:
	    viewer = self.viewer.__class__
	    
	return CellVariable(
	    mesh = self.mesh, 
	    name = self.name + "_old", 
	    value = self.getValue(),
	    scaling = self.scaling,
	    viewer = viewer,
	    hasOld = 0)

    def plot(self):
	if self.viewer != None:
	    self.viewer.plot()
	
    def getGridArray(self):
	return self.mesh.makeGridData(self.value)
	
    def setValue(self,value,cells = ()):
	if cells == ():
	    if type(value) == type(Numeric.array((1.))):
		self[:] = value[:]
	    elif type(value) in [type(1.),type(1)]:
		self[:] = value
	    else:
		raise TypeError, str(value) + " is not numeric or a Numeric.array"
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
	if self.old == None:
	    return self
	else:
	    return self.old

    def updateOld(self):
	if self.old != None:
	    self.old.setValue(self.value)
	    
