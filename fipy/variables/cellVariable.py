#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "cellVariable.py"
 #                                    created: 12/9/03 {2:03:28 PM} 
 #                                last update: 3/9/04 {12:01:03 PM} 
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

import Numeric

from fivol.variables.variable import Variable
import fivol.tools.array

class CellVariable(Variable):
    def __init__(self, mesh, name = '', value=0., unit = None, hasOld = 0):
	array = Numeric.zeros([mesh.getNumberOfCells()],'d')
# 	array[:] = value
	
	Variable.__init__(self, mesh = mesh, name = name, value = value, unit = unit, array = array)

	if hasOld:
	    self.old = self.copy()
	else:
            self.old = None
	    
	self.arithmeticFaceValue = self.harmonicFaceValue = self.grad = self.faceGrad = None
	
    def getVariableClass(self):
	return CellVariable

    def copy(self):

        return self.__class__(
            mesh = self.mesh, 
	    name = self.name + "_old", 
	    value = self.getValue(),
	    hasOld = 0)
	    
##	return CellVariable(
##	    mesh = self.mesh, 
##	    name = self.name + "_old", 
##	    value = self.getValue(),
##	    hasOld = 0)

    def getGridArray(self):
	return self.mesh.makeGridData(self.value)
	
    def __call__(self, point, order = 0):
	d = self.mesh.getPointToCellDistances(point)
	i = Numeric.argsort(d)
	return self[i[0]]
## 	return (self[i[0]] * self[i[1]] * (d[i[0]] + d[i[1]])) / (self[i[0]] * d[i[0]] + self[i[1]] * d[i[1]])
	
    def getValue(self, points = (), cells = ()):
	if points == () and cells == ():
	    return Variable.getValue(self)
	elif cells != ():
	    return fivol.tools.array.take(Variable.getValue(self), [cell.getID() for cell in cells])
	else:
	    return [self(point) for point in points]
	
    def setValue(self,value,cells = ()):
	if cells == ():
	    self[:] = value
	else:
## 	    return fivol.tools.array.put(self.getValue(), [cell.getID() for cell in cells], value)
## 	    self.markStale()

	    for cell in cells:
		self[cell.getID()] = value
	
    def getGrad(self):
	if self.grad is None:
	    from cellGradVariable import CellGradVariable
	    self.grad = CellGradVariable(self)
        
	return self.grad

    def getArithmeticFaceValue(self):
	if self.arithmeticFaceValue is None:
	    from arithmeticCellToFaceVariable import ArithmeticCellToFaceVariable
	    self.arithmeticFaceValue = ArithmeticCellToFaceVariable(self)

	return self.arithmeticFaceValue

    def getHarmonicFaceValue(self):
	if self.harmonicFaceValue is None:
	    from harmonicCellToFaceVariable import HarmonicCellToFaceVariable
	    self.harmonicFaceValue = HarmonicCellToFaceVariable(self)

	return self.harmonicFaceValue

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
	    
    def resetToOld(self):
	if self.old is not None:
	    self.setValue(self.old.value)
