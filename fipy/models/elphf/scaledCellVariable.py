## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "scaledCellVariable.py"
 #                                    created: 6/2/04 {6:17:03 PM} 
 #                                last update: 7/28/04 {9:12:16 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  scaledCellVariable.py
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
 #  See the file "license.terms" for information on usage and  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2004-06-02 JEG 1.0 original
 # ###################################################################
 ##

 
from fipy.variables.cellVariable import CellVariable
from fipy.tools.dimensions import physicalField

class ScaledCellVariable(CellVariable):
    def __init__(self, mesh, name, value=0., hasOld = 1, scale = 1):
	self.scale = scale
	value = physicalField.Scale(value, self.scale)
	CellVariable.__init__(self, mesh = mesh, name = name, value = value, hasOld = hasOld)
	self.scaled = self * self.getScale()

    def getScale(self):
	return self.scale
	
    def setValue(self, value, cells = ()):
	return CellVariable.setValue(self, physicalField.Scale(value, self.getScale()), cells)
	
    def getScaled(self):
	return self.scaled

