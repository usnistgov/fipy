#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "term.py"
 #                                    created: 11/12/03 {10:54:37 AM} 
 #                                last update: 5/13/04 {10:20:23 AM} 
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
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

import Numeric

from fipy.variables.variable import Variable
from fipy.tools.dimensions.physicalField import PhysicalField

class Term:
    def __init__(self,mesh, weight):
	self.mesh = mesh
	self.weight = weight
	self.calcCoeffScale()
        
    def buildMatrix(self, L, oldArray, b):
	pass
	
    def getMesh(self):
	return self.mesh
	
    def getCoeffScale(self):
##        self.calcCoeffScale()
        return self.coeffScale

    def calcCoeffScale(self):
        if isinstance(self.coeff, PhysicalField) or isinstance(self.coeff, Variable):
	    self.coeffScale = PhysicalField(1, self.coeff.getUnit())
	else:
	    self.coeffScale = 1
    
