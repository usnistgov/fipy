#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "vectorCellToFaceVariable.py"
 #                                    created: 7/26/04 {2:23:41 PM} 
 #                                last update: 7/26/04 {4:00:53 PM} 
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

from fipy.variables.vectorFaceVariable import VectorFaceVariable
from fipy.tools import array
from fipy.tools.inline import inline

class VectorCellToFaceVariable(VectorFaceVariable):
    def __init__(self, var):
	VectorFaceVariable.__init__(self, var.getMesh())
	self.var = self.requires(var)

    def calcValue(self):
	alpha = self.mesh.getFaceToCellDistanceRatio()
	id1, id2 = self.mesh.getAdjacentCellIDs()
	inline.optionalInline(self._calcValueIn, self._calcValuePy, alpha, id1, id2)

	
