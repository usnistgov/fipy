#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "harmonicCellToFaceVariable.py"
 #                                    created: 2/20/04 {11:15:10 AM} 
 #                                last update: 7/24/04 {9:01:59 AM} 
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

from fipy.variables.cellToFaceVariable import CellToFaceVariable
from fipy.tools import array
from fipy.tools.inline import inline

class HarmonicCellToFaceVariable(CellToFaceVariable):
    def _calcValuePy(self, alpha, id1, id2):
	cell1 = array.take(self.var,id1)
	cell2 = array.take(self.var,id2)
	self.value = ((cell2 - cell1) * alpha + cell1)
	if self.value != 0:
	    self.value = cell1 * cell2 / self.value
	
    def _calcValueIn(self, alpha, id1, id2):
	inline.runInlineLoop1("""
	    double	cell1 = var(id1(i));
	    double	cell2 = var(id2(i));
	    double	tmp = ((cell2 - cell1) * alpha(i) + cell1);
	    if (tmp != 0) {
		val(i) = cell1 * cell2 / tmp;
	    } else {
		val(i) = tmp;
	    }
	""",
	var = self.var.getNumericValue(),
	val = self._getArray(), 
	alpha = alpha,
	id1 = id1, id2 = id2,
	ni = len(self.mesh.getFaces())
	)
