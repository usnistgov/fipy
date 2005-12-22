#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "vectorArithmeticCellToFaceVariable.py"
 #                                    created: 7/26/04 {11:14:05 AM} 
 #                                last update: 12/22/05 {3:59:35 PM} 
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

from fipy.variables.vectorCellToFaceVariable import _VectorCellToFaceVariable
from fipy.tools import numerix
from fipy.tools.inline import inline

class _VectorArithmeticCellToFaceVariable(_VectorCellToFaceVariable):
    def _calcValuePy(self, alpha, id1, id2):
	cell1 = numerix.take(self.var,id1)
	cell2 = numerix.take(self.var,id2)
        return (cell1 - cell2) * alpha[:,numerix.NewAxis] + cell2
	
    def _calcValueIn(self, alpha, id1, id2):
##        print numerix.shape(self.value.value)
##        print numerix.shape(self.var.getNumericValue())
        val = self._getArray().copy()
        
	inline._runInlineLoop2("""
	    double cell2 = var(id2(i),j);
	    val(i,j) = (var(id1(i),j) - cell2) * alpha(i) + cell2;
	""",
	var = self.var.getNumericValue(),
	val = val,
	alpha = alpha,
	id1 = id1, id2 = id2,
	ni = self.mesh._getNumberOfFaces(),
        nj = self.mesh.getDim())

        return self._makeValue(value = val)
##         return self._makeValue(value = val, unit = self.getUnit())
