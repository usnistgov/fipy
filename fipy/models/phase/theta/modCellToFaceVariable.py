#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "modCellToFaceVariable.py"
 #                                    created: 12/18/03 {2:23:41 PM} 
 #                                last update: 9/3/04 {10:37:47 PM} 
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

from fipy.tools.inline import inline
from fipy.variables.arithmeticCellToFaceVariable import ArithmeticCellToFaceVariable

class ModCellToFaceVariable(ArithmeticCellToFaceVariable):
    def __init__(self, var, modIn):
	ArithmeticCellToFaceVariable.__init__(self,var)
        self.modIn = modIn
        
    def  _calcValueIn(self, alpha, id1, id2):
        
	inline._runInline(self.modIn + """
        int i;
        for(i = 0; i < ni; i++)
        {
	    double cell2 = var(id2(i));
	    val(i) = mod(var(id1(i)) - cell2) * alpha(i) + cell2;
        }
	""",var = self.var.getNumericValue(),
            val = self._getArray(), 
            alpha = alpha,
            id1 = id1, id2 = id2,
            ni = self.mesh._getNumberOfFaces())
