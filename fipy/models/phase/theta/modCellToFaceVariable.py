#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "modCellToFaceVariable.py"
 #                                    created: 12/18/03 {2:23:41 PM} 
 #                                last update: 2/17/04 {5:58:17 PM} 
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

from fivol.inline import inline
from fivol.variables.cellToFaceVariable import CellToFaceVariable

class ModCellToFaceVariable(CellToFaceVariable):
    def __init__(self, var = None, mod = None):
        CellToFaceVariable.__init__(self,var)
        self.mod = mod

    def  _calcValueIn(self, alpha, id1, id2):
	inline.runInlineLoop1(self.mod + """
        double	cell2 = var(id2(i));
        val(i) = mod(var(id1(i)) - cell2) * alpha(i) + cell2;
	""",var = self.var.getNumericValue(),
            val = self.value.value, 
            alpha = alpha,
            id1 = id1, id2 = id2,
            ni = len(self.mesh.getFaces()))
