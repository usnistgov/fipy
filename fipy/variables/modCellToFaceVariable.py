#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "modCellToFaceVariable.py"
 #
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

__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix
from fipy.tools import inline
from fipy.variables.arithmeticCellToFaceVariable import _ArithmeticCellToFaceVariable

class _ModCellToFaceVariable(_ArithmeticCellToFaceVariable):
    def __init__(self, var, modIn):
        _ArithmeticCellToFaceVariable.__init__(self,var)
        self.modIn = modIn
        
    if inline.doInline:
        def  _calcValue_(self, alpha, id1, id2):
            val = self._array.copy()
            
            inline._runInline(self.modIn + """
            int ID1 = id1[i];
            int ID2 = id2[i];
            double cell2 = var[ID2];
            val[i] = mod(cell2 - var[ID1]) * alpha[i] + var[ID1];
            """,var = self.var.numericValue,
                val = val, 
                alpha = alpha,
                id1 = id1, id2 = id2,
                ni = self.mesh.numberOfFaces)
                
            return self._makeValue(value = val)
