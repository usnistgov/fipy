#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "arithmeticCellToFaceVariable.py"
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

__all__ = []

from fipy.variables.cellToFaceVariable import _CellToFaceVariable
from fipy.tools import numerix
from fipy.tools import inline

class _ArithmeticCellToFaceVariable(_CellToFaceVariable):
    if inline.doInline:
        def _calcValue_(self, alpha, id1, id2):
            val = self._array.copy()
            
            inline._runIterateElementInline("""
                int ID1 = ITEM(id1, i, NULL);
                int ID2 = ITEM(id2, i, NULL);
                double cell1 = ITEM(var, ID1, vec);
                double cell2 = ITEM(var, ID2, vec);
                ITEM(val, i, vec) = (cell2 - cell1) * ITEM(alpha, i, NULL) + cell1;
            """,
            var = self.var.numericValue,
            val = val, 
            alpha = alpha,
            id1 = id1, id2 = id2,
            shape=numerix.array(numerix.shape(val)),
            ni = self.mesh.numberOfFaces)

            return self._makeValue(value = val)
    else:
        def _calcValue_(self, alpha, id1, id2):
            cell1 = numerix.take(self.var, id1, axis=-1)
            cell2 = numerix.take(self.var, id2, axis=-1)
            return (cell2 - cell1) * alpha + cell1

        
