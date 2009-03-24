#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cellToFaceVariable.py"
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

from fipy.variables.faceVariable import FaceVariable
from fipy.tools import inline
from fipy.tools import numerix

class _CellToFaceVariable(FaceVariable):
    def __init__(self, var):
        FaceVariable.__init__(self, mesh=var.getMesh(), elementshape=var.shape[:-1])
        self.var = self._requires(var)

    def _calcValue(self):
        alpha = self.mesh._getFaceToCellDistanceRatio()
        id1, id2 = self.mesh._getAdjacentCellIDs()
        return inline._optionalInline(self._calcValueIn, self._calcValuePy, alpha, id1, id2)

    def __getstate__(self):
        return {
            'var': self.var
        }
        
    def __setstate__(self, dict):
        self.__init__(**dict)

