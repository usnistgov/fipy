#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "scSourceVariable.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 7/24/04 {9:01:50 AM} 
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
from fipy.tools.inline import inline
from fipy.variables.cellVariable import CellVariable

class ScSourceVariable(CellVariable):
    def __init__(self, mPhi = None, phase = None, anisotropy = None):
        CellVariable.__init__(self, mesh = phase.getMesh())
        self.mPhi = self.requires(mPhi)
        self.phase = self.requires(phase)
        self.anisotropy = self.requires(anisotropy)

    def _calcValue(self):
        inline.optionalInline(self._calcValueIn, self._calcValuePy)
    
    def _calcValuePy(self):
        self.value = (self.mPhi[:] > 0.) * self.mPhi[:] * self.phase[:] + self.anisotropy[:]

    def _calcValueIn(self):
        inline.runInlineLoop1("""
            if (mPhi(i) > 0.)
                value(i) = mPhi(i) * phase(i) + anisotropy(i);                
            else
                value(i) = anisotropy(i);
        """,mPhi = self.mPhi.getNumericValue(),
            phase = self.phase.getNumericValue(),
            anisotropy = self.anisotropy.getNumericValue(),
            value = self._getArray(),
            ni = len(self._getArray()))
