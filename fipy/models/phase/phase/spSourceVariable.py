#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "spSourceVariable.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 9/3/04 {10:33:21 PM} 
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

class SpSourceVariable(CellVariable):
    def __init__(self, theta = None, mPhi = None, phase = None, parameters = None):
        CellVariable.__init__(self, mesh = theta.getMesh())
        self.theta = self.requires(theta)
        self.phase = self.requires(phase)                                   
        self.parameters = parameters
        self.mPhi = self.requires(mPhi)

    def _calcValue(self):
        inline.optionalInline(self._calcValueIn, self._calcValuePy)
    
    def  _calcValuePy(self):
        s = self.parameters['s']
	epsilon = self.parameters['epsilon']

        thetaMag = self.theta.getGrad().getMag()[:]
        
        spSourceCoeff = self.mPhi[:] * (self.phase[:] - (self.mPhi[:] < 0.))
        spSourceCoeff += (2*s + epsilon**2 * thetaMag) * thetaMag

	self.value = spSourceCoeff

    def _calcValueIn(self):
        inline.runInlineLoop1("""
            double tmp = 0.;
            if(mPhi(i) < 0.)
              tmp = 1.;
            value(i) = mPhi(i) * (phase(i) - tmp);
            value(i) += (2 * s + epsilonSq * thetaMag(i)) * thetaMag(i);
            """,mPhi = self.mPhi.getNumericValue(),
                phase = self.phase.getNumericValue(),
                value = self._getArray(),
                s = self.parameters['s'],
                epsilonSq = self.parameters['epsilon']**2,
                thetaMag = self.theta.getGrad().getMag().getNumericValue(),
                ni = len(self._getArray()))
