#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "phaseHalfAngleVariable.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 4/1/05 {11:02:45 AM}
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
from fipy.variables.faceVariable import FaceVariable

import fipy.tools.array as array

class PhaseHalfAngleVariable(FaceVariable):
    def __init__(self, parameters = None, phase = None, theta = None):
        FaceVariable.__init__(self, phase.getMesh())
	self.parameters = parameters
	self.phase = self._requires(phase)
        self.theta = self._requires(theta)
        
    def _calcValue(self):
        inline._optionalInline(self._calcValueInline, self._calcValuePy)

    def _calcValuePy(self):
	N = self.parameters['symmetry']
        dphi = self.phase.getFaceGrad()[:,:]
        thetaFace = self.theta.getArithmeticFaceValue()[:]
	z = array.arctan2(dphi[:,1],dphi[:,0])
        z = N * (z - thetaFace)
        self.value = array.tan(z / 2.)

    def _calcValueInline(self):
        inline._runInlineLoop1("""
        z = atan2(dphi(i,1), dphi(i,0));
        z = symmetry * (z - thetaFace(i));
        value(i) = tan(z / 2.);""",z = 0.,
                              dphi = self.phase.getFaceGrad().getNumericValue(),
                              symmetry = self.parameters['symmetry'],
                              thetaFace = self.theta.getArithmeticFaceValue().getNumericValue(),
                              value = self._getArray(),
                              ni = len(self._getArray()))                              
            


