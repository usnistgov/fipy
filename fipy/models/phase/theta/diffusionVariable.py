#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "diffusionVariable.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 2/20/04 {2:05:32 PM}
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
from fipy.variables.faceVariable import FaceVariable

class DiffusionVariable(FaceVariable):

    def __init__(self, phase = None, theta = None, parameters = None):

        FaceVariable.__init__(self, phase.getMesh())

        self.parameters = parameters
        self.phase = self.requires(phase)
        self.theta = self.requires(theta)

    def calcValue(self):
        inline.optionalInline(self._calcValueInline, self._calcValue)

    def _calcValue(self):

        gamma = self.parameters['gamma']

        
        phaseFace = self.phase.getArithmeticFaceValue()[:]
        phaseSq = phaseFace * phaseFace
        gradMag = self.theta.getFaceGrad().getMag()[:]
        IGamma = Numeric.where(gradMag > 1. / gamma, 1 / gradMag, gamma)
## 	IGamma = (gradMag > 1. / gamma) * (1 / gradMag - gamma) + gamma

        self.value = phaseSq * (self.parameters['s'] * IGamma + self.parameters['epsilon']**2)
    
    def _calcValueInline(self):

        inline.runInlineLoop1("""
        phaseSq = phaseFace(i) * phaseFace(i);
        IGamma = gamma;
        if(gradMag(i) > 1. / gamma)
          IGamma = 1 / gradMag(i);
        value(i) = phaseSq * (s * IGamma + epsilon * epsilon);""",
                              phaseSq = 0.,
                              phaseFace = self.phase.getArithmeticFaceValue().getNumericValue(),
                              gradMag = self.theta.getFaceGrad().getMag().getNumericValue(),
                              gamma = self.parameters['gamma'],
                              IGamma = 0.,
                              s = self.parameters['s'],
                              epsilon = self.parameters['epsilon'],
                              value = self.value.value,
                              ni = self.mesh.getNumberOfFaces())
                              
                              
