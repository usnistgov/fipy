"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "phaseHalfAngleVariable.py"
                                   created: 11/12/03 {10:39:23 AM} 
                               last update: 01/05/04 { 5:40:18 PM}
 Author: Jonathan Guyer
 E-mail: guyer@nist.gov
 Author: Daniel Wheeler
 E-mail: daniel.wheeler@nist.gov
   mail: NIST
    www: http://ctcms.nist.gov
 
========================================================================
This software was developed at the National Institute of Standards
and Technology by employees of the Federal Government in the course
of their official duties.  Pursuant to title 17 Section 105 of the
United States Code this software is not subject to copyright
protection and is in the public domain.  PFM is an experimental
system.  NIST assumes no responsibility whatsoever for its use by
other parties, and makes no guarantees, expressed or implied, about
its quality, reliability, or any other characteristic.  We would
appreciate acknowledgement if the software is used.

This software can be redistributed and/or modified freely
provided that any derivative works bear some notice that they are
derived from it, and any modified versions bear some notice that
they have been modified.
========================================================================
 
 Description: 

 History

 modified   by  rev reason
 ---------- --- --- -----------
 2003-11-12 JEG 1.0 original
###################################################################
"""

from variables.faceVariable import FaceVariable
import Numeric

class PhaseHalfAngleVariable(FaceVariable):
    def __init__(self, parameters = None, phase = None, theta = None):
        FaceVariable.__init__(self, phase.getMesh())
	self.parameters = parameters
	self.phase = self.requires(phase)
        self.theta = self.requires(theta)

    def calcValue(self):
	N = self.parameters['symmetry']
        dphi = self.phase.getFaceGrad()[:,:]
        thetaFace = self.theta.getFaceValue()[:]

        z = Numeric.arctan2(dphi[:,1],dphi[:,0])
	z = N * (z - thetaFace)
	self.value = Numeric.tan(z / 2.)

