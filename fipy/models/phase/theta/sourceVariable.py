"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "sourceVariable.py"
                                   created: 11/12/03 {10:39:23 AM} 
                               last update: 01/07/04 { 4:30:43 PM}
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

from variables.cellVariable import CellVariable

class SourceVariable(CellVariable):

    def __init__(self,
                 phase = None,
                 theta = None,
                 diffCoeff = None,
                 halfAngleVariable = None,
                 parameters = None):

        CellVariable.__init__(self, phase.getMesh())

        self.parameters = parameters
        self.phase = self.required(phase)
        self.theta = self.required(theta)
        self.diffCoeff = self.required(diffCoeff)
        self.halfAngleVariable = self.required(halfAngleVariable)

    def calcValue(self):

        mesh = self.phase.getMesh
        c2 = self.['anisotropy']
        
        thetaGradDiff = self.theta.getFaceGrad()[:] - self.theta.getFaceGradNoMod()[:]

        correctionTerm = tools.addOverFaces(faceGradient = thetaGradDiff,
                                            faceVariable = self.diffCoeff[:],
                                            mesh = mesh,
                                            NCells = len(self.phase[:]))

        halfAngleSq = self.halfAngleVariable[:] * self.halfAngleVariable[:]
        beta = (1. - halfAngleSq) / (1. + halfAngleSq)
        dbeta = N * 2. * self.halfAngleVariable[:] / (1. + halfAngleSq)
        
        self.value = correctionTerm
        self.value += (self.parameters['alpha']**2 * self.parameters['anisotropy'] * dbeta *
                       self.phase.getGradMag() * (1. + c2 * beta))
        
    

