#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "phaseEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 9/3/04 {10:35:49 PM} 
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

from fipy.terms.transientTerm import TransientTerm
from fipy.terms.explicitDiffusionTerm import ExplicitDiffusionTerm
from fipy.terms.dependentSourceTerm import DependentSourceTerm

from phaseDiffusionVariable import PhaseDiffusionVariable
from anisotropyVariable import AnisotropyVariable
from spSourceVariable import SpSourceVariable
from phaseHalfAngleVariable import PhaseHalfAngleVariable
from scSourceVariable import ScSourceVariable

def buildPhaseEquation(phase = None, theta = None, temperature = None, parameters = {}, mPhi = None):

    thetaOld = theta.getOld()
        	
    mPhiVar = mPhi(phase = phase, temperature = temperature, parameters = parameters)

    halfAngle = PhaseHalfAngleVariable(
        parameters = parameters,
        phase = phase,
        theta = thetaOld
        )

    diffTerm = ExplicitDiffusionTerm(
        diffCoeff = PhaseDiffusionVariable(parameters = parameters,
                                           halfAngle = halfAngle))
	
    spTerm = DependentSourceTerm(
        sourceCoeff = SpSourceVariable(theta = thetaOld,
                                       mPhi = mPhiVar,
                                       phase = phase,
                                       parameters = parameters))

    anisotropy = AnisotropyVariable(parameters = parameters, phase = phase, halfAngle = halfAngle)

    sourceCoeff = ScSourceVariable(mPhi = mPhiVar,
                                   phase = phase,
                                   anisotropy = anisotropy)

    transientCoeff = parameters['tau']
    

    return TransientTerm(tranCoeff = transientCoeff) - diffTerm  + spTerm - sourceCoeff
