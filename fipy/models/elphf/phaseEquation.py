"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "phaseEquation.py"
                                   created: 11/12/03 {10:39:23 AM} 
                               last update: 12/23/03 {6:21:43 PM} 
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

from equations.matrixEquation import MatrixEquation
from terms.transientTerm import TransientTerm
from substitutionalSumVariable import SubstitutionalSumVariable
from terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from terms.scSourceTerm import ScSourceTerm
from terms.spSourceTerm import SpSourceTerm
import Numeric

class PhaseEquation(MatrixEquation):
    def __init__(self,
                 phase,
		 timeStepDuration,
		 fields = {},
                 phaseMobility = 1.,
		 phaseGradientEnergy = 1.,
                 solver='default_solver',
                 boundaryConditions=()):
		     
        mesh = phase.getMesh()
	
	diffusionTerm = ImplicitDiffusionTerm(
	    diffCoeff = phaseMobility * phaseGradientEnergy,
	    mesh = mesh,
	    boundaryConditions = boundaryConditions)
	    
	enthalpy = fields['solvent'].getStandardPotential()
	barrier = fields['solvent'].getBarrierHeight()
	
	for component in fields['substitutionals']:
	    enthalpy = enthalpy + component * component.getStandardPotential() #.getOld()
	    barrier = barrier + component * component.getBarrierHeight() #.getOld()
	    
	self.mPhi = -phaseMobility * (30. * phase * (1. - phase) * enthalpy + (1. - 2 * phase) * barrier)
	
	self.spTerm = SpSourceTerm(
	    sourceCoeff = self.mPhi * (phase - (self.mPhi < 0.)),
	    mesh = mesh)
	    
	self.scTerm = ScSourceTerm(
	    sourceCoeff = (self.mPhi > 0.) * self.mPhi * phase,
	    mesh = mesh)
	    
	terms = (
	    TransientTerm(tranCoeff = 1. / timeStepDuration,mesh = mesh),
	    diffusionTerm,
	    self.scTerm,
	    self.spTerm
            )
	    
	MatrixEquation.__init__(
            self,
            var = phase,
            terms = terms,
            solver = solver,
            solutionTolerance = 1e-10)

