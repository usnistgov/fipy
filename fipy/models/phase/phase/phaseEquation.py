"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "phaseEquation.py"
                                   created: 11/12/03 {10:39:23 AM} 
                               last update: 12/22/03 {5:53:13 PM} 
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
from terms.explicitDiffusionTerm import ExplicitDiffusionTerm
from terms.scSourceTerm import ScSourceTerm
from terms.spSourceTerm import SpSourceTerm
from phaseDiffusionVariable import PhaseDiffusionVariable
from anisotropyVariable import AnisotropyVariable
from spSourceVariable import SpSourceVariable
import Numeric

class PhaseEquation(MatrixEquation):
    """
    Diffusion equation is implicit.
    """    
    def __init__(self,
                 var,
                 solver = 'default_solver',
                 boundaryConditions = (),
                 fields = {},
                 parameters = {},
                 mPhi = 0.):
        
        mesh = var.getMesh()
	
	self.parameters = parameters

	self.var = var
	self.temp = fields['temperature']
        self.thetaOld = fields['theta'].getOld()
        	
	self.mPhi = mPhi

        self.diffTerm = ExplicitDiffusionTerm(
	    diffCoeff = PhaseDiffusionVariable(self.parameters, self.var, self.thetaOld),
##            diffCoeff = self.getPhaseDiffCoeff(),
	    mesh = mesh,
	    boundaryConditions = boundaryConditions)
	
        self.spTerm = SpSourceTerm(
##	    sourceCoeff = self.getSpSourceCoeff(),
            sourceCoeff = SpSourceVariable(theta = self.thetaOld, mPhi = self.mPhi, phase = self.var, parameters = self.parameters),
	    mesh = mesh)

        anisotropy = AnisotropyVariable(parameters, self.var, self.thetaOld)

        self.scTerm = ScSourceTerm(
            sourceCoeff = (self.mPhi > 0.) * self.mPhi * self.var + anisotropy,
	    mesh = mesh)
	
	transientCoeff = parameters['tau'] / parameters['time step duration']
        
	terms = (
	    TransientTerm(transientCoeff, mesh),
	    self.diffTerm,
            self.scTerm,
            self.spTerm
	)

	MatrixEquation.__init__(
            self,
            var,
            terms,
            solver)

    
	    
##    def getPhaseDiffCoeff(self):
##	alpha = self.parameters['alpha']
##	N = self.parameters['symmetry']
##	c2 = self.parameters['anisotropy']
##	dphi = self.var.getFaceGrad()
##	thetaFace = self.thetaOld.getFaceValue()

##        z = Numeric.arctan2(dphi[:,1],dphi[:,0])
##	z = N * (z - thetaFace)
##	z = Numeric.tan(z / 2.)
##	z = z * z
##	z = (1. - z) / (1. + z)
##	z = (1.+ c2 * z)

##	return alpha**2 * z * z

    def getSpSourceCoeff(self):
        s = self.parameters['s']
	epsilon = self.parameters['epsilon']
	thetaMag = self.thetaOld.getGrad().getMag()
        
        spSourceCoeff = self.mPhi * (self.var - (self.mPhi < 0.)) 
	spSourceCoeff += (2*s + epsilon**2 * thetaMag) * thetaMag
	
	return spSourceCoeff
