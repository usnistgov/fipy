"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "phaseEquation.py"
                                   created: 11/12/03 {10:39:23 AM} 
                               last update: 12/10/03 {12:06:40 AM} 
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
from phaseMVariable import PhaseMVariable
from phaseSpSourceVariable import PhaseSpSourceVariable
from phaseScSourceVariable import PhaseScSourceVariable
from phaseDiffCoeffVariable import PhaseDiffCoeffVariable
import Numeric

class PhaseEquation(MatrixEquation):
    """
    Diffusion equation is implicit.
    """    
    def __init__(self,
                 var,
                 solver = 'default_solver',
                 boundaryConditions = (),
                 parameters = {}):
        
        mesh = var.getMesh()
	
	parameters['mPhi'] = PhaseMVariable(
	    mesh = mesh,
	    phi = parameters['phi'],
	    temperature = parameters['temperature'])
	
        self.diffTerm = ExplicitDiffusionTerm(
	    diffCoeff = PhaseDiffCoeffVariable(
		mesh = mesh,
		parameters = parameters),
	    mesh = mesh,
	    boundaryConditions = boundaryConditions)
	    
        self.spTerm = SpSourceTerm(
	    sourceCoeff = PhaseSpSourceVariable(
		mesh = mesh, 
		parameters = parameters),
	    mesh = mesh)
	    
        self.scTerm = ScSourceTerm(
	    sourceCoeff = PhaseScSourceVariable(
		mesh = mesh, 
		parameters = parameters),
	    mesh = mesh)
	
	transientCoeff = parameters['tau']
	terms = (
	    TransientTerm(transientCoeff,mesh),
	    self.diffTerm,
            self.scTerm,
            self.spTerm
	)

	MatrixEquation.__init__(
            self,
            var,
            terms,
            solver)
