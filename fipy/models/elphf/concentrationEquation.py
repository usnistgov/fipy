"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "concentrationEquation.py"
                                   created: 11/12/03 {10:39:23 AM} 
                               last update: 12/29/03 {11:29:58 AM} 
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
from terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from terms.powerLawConvectionTerm import PowerLawConvectionTerm
from terms.centralDiffConvectionTerm import CentralDifferenceConvectionTerm

class ConcentrationEquation(MatrixEquation):
    """
    Diffusion equation is implicit.
    """    
    def __init__(self,
                 Cj,
		 timeStepDuration,
		 fields = {},
                 diffusivity = 1.,
		 convectionScheme = PowerLawConvectionTerm,
                 solver='default_solver',
                 boundaryConditions=()):
		     
        mesh = Cj.getMesh()
	
	diffusionTerm = ImplicitDiffusionTerm(
	    diffCoeff = diffusivity,
	    mesh = mesh,
	    boundaryConditions = boundaryConditions)
	    
	convectionTerm = convectionScheme(
	    convCoeff = self.getConvectionCoeff(Cj, fields, diffusivity),
	    mesh = mesh, 
	    boundaryConditions = boundaryConditions,
	    diffusionTerm = diffusionTerm)

	terms = (
	    TransientTerm(tranCoeff = 1. / timeStepDuration,mesh = mesh),
	    diffusionTerm,
	    convectionTerm
            )
	    
	MatrixEquation.__init__(
            self,
            var = Cj,
            terms = terms,
            solver = solver,
            solutionTolerance = 1e-10)

    def getConvectionCoeff(self, Cj, fields, diffusivity):
	Cj.pConvCoeff = diffusivity * Cj.getStandardPotential() * fields['phase'].get_p().getFaceGrad() 
# 	Cj.pConvCoeff = Cj.weightedDiffusivity * Cj.getStandardPotential() * 30 * fields['phase'].get_gFace().transpose() * fields['phase'].getFaceGrad()
	Cj.gConvCoeff = diffusivity * Cj.getBarrierHeight() * fields['phase'].get_g().getFaceGrad() 
	Cj.electromigrationCoeff = diffusivity * Cj.getValence() * fields['potential'].getFaceGrad() 
	
	return Cj.pConvCoeff + Cj.gConvCoeff + Cj.electromigrationCoeff

