#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "concentrationEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 9/3/04 {10:29:48 PM} 
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

## from fipy.equations.matrixEquation import MatrixEquation
from fipy.equations.preRelaxationEquation import PreRelaxationEquation
from fipy.equations.postRelaxationEquation import PostRelaxationEquation
from fipy.equations.relaxationEquation import RelaxationEquation
from fipy.terms.transientTerm import TransientTerm
from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fipy.terms.powerLawConvectionTerm import PowerLawConvectionTerm
from fipy.terms.centralDiffConvectionTerm import CentralDifferenceConvectionTerm

class ConcentrationEquation(RelaxationEquation):
    def __init__(self,
                 Cj,
		 fields = {},
		 convectionScheme = PowerLawConvectionTerm,
                 solver='default_solver',
		 relaxation = 0.,
		 phaseRelaxation = 1.,
		 solutionTolerance = 1e-10,
                 boundaryConditions=()):
		     
	self.phaseRelaxation = phaseRelaxation
		     
        mesh = Cj.getMesh()
	
	diffusionTerm = ImplicitDiffusionTerm(
	    diffCoeff = Cj.getDiffusivity(),
	    mesh = mesh,
	    boundaryConditions = boundaryConditions)
        
	convectionTerm = convectionScheme(
	    convCoeff = self.getConvectionCoeff(Cj, fields, relaxation),
	    mesh = mesh, 
	    boundaryConditions = boundaryConditions,
	    diffusionTerm = diffusionTerm)
	    
	terms = (
	    TransientTerm(tranCoeff = 1., mesh = mesh),
	    diffusionTerm,
	    convectionTerm
	)
	
	RelaxationEquation.__init__(
            self,
            var = Cj,
            terms = terms,
            solver = solver,
	    solutionTolerance = solutionTolerance,
	    relaxation = relaxation)
	    
    def getConvectionCoeff(self, Cj, fields, diffusivity = None):
	if diffusivity is None:
	    diffusivity = Cj.getDiffusivity()
	    
## 	diffusivity = diffusivity / "1 ENERGY"
## 	Cj.pConvCoeff = diffusivity * Cj.getStandardPotential() * fields['phase'].get_p().getFaceGrad() 
## 	Cj.gConvCoeff = diffusivity * Cj.getBarrierHeight() * fields['phase'].get_g().getFaceGrad() 
	Cj.pConvCoeff = diffusivity * Cj.getStandardPotential() * fields['phase'].get_pPrime().getArithmeticFaceValue().transpose()
	Cj.gConvCoeff = diffusivity * Cj.getBarrierHeight() * fields['phase'].get_gPrime().getArithmeticFaceValue().transpose()
## 	Cj.electromigrationCoeff = diffusivity * "1 Faraday" * Cj.getValence() * fields['potential'].getFaceGrad() 
	Cj.electromigrationCoeff = diffusivity * Cj.getValence() * fields['potential'].getFaceGrad() 
	
## 	self.phaseRelaxation
	
## 	return (Cj.pConvCoeff + Cj.gConvCoeff) + Cj.electromigrationCoeff
	return (Cj.pConvCoeff + Cj.gConvCoeff) * fields['phase'].getFaceGrad() + Cj.electromigrationCoeff

