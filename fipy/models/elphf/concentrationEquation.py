#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "concentrationEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 2/18/05 {3:06:08 PM} 
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
from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fipy.terms.powerLawConvectionTerm import PowerLawConvectionTerm

from equationFactory import EquationFactory

class ConcentrationEquationFactory(EquationFactory):
    def make(self, Cj, fields, convectionScheme = PowerLawConvectionTerm):
				    
	diffusionTerm = ImplicitDiffusionTerm(coeff = Cj.getDiffusivity())

	return TransientTerm() - diffusionTerm \
	    - convectionScheme(coeff = self.getConvectionCoeff(Cj, fields), 
			       diffusionTerm = diffusionTerm)
	    
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
	
## 	return (Cj.pConvCoeff + Cj.gConvCoeff) + Cj.electromigrationCoeff
	return (Cj.pConvCoeff + Cj.gConvCoeff) * fields['phase'].getFaceGrad() + Cj.electromigrationCoeff

factory = ConcentrationEquationFactory()