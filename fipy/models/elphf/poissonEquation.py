#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "poissonEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 1/24/04 {11:59:37 PM} 
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


## from fivol.equations.matrixEquation import MatrixEquation
from fivol.equations.preRelaxationEquation import PreRelaxationEquation
from fivol.equations.postRelaxationEquation import PostRelaxationEquation
from fivol.equations.relaxationEquation import RelaxationEquation
from fivol.terms.transientTerm import TransientTerm
from fivol.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fivol.terms.scSourceTerm import ScSourceTerm
from fivol.terms.spSourceTerm import SpSourceTerm

from fivol.tools.dimensions import physicalField

from substitutionalSumVariable import SubstitutionalSumVariable

class PoissonEquation(PreRelaxationEquation):
    def __init__(self,
                 potential,
		 parameters,
		 fields = {},
                 solver='default_solver',
		 solutionTolerance = 1e-10,
		 relaxation = 1.,
                 boundaryConditions=()):
		     
        mesh = potential.getMesh()
	
	dielectric = physicalField.PhysicalField(value = "eps0 / (Faraday**2 * LENGTH**2 / (ENERGY * MOLARVOLUME))")
	# LENGTH, ENERGY, or MOLARVOLUME might not have correct units
	# in simple problems
	dielectric /= physicalField.PhysicalField(value = 1, unit = dielectric.inBaseUnits().getUnit())
	dielectric *= physicalField.Scale(parameters['dielectric'], "eps0") 
	
	diffusionTerm = ImplicitDiffusionTerm(
	    diffCoeff = dielectric,
	    mesh = mesh,
	    boundaryConditions = boundaryConditions)
	    
	fields['charge'] = fields['solvent'].getValence()
	for component in list(fields['interstitials']) + list(fields['substitutionals']):
	    fields['charge'] = fields['charge'] + self.getConcentration(component) * component.getValence() #.getOld()
	
## 	charge = charge * "1 Faraday/MOLARVOLUME"
	
	self.scTerm = ScSourceTerm(
	    sourceCoeff = fields['charge'],
	    mesh = mesh)
	    
	terms = (
	    diffusionTerm,
	    self.scTerm
	)
	    
	PreRelaxationEquation.__init__(
            self,
            var = potential,
            terms = terms,
            solver = solver,
	    solutionTolerance = solutionTolerance,
	    relaxation = relaxation)
	    
    def getConcentration(self, component):
	return component
	

