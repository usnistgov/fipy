#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "poissonEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 11/1/04 {11:33:36 AM} 
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

__docformat__ = 'restructuredtext'

import Numeric
 
## from fipy.equations.matrixEquation import MatrixEquation
from fipy.equations.preRelaxationEquation import PreRelaxationEquation
from fipy.equations.postRelaxationEquation import PostRelaxationEquation
from fipy.equations.relaxationEquation import RelaxationEquation
from fipy.terms.transientTerm import TransientTerm
from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fipy.terms.scSourceTerm import ScSourceTerm
from fipy.terms.spSourceTerm import SpSourceTerm

from fipy.tools.dimensions import physicalField

class PoissonEquation(RelaxationEquation):
    r"""
    Represents the Poisson equation
    
    .. raw:: latex
    
       \[ 
	   \underbrace{
	       \nabla\cdot\left(\epsilon\nabla\phi\right) 
	   }_{\text{diffusion}}
	   +
	   \underbrace{
	       \rho
	   }_{\text{source}}
	   = 0
       \]

       where \( \phi \) is the electrostatic potential, 
       \( \epsilon \) is the dielectric constant
       \( \rho \equiv \sum_{j=1}^n z_j C_j \), is the total charge,
       \( C_j \) is the concentration of the \( j^\text{th} \)
       species, and \( z_j \) is the valence of that species.
    """
    def __init__(self,
                 potential,
		 parameters,
		 fields = {},
                 solver='default_solver',
		 solutionTolerance = 1e-10,
		 relaxation = 1.,
                 boundaryConditions=()):
		     
        mesh = potential.getMesh()
	
	from elphf import constant as k
	
	permittivity = physicalField.Scale(parameters['permittivity'], k['Faraday']**2 * k['LENGTH']**2 / (k['ENERGY'] * k['MOLARVOLUME'])) 
	
	diffusionTerm = ImplicitDiffusionTerm(
	    diffCoeff = permittivity,
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
	    
	RelaxationEquation.__init__(
            self,
            var = potential,
            terms = terms,
            solver = solver,
	    solutionTolerance = solutionTolerance,
	    relaxation = relaxation)
	    
    def getConcentration(self, component):
	return component
	
##     def getResidual(self):
## 	return Numeric.array((1e-16,))
	

