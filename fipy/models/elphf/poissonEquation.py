#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "poissonEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 12/8/04 {5:12:52 PM} 
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
__docformat__ = 'restructuredtext'

from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fipy.tools.dimensions import physicalField

from equationFactory import EquationFactory

class PoissonEquationFactory(EquationFactory):
    def make(self, fields, parameters):
	from elphf import constant

	permittivity = physicalField.Scale(parameters['permittivity'], 
					   constant['Faraday']**2 * constant['LENGTH']**2 
					   / (constant['ENERGY'] * constant['MOLARVOLUME'])) 

	return ImplicitDiffusionTerm(diffCoeff = permittivity) + fields['charge']

factory = PoissonEquationFactory()
						   
