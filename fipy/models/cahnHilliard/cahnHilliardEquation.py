#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cahnHilliardEquation.py"
 #                                    created: 7/2/04 {10:39:23 AM} 
 #                                last update: 9/3/04 {10:40:23 PM} 
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

"""

The Cahn-Hilliard equation is given by:

.. raw:: latex

    $$ \\frac{\\partial \\phi}{\\partial t} = \\nabla \\cdot D \\nabla \\left( \\frac{\\partial f}{\\partial \\phi} - \\epsilon^2 \\nabla^2 \\phi \\right) $$

where the free energy functional is given by,

.. raw:: latex

    $$ f = \\frac{a^2}{2} \\phi^2 (1 - \\phi)^2 $$

In the `CahnHilliardEquation` object the equation is transformed into
the following form,

.. raw:: latex

    $$ \\frac{\\partial \\phi}{\\partial t} = \\nabla \\cdot D \\frac{\\partial^2 f}{\\partial \\phi^2} \\nabla \\phi - \\nabla \\cdot D \\nabla \\epsilon^2 \\nabla^2 \\phi $$

This form of the equation allows the `CahnHilliardEquation` to be
constructed from a transient term, a diffusion term, and a fourth
order diffusion term. Notice that the diffusion coefficient for the
diffusion term does not always remain positive since,

.. raw:: latex

    $$ \\frac{\\partial^2 f}{\\partial \\phi^2} = a^2 (1 - 6 \\phi (1 - \\phi)) $$

can be less than zero and thus unstable. The fourth order diffusion
term acts to stabilize the problem. 

"""
__docformat__ = 'restructuredtext'

import Numeric

from fipy.equations.matrixEquation import MatrixEquation
from fipy.terms.transientTerm import TransientTerm
from fipy.terms.nthOrderDiffusionTerm import NthOrderDiffusionTerm

class CahnHilliardEquation(MatrixEquation):

    def __init__(self,
                 var,
                 solver = 'default_solver',
                 transientCoeff = 1.0,
                 boundaryConditions = (),
                 parameters = {}):
        
	self.parameters = parameters
	self.var = var
        diffusionCoeff = self.parameters['diffusionCoeff']

        faceVar = var.getArithmeticFaceValue()
        doubleWellDerivative = self.parameters['asq'] * ( 1 - 6 * faceVar * (1 - faceVar))
        
	terms = (
	    TransientTerm(
                transientCoeff,
                mesh = var.getMesh()
                )
            ,
	    NthOrderDiffusionTerm(
                coeffs = (diffusionCoeff, -self.parameters['epsilon']**2),
                mesh = var.getMesh(),
                boundaryConditions = boundaryConditions
                )
            ,
            NthOrderDiffusionTerm(
                coeffs = (diffusionCoeff * doubleWellDerivative,),
                mesh = var.getMesh(),
                boundaryConditions = boundaryConditions
                )
            )

	MatrixEquation.__init__(
            self,
            var,
            terms,
            solver)
