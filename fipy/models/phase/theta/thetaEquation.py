#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "thetaEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 9/3/04 {10:41:52 PM} { 4:13:57 PM}
 # 
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

from fipy.terms.transientTerm import TransientTerm
from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm

from sourceVariable import SourceVariable
from diffusionVariable import DiffusionVariable
from transientVariable import TransientVariable
from thetaHalfAngleVariable import ThetaHalfAngleVariable
from noModularVariable import NoModularVariable

def buildThetaEquation(phase = None,
                       theta = None,
                       parameters = {}):
	r"""
	Builds a crystal orientation equation of the following form,

	.. raw:: latex

	    $$ P(\epsilon | \nabla \theta |) \tau_{\theta} \phi^2 
	    \frac{\partial \theta}{\partial t} 
	    = \nabla \cdot \left[ \phi^2 \left( \frac{s}{| \nabla \theta |} 
	    + \epsilon^2 \right) \nabla \theta \right] $$
	    where
	    $$ P(w) = 1 - \exp{(-\beta w)} + \frac{\mu}{\epsilon} \exp{(-\beta w)}. $$
        

        :Parameters:
	  - `phase` : The phase field.
	  - `theta` : The crystal orientation.
	  - `parameters` : A dictionary with keys `'small value'`, `'beta'`,
	    `'mu'`, `'tau'`, `'gamma'`, `'epsilon'`, `'s'`, `'anisotropy'`, `'alpha'`,
	    `'symmetry'`.

	"""

        transientCoeff = TransientVariable(phase = phase, theta = theta, parameters = parameters)

        diffusionCoeff = DiffusionVariable(phase = phase, theta = theta, parameters = parameters)

        halfAngleVariable = ThetaHalfAngleVariable(phase = phase, theta = theta, parameters = parameters)
        
        sourceCoeff = SourceVariable(phase = phase,
                                     theta = theta,
                                     diffCoeff = diffusionCoeff,
                                     halfAngleVariable = halfAngleVariable,
                                     parameters = parameters)
        	
        transientTerm = TransientTerm(transientCoeff)
        diffusionTerm = ImplicitDiffusionTerm(diffusionCoeff)

        return transientTerm - diffusionTerm - sourceCoeff


