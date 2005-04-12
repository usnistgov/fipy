#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "phaseEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 2/18/05 {3:17:07 PM} 
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
from fipy.terms.explicitDiffusionTerm import ExplicitDiffusionTerm
from fipy.terms.implicitSourceTerm import ImplicitSourceTerm

from phaseDiffusionVariable import _PhaseDiffusionVariable
from anisotropyVariable import _AnisotropyVariable
from spSourceVariable import _SpSourceVariable
from phaseHalfAngleVariable import _PhaseHalfAngleVariable
from scSourceVariable import _ScSourceVariable

def buildPhaseEquation(phase = None, theta = None, temperature = None, parameters = {}, mPhi = None):
    r"""
    Creates a phase field equation of the following form,

    .. raw:: latex
    
        $$ \tau_{\phi} \frac{\partial \phi}{\partial t} = \nabla \cdot
        \left[ D \nabla \phi + A \nabla \xi \right] + \phi ( 1 - \phi
        ) m_1 ( \phi , T) - 2 s \phi | \nabla \theta | - \epsilon^2
        \phi | \nabla \theta |^2. $$

        The coefficients $D$ and $A$ are given by,

        $$ D = \alpha^2 \left[ 1 + c \beta \right]^2 $$

        and

        $$ A = \alpha^2 c \left[ 1 + c \beta \right] \Phi_\psi $$

        where $ \beta = \frac{ 1 - \Phi^2 } { 1 + \Phi^2} $,
        $ \Phi = \tan \left( \frac{ \theta } { 2 } + \frac{ N } { 2 } \arctan \psi \right) $,
        $ \psi = \frac{ \phi_y } { \phi_x } $ and
        $ \xi_x = -\phi_y $ and $ \xi_y = \phi_x $.



    :Parameters:
      - `phase`: The phase field.
      - `theta`: The crystal orientation.
      - `temperature`: The system temperature.
      - `parameters`: A dictionary that includes the following keys,
        `'tau'`, `'epsilon'`, `'s'`, `'anisotropy'`, `'alpha'`, `'symmetry'`, `'c2'`. 
        
    """
            
      

    thetaOld = theta.getOld()
        	
    mPhiVar = mPhi(phase = phase, temperature = temperature, parameters = parameters)

    halfAngle = _PhaseHalfAngleVariable(
        parameters = parameters,
        phase = phase,
        theta = thetaOld
        )

    diffTerm = ExplicitDiffusionTerm(
        coeff = _PhaseDiffusionVariable(parameters = parameters,
                                           halfAngle = halfAngle))
	
    spTerm = ImplicitSourceTerm(
        coeff = _SpSourceVariable(theta = thetaOld,
                                       mPhi = mPhiVar,
                                       phase = phase,
                                       parameters = parameters))

    anisotropy = _AnisotropyVariable(parameters = parameters, phase = phase, halfAngle = halfAngle)

    sourceCoeff = _ScSourceVariable(mPhi = mPhiVar,
                                   phase = phase,
                                   anisotropy = anisotropy)

    transientCoeff = parameters['tau']
    

    return TransientTerm(coeff = transientCoeff) - diffTerm  + spTerm - sourceCoeff
