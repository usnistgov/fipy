#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "temperatureEquation.py"
 #                                    created: 12/17/03 {10:39:23 AM} 
 #                                last update: 9/3/04 {10:35:59 PM} 
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

def buildTemperatureEquation(phase, parameters = {}):
    r"""
    Creates a temperature equation of the following form,
    
    .. raw:: latex

        $$ \frac{\partial T}{\partial t} = D_T \nabla^2 T + \frac{\partial \phi}{\partial t}. $$

    :Parameters:
      - `phase` : The phase field.
      - `parameters` : A dictionary with keys `'latent heat'`, `'heat capacity'`,
        `'timeStepDuration'` and `'temperature diffusion'`.

    """
    latentHeat = parameters['latent heat']
    heatCapacity = parameters['heat capacity']
    phaseOld =  phase.getOld()
    dt = parameters['timeStepDuration']
        
    return TransientTerm(1.) - ImplicitDiffusionTerm(parameters['temperature diffusion']) - latentHeat / heatCapacity * (phase - phaseOld) / dt
