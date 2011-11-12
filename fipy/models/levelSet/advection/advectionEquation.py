#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "advectionEquation.py"
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
 # ###################################################################
 ##
 
__docformat__ = 'restructuredtext'

from fipy.terms.transientTerm import TransientTerm

from fipy.models.levelSet.advection.advectionTerm import _AdvectionTerm

__all__ = ["buildAdvectionEquation"]

def buildAdvectionEquation(advectionCoeff=None,
                           advectionTerm=None):
    r"""

    The `buildAdvectionEquation` function constructs and returns an
    advection equation. The advection equation is given by:

    .. math::

       \frac{\partial \phi}{\partial t} + u \abs{\nabla \phi} = 0.

    This solution method for the `_AdvectionTerm` is set up specifically to
    evolve `var` while preserving `var` as a distance function. This
    equation is used in conjunction with the `DistanceFunction`
    object. Further details of the numerical method can be found in "Level
    Set Methods and Fast Marching Methods" by J.A. Sethian, Cambridge
    University Press, 1999. Testing for the advection equation is in
    :mod:`examples.levelSet.advection`

    :Parameters:
      - `advectionCoeff`: The coeff to pass to the `advectionTerm`.
      - `advectionTerm`: An advection term class.

    """
    if advectionTerm == None:
        advectionTerm = _AdvectionTerm

    return TransientTerm(1.) + advectionTerm(advectionCoeff)
        

        
        
        
        
