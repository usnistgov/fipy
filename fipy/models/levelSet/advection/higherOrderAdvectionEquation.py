#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "higherOrderAdvectionEquation.py"
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
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.models.levelSet.advection.higherOrderAdvectionTerm import _HigherOrderAdvectionTerm
from fipy.models.levelSet.advection.advectionEquation import buildAdvectionEquation

__all__ = ["buildHigherOrderAdvectionEquation"]

def buildHigherOrderAdvectionEquation(advectionCoeff = None):
    r"""

    The `buildHigherOrderAdvectionEquation` function returns an
    advection equation that uses the `_HigherOrderAdvectionTerm`. The
    advection equation is given by,

    .. math::

       \frac{\partial \phi}{\partial t} + u \abs{\nabla \phi} = 0.

    :Parameters:
      - `advectionCoeff`: The `coeff` to pass to the `_HigherOrderAdvectionTerm`

    """
    return buildAdvectionEquation(advectionCoeff = advectionCoeff, advectionTerm = _HigherOrderAdvectionTerm)
