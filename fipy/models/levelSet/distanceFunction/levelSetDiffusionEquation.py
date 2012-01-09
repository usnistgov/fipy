#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "levelSetDiffusionEquation.py"
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

__all__ = []

from fipy.models.levelSet.distanceFunction.levelSetDiffusionVariable import _LevelSetDiffusionVariable
from fipy.terms.transientTerm import TransientTerm
from fipy.terms.diffusionTerm import DiffusionTermNoCorrection

def _buildLevelSetDiffusionEquation(ionVar = None,
                                   distanceVar = None,
                                   transientCoeff = 1.,
                                   diffusionCoeff = 1.):
    r"""

    The `LevelSetDiffusionEquation` solves the diffusion of a species in
    conjunction with the level set equation. Essentially the species is
    only transported in the electrolyte. The governing equation is given
    by,

    .. math::

       \frac{\partial c}{\partial t} = \nabla \cdot D \nabla  c

    where,

    .. math::

       D &= \begin{cases}
           D_c & \text{when $\phi > 0$} \\
           0 & \text{when $\phi \le 0$}
      \end{cases}

    :Parameters:    
      - `ionVar` : The species concentration variable.    
      - `distanceVar` : A `DistanceVariable` object.
      - `transientCoeff` : The coefficient for the `TransientTerm`
      - `diffusionCoeff` : The coefficient for the `DiffusionTerm`

    """

    diffusionCoeff = _LevelSetDiffusionVariable(distanceVar,
                                               diffusionCoeff)
        
    return TransientTerm(transientCoeff) - DiffusionTermNoCorrection(diffusionCoeff)
        
def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test() 
