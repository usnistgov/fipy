#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "levelSetDiffusionEquation.py"
 #                                    created: 9/8/04 {10:39:23 AM} 
 #                                last update: 9/8/04 {4:00:40 PM} 
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

The `LevelSetDiffusionEquation` solves the diffusion of a species in
conjunction with the level set equation. Essentially the spcies is
only transported in the electrolyte. The governing equation is given
by,

.. raw:: latex

    $$ \\frac{\\partial c}{\\partial t} = \\nabla \\cdot D \\nabla  c $$

where,

.. raw:: latex

    $$ D = D_c \\;\\; \\text{when} \\;\\; \\phi > 0 $$
    $$ D = 0   \\;\\; \\text{when} \\;\\; \\phi \\le 0 $$
   
"""

__docformat__ = 'restructuredtext'


from fipy.equations.diffusionEquation import DiffusionEquation
from fipy.models.levelSet.distanceFunction.levelSetDiffusionVariable import LevelSetDiffusionVariable
from fipy.solvers.linearPCGSolver import LinearPCGSolver

class LevelSetDiffusionEquation(DiffusionEquation):
    
    def __init__(self,
                 var,
                 distanceVar = None,
                 transientCoeff = 1.,
                 diffusionCoeff = 1.,
                 solver = LinearPCGSolver(tolerance = 1.e-15,
                                          steps = 1000),
                 boundaryConditions = (),
                 otherTerms = ()):
        """
        
        A `LevelSetDiffusionEquation` is instantiated with the
        following arguments,

        `var` - The species concentration variable.

        `distanceVar` - A `DistanceVariable` object

        `transientCoeff` - In general 1 is used.

        `solver` - A given solver.

        `boundaryConditions` - A tuple of `BoundaryCondition` objects.

        """

	DiffusionEquation.__init__(
            self,
            var,
            transientCoeff = transientCoeff,
            diffusionCoeff = LevelSetDiffusionVariable(distanceVar,
                                                       diffusionCoeff),
            solver = solver,
            boundaryConditions = boundaryConditions,
            otherTerms = otherTerms
            )
        
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
