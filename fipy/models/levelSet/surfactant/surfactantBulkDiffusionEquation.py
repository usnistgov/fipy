#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "surfactantBulkDiffusionEquation.py"
 #                                    created: 8/31/04 {10:39:23 AM} 
 #                                last update: 4/1/05 {11:02:34 AM} 
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

The `SurfactantBulkDiffusionEquation` solves the bulk diffusion of a
species with a source term for the jump from the bulk to an interface.
The governing equation is given by,

.. raw:: latex

    $$ \\frac{\\partial c}{\\partial t} = \\nabla \\cdot D \\nabla  c $$

where,

.. raw:: latex

    $$ D = D_c \\;\\; \\text{when} \\;\\; \\phi > 0 $$
    $$ D = 0   \\;\\; \\text{when} \\;\\; \\phi \\le 0 $$

The jump condition at the interface is defined by Langmuir
adsorption. Langmuir adsorption essentially states that the ability for
a species to jump from an electrolyte to an interface is proportional to
the concentration in the electrolyte, available site density and a
jump coefficient. The boundary condition at the interface is given by

.. raw:: latex

    $$ D \\hat{n} \\cdot \\nabla c = -k c (1 - \\theta) \;\; \\text{at} \;\; \\phi = 0$$

"""

__docformat__ = 'restructuredtext'

import Numeric

from fipy.models.levelSet.distanceFunction.levelSetDiffusionEquation import buildLevelSetDiffusionEquation
from fipy.terms.implicitSourceTerm import ImplicitSourceTerm
from fipy.variables.cellVariable import CellVariable

class AdsorptionCoeff(CellVariable):
    def __init__(self, rateConstant = None, distanceVar = None):
        CellVariable.__init__(self, mesh = distanceVar.getMesh())
        self.distanceVar = self._requires(distanceVar)
        self.rateConstant = rateConstant

    def _calcValue(self):
        self.value = self.rateConstant * self.distanceVar.getCellInterfaceAreas() / self.mesh.getCellVolumes()

class ScAdsorptionCoeff(AdsorptionCoeff):
    def __init__(self, bulkVar = None, surfactantVar = None, rateConstant = None, distanceVar = None):
        AdsorptionCoeff.__init__(self, rateConstant = rateConstant, distanceVar = distanceVar)
        self.bulkVar = self._requires(bulkVar)
        self.surfactantVar = self._requires(surfactantVar)
    
    def _calcValue(self):
        AdsorptionCoeff._calcValue(self)
        bulk = Numeric.array(self.bulkVar)
        val = Numeric.array(self.value)
        self.value = val * bulk * Numeric.array(self.surfactantVar.getInterfaceVar())

def buildSurfactantBulkDiffusionEquation(bulkVar = None,
                                         distanceVar = None,
                                         surfactantVar = None,
                                         otherSurfactantVar = None,
                                         diffusionCoeff = None,
                                         transientCoeff = 1.,
                                         rateConstant = None):

        """
        
        A `SurfactantBulkDiffusionEquation` is instantiated with the
        following arguments,

        `var` - The bulk surfactant concentration variable.

        `distanceVar` - A `DistanceVariable` object

        `surfactantVariable` - A `SurfactantVariable` object

        `diffusionCoeff` - A float or a `FaceVariable`.

        `transientCoeff` - In general 1 is used.

        `jumpRate` - The adsorption coefficient.

        `boundaryConditions` - A tuple of `BoundaryCondition` objects.

        """
        
        spSourceTerm = ImplicitSourceTerm(AdsorptionCoeff(rateConstant = rateConstant,
                                                    distanceVar = distanceVar))

        coeff = ScAdsorptionCoeff(bulkVar = bulkVar,
                                  surfactantVar = surfactantVar,
                                  rateConstant = rateConstant,
                                  distanceVar = distanceVar)
                                  
        eq = buildLevelSetDiffusionEquation(ionVar = bulkVar,
                                            distanceVar = distanceVar,
                                            diffusionCoeff = diffusionCoeff,
                                            transientCoeff = transientCoeff)

        if otherSurfactantVar is not None:
            otherCoeff = ScAdsorptionCoeff(bulkVar = bulkVar,
                                           surfactantVar = otherSurfactantVar,
                                           rateConstant = rateConstant,
                                           distanceVar = distanceVar)
        else:
            otherCoeff = 0
            
        return eq - coeff + spSourceTerm - otherCoeff
