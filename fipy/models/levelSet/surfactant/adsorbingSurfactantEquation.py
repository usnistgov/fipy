#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "adsorbingSurfactantEquation.py"
 #                                    created: 8/31/04 {10:39:23 AM} 
 #                                last update: 8/31/04 {4:00:26 PM} 
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

The `AdsorbingSurfactantEquation` object solves the
`SurfactantEquation but with an adsorbing species from some bulk
value. The equation that describes the surfactant adsorbing is given
by,

.. raw:: latex

    $$ \\dot{\\theta} = J v \\theta + k c (1 - \\theta) $$

This last term in this equation accounts for Langmuir type adsorption
from the bulk. It assumes a vacant proportion of surface sites. The adsorption term
is added to the source in teh follwoing way,

.. raw:: latex

    $$ S_c = k c \;\; \\text{and} \;\; S_p = k c $$

The following is a test case:

   >>> from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
   >>> from fipy.models.levelSet.surfactant.surfactantVariable import SurfactantVariable
   >>> from fipy.meshes.grid2D import Grid2D
   >>> dx = .5
   >>> dy = 2.3
   >>> dt = 0.25
   >>> k = 0.56
   >>> initialValue = 0.1
   >>> c = 0.2
   >>> mesh = Grid2D(dx = dx, dy = dy, nx = 5, ny = 1)
   >>> distanceVar = DistanceVariable(mesh = mesh, value = (-dx*3/2, -dx/2, dx/2, 3*dx/2 ,5*dx/2))
   >>> var = SurfactantVariable(value = (0, 0, initialValue, 0 ,0), distanceVar = distanceVar)
   >>> bulkVar = CellVariable(mesh = mesh, value = (c , c, c, c, c))
   >>> eqn = AdsorbingSurfactantEquation(var, distanceVar, bulkVar, k)
   >>> eqn.solve(dt = dt)
   >>> answer = (initialValue + dt * k * c) / (1 + dt * k * c)
   >>> Numeric.allclose(var.getInterfaceValue(), Numeric.array((0, 0, answer, 0, 0)))
   1

"""

import Numeric

from fipy.variables.cellVariable import CellVariable
from surfactantEquation import SurfactantEquation
from fipy.terms.spSourceTerm import SpSourceTerm
from fipy.terms.scSourceTerm import ScSourceTerm

 
class AdsorptionCoeff(CellVariable):
    def __init__(self, distanceVar, bulkVar, rateConstant):
        CellVariable.__init__(self, mesh = distanceVar.getMesh())

        self.distanceVar = self.requires(distanceVar)
        self.bulkVar = self.requires(bulkVar)
        self.rateConstant = rateConstant
        self.dt = 0

    def _calcValue(self):
        self.value = self.dt * Numeric.array(self.bulkVar) * self.rateConstant * self.multiplier()

    def updateDt(self, dt):
        self.dt = dt
        self.markStale()

class SpAdsorptionCoeff(AdsorptionCoeff):
    def multiplier(self):
        return self.distanceVar.getCellInterfaceFlag()
    
class ScAdsorptionCoeff(AdsorptionCoeff):
    def multiplier(self):
        return self.distanceVar.getCellInterfaceAreas() / self.mesh.getCellVolumes() 
 
class AdsorbingSurfactantEquation(SurfactantEquation):
    def __init__(self,
                 var,
                 distanceVar,
                 bulkVar,
                 rateConstant):
        
        SurfactantEquation.__init__(self, var, distanceVar)
        
        self.spCoeff = SpAdsorptionCoeff(distanceVar, bulkVar, rateConstant)
        self.scCoeff = ScAdsorptionCoeff(distanceVar, bulkVar, rateConstant)

        self.terms += (
            SpSourceTerm(self.spCoeff, self.var.getMesh()),
            ScSourceTerm(self.scCoeff, self.var.getMesh())
            )

    def solve(self, dt):
        self.dt = dt
        self.scCoeff.updateDt(dt)
        self.spCoeff.updateDt(dt)
        SurfactantEquation.solve(self, dt)
        
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
