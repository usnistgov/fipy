#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "metalIonSourceVariable.py"
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

from fipy.tools import numerix

from fipy.variables.cellVariable import CellVariable

class _MetalIonSourceVariable(CellVariable):
    r"""
    The `_MetalIonSourceVariable` object evaluates the source
    coefficient for the `MetalIonEquation`. The source only exists
    on the cells in front of the interface and simulates the loss
    of material from bulk diffusion as the interface advances. The
    source is given by,

    .. math::

       D \hat{n} \cdot \nabla c = \frac{v(c)}{\Omega} \qquad \text{at $\phi = 0$}

    Here is a test,

    >>> from fipy.meshes import Grid2D
    >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 2)
    >>> from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
    >>> distance = DistanceVariable(mesh = mesh, value = (-.5, .5, .5, 1.5))
    >>> ionVar = CellVariable(mesh = mesh, value = (1, 1, 1, 1))
    >>> depositionRate = CellVariable(mesh=mesh, value=(1, 1, 1, 1))
    >>> arr = _MetalIonSourceVariable(ionVar, distance, depositionRate, 1)
    >>> sqrt = numerix.sqrt(2)
    >>> ans = CellVariable(mesh=mesh, value=(0, 1 / sqrt, 1 / sqrt, 0))
    >>> print numerix.allclose(arr, ans)
    True
    """
    def __init__(self, ionVar = None, distanceVar = None, depositionRate = None, metalIonMolarVolume = None):
        """
        Creates a `_MetalIonSourceVariable` object.

        :Parameters:
          - `ionVar` : The metal ion concentration.
          - `distanceVar` : A `DistanceVariable` object.
          - `depositionRate` : The deposition rate.
          - `metalIonMolarVolume` : Molar volume of the metal ions.
       
        """
        
        CellVariable.__init__(self, distanceVar.mesh, hasOld = 0)
        self.ionVar = self._requires(ionVar)
        self.distanceVar = self._requires(distanceVar)
        self.depositionRate = self._requires(depositionRate)
        self.metalIonMolarVolume = metalIonMolarVolume
        
    def _calcValue(self):
        ionVar = numerix.array(self.ionVar)
        ionVar = numerix.where(ionVar > 1e-20, ionVar, 1e-20)
        return numerix.array(self.depositionRate) * self.distanceVar.cellInterfaceAreas / (self.mesh.cellVolumes * self.metalIonMolarVolume) / ionVar
    
        
def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test() 
