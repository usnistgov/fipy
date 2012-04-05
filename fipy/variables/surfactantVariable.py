#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "surfactantVariable.py"
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

from fipy.variables.cellVariable import CellVariable
from fipy.tools import numerix
from fipy.tools.decorators import getsetDeprecated

__all__ = ["SurfactantVariable"]

class SurfactantVariable(CellVariable):
    """

    The `SurfactantVariable` maintains a conserved volumetric
    concentration on cells adjacent to, but in front of, the
    interface. The `value` argument corresponds to the initial
    concentration of surfactant on the interface (moles divided by
    area). The value held by the `SurfactantVariable` is actually a
    volume density (moles divided by volume).

    """
    
    def __init__(self, value = 0., distanceVar = None, name = 'surfactant variable', hasOld=False):
        """

        A simple 1D test:

        >>> from fipy.meshes import Grid1D
        >>> mesh = Grid1D(dx = 1., nx = 4)
        >>> from fipy.variables.distanceVariable import DistanceVariable
        >>> distanceVariable = DistanceVariable(mesh = mesh, 
        ...                                     value = (-1.5, -0.5, 0.5, 941.5))
        >>> surfactantVariable = SurfactantVariable(value = 1, 
        ...                                         distanceVar = distanceVariable)
        >>> print numerix.allclose(surfactantVariable, (0, 0., 1., 0))
        1

        A 2D test case:

        >>> from fipy.meshes import Grid2D
        >>> mesh = Grid2D(dx = 1., dy = 1., nx = 3, ny = 3)
        >>> distanceVariable = DistanceVariable(mesh = mesh,
        ...                                     value = (1.5, 0.5, 1.5,
        ...                                              0.5,-0.5, 0.5,
        ...                                              1.5, 0.5, 1.5))
        >>> surfactantVariable = SurfactantVariable(value = 1, 
        ...                                         distanceVar = distanceVariable)
        >>> print numerix.allclose(surfactantVariable, (0, 1, 0, 1, 0, 1, 0, 1, 0))
        1

        Another 2D test case:

        >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
        >>> distanceVariable = DistanceVariable(mesh = mesh, 
        ...                                     value = (-0.5, 0.5, 0.5, 1.5))
        >>> surfactantVariable = SurfactantVariable(value = 1, 
        ...                                         distanceVar = distanceVariable)
        >>> print numerix.allclose(surfactantVariable, 
        ...                  (0, numerix.sqrt(2), numerix.sqrt(2), 0))
        1

        :Parameters:
          - `value`: The initial value.
          - `distanceVar`: A `DistanceVariable` object.
          - `name`: The name of the variable.
          
        """


        CellVariable.__init__(self, mesh = distanceVar.mesh, name = name, hasOld=False)

        self.distanceVar = self._requires(distanceVar)
        self._value = numerix.array(distanceVar.cellInterfaceAreas) * value / self.mesh.cellVolumes

        if hasOld:
            self._old = self.copy()
        else:
            self._old = None

        self.interfaceSurfactantVariable = None

    @getsetDeprecated
    def getInterfaceVar(self):
        return self.interfaceVar

    @property
    def interfaceVar(self):
        """
        
        Returns the `SurfactantVariable` rendered as an
        `_InterfaceSurfactantVariable` which evaluates the surfactant
        concentration as an area concentration the interface rather
        than a volumetric concentration.

        """
        if self.interfaceSurfactantVariable is None:
            self.interfaceSurfactantVariable = _InterfaceSurfactantVariable(self)

        return self.interfaceSurfactantVariable

    @getsetDeprecated(new_name="distanceVar")
    def _getDistanceVar(self):
        return self.distanceVar
   
    def _calcValue(self):
        return self._value

    def copy(self):
        return self.__class__(
            distanceVar=self.distanceVar,
            name=self.name + "_old", 
            value=self.value.copy(),
            hasOld=False)
    
class _InterfaceSurfactantVariable(CellVariable):
    def __init__(self, surfactantVar):
        CellVariable.__init__(self, name = surfactantVar.name + "_interface", mesh = surfactantVar.mesh)
        self.surfactantVar = self._requires(surfactantVar)

    def _calcValue(self):
        areas = numerix.array(self.surfactantVar.distanceVar.cellInterfaceAreas)
        areas = numerix.where(areas > 1e-20, areas, 1)
        return numerix.array(self.surfactantVar) * self.mesh.cellVolumes / areas

def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test() 
