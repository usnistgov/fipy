#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "surfactantVariable.py"
 #                                    created: 7/29/04 {10:39:23 AM} 
 #                                last update: 10/19/04 {4:41:41 PM} 
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

from fipy.variables.cellVariable import CellVariable
import Numeric

class SurfactantVariable(CellVariable):
    """

    `SurfactantVariable` initializes itself across the zero level set
    interface using the delta function. The `value` argument
    corresponds to the initial concentration of surfactant on the
    interface (moles divided by area). The value held by the
    `SurfactantVariable` is actually a volume density (moles divided by
    volume).

    """
    
    def __init__(self, value = 0., distanceVar = None, name = 'surfactant variable'):
        """

        A simple 1D test:

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(dx = 1., dy = 1., nx = 4, ny = 1)
           >>> from fipy.models.levelSet.distanceFunction.distanceVariable \\
           ...     import DistanceVariable
           >>> distanceVariable = DistanceVariable(mesh = mesh, 
           ...                                     value = (-1.5, -0.5, 0.5, 941.5))
           >>> surfactantVariable = SurfactantVariable(value = 1, 
           ...                                         distanceVar = distanceVariable)
           >>> Numeric.allclose(surfactantVariable, (0, 0., 1., 0))
           1

        A 2D test case:

           >>> mesh = Grid2D(dx = 1., dy = 1., nx = 3, ny = 3)
           >>> distanceVariable = DistanceVariable(mesh = mesh,
           ...                                     value = (1.5, 0.5, 1.5,
           ...                                              0.5,-0.5, 0.5,
           ...                                              1.5, 0.5, 1.5))
           >>> surfactantVariable = SurfactantVariable(value = 1, 
           ...                                         distanceVar = distanceVariable)
           >>> Numeric.allclose(surfactantVariable, (0, 1, 0, 1, 0, 1, 0, 1, 0))
           1

        Another 2D test case:

           >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
           >>> distanceVariable = DistanceVariable(mesh = mesh, 
           ...                                     value = (-0.5, 0.5, 0.5, 1.5))
           >>> surfactantVariable = SurfactantVariable(value = 1, 
           ...                                         distanceVar = distanceVariable)
           >>> Numeric.allclose(surfactantVariable, 
           ...                  (0, Numeric.sqrt(2), Numeric.sqrt(2), 0))
           1
           
        """
        
        CellVariable.__init__(self, mesh = distanceVar.getMesh(), name = name)

        self.value = distanceVar.getCellInterfaceAreas() * value / self.mesh.getCellVolumes()

        self.distanceVar = distanceVar

    def getInterfaceValue(self):
        areas = self.distanceVar.getCellInterfaceAreas()        
        areas = Numeric.where(areas > 1e-20, areas, 1)
        value = Numeric.array(self) * self.mesh.getCellVolumes() / areas
        return value

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
