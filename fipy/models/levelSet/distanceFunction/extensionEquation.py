#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "distanceFunctionEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 6/3/04 {2:53:54 PM} 
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

A `ExtensionEquation` object solves the equation,

.. raw:: latex

    $$ \\nabla u \\cdot \\nabla \\phi = 0 $$

using the fast marching method with an initial condition defined at
the zero level set. Essentially the equation solves a fake distance
function equation to march out the velocity from the interface.

   >>> from fipy.meshes.grid2D import Grid2D
   >>> from distanceVariable import DistanceVariable
   >>> from fipy.variables.cellVariable import CellVariable
   >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 2)
   >>> var = DistanceVariable(mesh = mesh, value = (-1, 1, 1, 1))
   >>> extensionVar = CellVariable(mesh = mesh, value = (-1, .5, 2, -1))
   >>> setValueFlag = ExtensionEquation(var, extensionVar).solve()
   >>> Numeric.allclose((-1, 1, 1, 1), Numeric.array(var))
   1
   >>> Numeric.allclose((1.25, .5, 2, 1.25), Numeric.array(extensionVar))
   1

   >>> mesh = Grid2D(dx = 1., dy = 1., nx = 3, ny = 3)
   >>> var = DistanceVariable(mesh = mesh, value = (-1, 1, 1,
   ...                                               1, 1, 1,
   ...                                               1, 1, 1))
   >>> extensionVar = CellVariable(mesh = mesh, value = (-1, .5, -1,
   ...                                                    2, -1, -1,
   ...                                                   -1, -1, -1))
   >>> setValueFlag = ExtensionEquation(var, extensionVar).solve()
   >>> answer = (1.25, .5, .5, 2, 1.25, 0.9544, 2, 1.5456, 1.25)
   >>> Numeric.allclose(answer, Numeric.array(extensionVar), atol = 1e-5)
   1

"""
__docformat__ = 'restructuredtext'

import Numeric

from distanceEquation import DistanceEquation

class ExtensionEquation(DistanceEquation):

    def __init__(self, var = None, extensionVar = None):
        """

        The `var` argument must contain both positive and negative
        values to define the zero level set.

        The `extensionVar` must be defined on the positive interface
        cells.
        
        """
        self.extensionVar = extensionVar
        DistanceEquation.__init__(self, var.copy())

    def _calcInterfaceValues(self):
        """

        Sets the values in cells at the interface (cells that have a neighbour
        of the opposite sign) to the shortest signed distance.

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 1)
           >>> from distanceVariable import DistanceVariable
           >>> from fipy.variables.cellVariable import CellVariable
           >>> var = DistanceVariable(mesh = mesh, value = (-1, 1))
           >>> extensionVar = CellVariable(mesh = mesh, value = (0, 1))
           >>> ExtensionEquation(var, extensionVar)._calcInterfaceValues()
           [1,1,]
           >>> print extensionVar
           [ 1., 1.,]

           >>> mesh = Grid2D(dx = 1., dy = 1., nx = 3, ny = 3)
           >>> var = DistanceVariable(mesh = mesh, value =
           ...     (-1, -1, -1, -1, 1, -1, -1, -1, -1))
           >>> extensionVar = CellVariable(mesh = mesh, value =
           ...     (-1, -1, -1, -1, .5, -1, -1, -1, -1))
           >>> ExtensionEquation(var, extensionVar)._calcInterfaceValues()
           [0,1,0,1,1,1,0,1,0,]
           >>> answer = Numeric.array((-1, .5, -1, .5, .5, .5,-1, .5, -1))
           >>> Numeric.allclose(answer, Numeric.array(extensionVar))
           1

           >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 2)
           >>> var = DistanceVariable(mesh = mesh, value = (-1, 1, 1, 1))
           >>> extensionVar = CellVariable(mesh = mesh, value = (-1, .5, 2, -1))
           >>> setValueFlag = ExtensionEquation(var, extensionVar)._calcInterfaceValues()
           >>> Numeric.allclose((1.25, .5, 2, -1), Numeric.array(extensionVar))
           1

        """
        
        setValueFlag = DistanceEquation._calcInterfaceValues(self)
        
        positiveInterfaceCellFlag = Numeric.logical_and(Numeric.where(setValueFlag == 1, 1, 0), self.var > 0)
        negativeInterfaceCellIDs = Numeric.nonzero(Numeric.logical_and(Numeric.where(setValueFlag == 1, 1, 0), self.var < 0))

        
        
        for cellID in negativeInterfaceCellIDs:
            DistanceEquation._calcTrialValue(self, cellID, positiveInterfaceCellFlag)
                        
        return setValueFlag


    def _calcLinear(self, phi1, d1, cellID, adjCellID):
        self.extensionVar[cellID] = self.extensionVar[adjCellID]
        return DistanceEquation._calcLinear(self, phi1, d1, cellID, adjCellID)

    def _calcQuadratic(self, phi1, phi2, n1, n2, d1, d2, area1, area2, cellID, adjCellID1, adjCellID2):
        val = DistanceEquation._calcQuadratic(self, phi1, phi2, n1, n2, d1, d2, area1, area2, cellID, adjCellID1, adjCellID2)

        if self.var[cellID] > 0:
            phi = val[0]
        else:
            phi = val[1]

        n1grad = (phi1 - phi) / d1
        n2grad = (phi2 - phi) / d2

        u1 = self.extensionVar[adjCellID1]
        u2 = self.extensionVar[adjCellID2]

        self.extensionVar[cellID] = (u1 * n1grad * area1 + u2 * n2grad * area2) / (area1 * n1grad + area2 * n2grad)

        return val

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 




