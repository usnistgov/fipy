#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "distanceVariable.py"
 #                                    created: 7/29/04 {10:39:23 AM} 
 #                                last update: 10/19/04 {4:38:29 PM} 
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

r"""

.. raw:: latex

    A `DistanceVariable` object calculates $\phi$ so it satisfies,

    $$ | \\nabla \\phi | = 1 $$

    using the fast marching method with an initial condition defined by
    the zero level set.

Currently the solution is first order, This suffices for initial
conditions with straight edges (e.g. trenches in
electrodeposition). The method should work for unstructured 2D grids
but testing on unstructured grids is untested thus far. This is a 2D
implementation as it stands. Extending to 3D should be relatively
simple.

Here we will define a few test cases. Firstly a 1D test case

   >>> from fipy.meshes.grid2D import Grid2D
   >>> mesh = Grid2D(dx = .5, dy = .2, nx = 8, ny = 1)
   >>> from distanceVariable import DistanceVariable
   >>> var = DistanceVariable(mesh = mesh, value = (-1, -1, -1, -1, 1, 1, 1, 1))
   >>> answer = (-1.75, -1.25, -.75, -0.25, 0.25, 0.75, 1.25, 1.75) 
   >>> Numeric.allclose(answer, var)
   1

A 1D test case with very small dimensions.

   >>> dx = 1e-10
   >>> mesh = Grid2D(dx = dx, dy = 1., nx = 8, ny = 1)
   >>> var = DistanceVariable(mesh = mesh, value = (-1, -1, -1, -1, 1, 1, 1, 1))
   >>> answer = Numeric.arange(8) * dx - 3.5 * dx
   >>> Numeric.allclose(answer, var)
   1


For future reference, the minimum distance for the interface cells can
be calculated with the following functions. The trial cell values will
also be calculated with these functions. In essence it is not
difficult to calculate the level set distance function on an
unstructured 3D grid. However a lot of testing will be required. The
minimum distance functions will take the following form.

.. raw:: latex

    $$ X_{\text{min}} = \frac{\left| \vec{s} \cross \vec{t} \right|}
    {\left| \vec{s} - \vec{t} \right|} $$

    and in 3D,

    $$ X_{\text{min}} = \frac{1}{3!} \left| \vec{s} \cdot \left(
    \vec{t} \cross \vec{u} \right) \left| $$

    where the vectors $\vec{s}$, $\vec{t}$ and $\vec{u}$ represent the
    vectors from the cell of interest to the neighboring cell.
    
"""
__docformat__ = 'restructuredtext'

import Numeric
import MA

from fipy.meshes.numMesh.mesh import MAtake
from fipy.variables.cellVariable import CellVariable

class DistanceVariable(CellVariable):
    def __init__(self, mesh, name = '', value = 0., unit = None, hasOld = 1, narrowBandWidth = 1e+10):
        CellVariable.__init__(self, mesh, name = name, value = value, unit = unit, hasOld = hasOld)
        self.markStale()
        self.narrowBandWidth = narrowBandWidth
        
    def _calcValue(self):

        ## calculate interface values

        cellToCellIDs = self.mesh.getCellToCellIDs()
        adjVals = MAtake(self.value, cellToCellIDs)
        adjInterfaceValues = MA.masked_array(adjVals, mask = (adjVals * self.value[:,Numeric.NewAxis]) > 0)
        dAP = self.mesh.getCellToCellDistances()
        distances = MA.sort(abs(self.value[:,Numeric.NewAxis] * dAP / (self.value[:,Numeric.NewAxis] - adjInterfaceValues)), 1)
        sign = (self.value > 0) * 2 - 1
        s = distances[:,0]
        t = distances[:,1]
        signedDistance = MA.where(s.mask(),
                                  self.value,
                                  MA.where(t.mask(),
                                           sign * s,
                                           sign * s * t / MA.sqrt(s**2 + t**2)))

        self.value = signedDistance

        ## calculate initial trial values
        interfaceFlag = Numeric.sum(Numeric.logical_not(distances.mask()), 1) > 0
        adjInterfaceFlag = MAtake(interfaceFlag, cellToCellIDs)
        hasAdjInterface = Numeric.sum(adjInterfaceFlag.filled(), 1) > 0
        trialFlag = Numeric.logical_and(Numeric.logical_not(interfaceFlag), hasAdjInterface) 
        trialIDs = list(Numeric.nonzero(trialFlag))
        evaluatedFlag = interfaceFlag
        
        for id in trialIDs:
            self.value[id] = self._calcTrialValue(id, evaluatedFlag)
            
        while len(trialIDs) > 0:

            id = trialIDs[Numeric.argmin(abs(Numeric.take(self.value, trialIDs)))]
            trialIDs.remove(id)
            evaluatedFlag[id] = 1
            if abs(self.value[id]) > self.narrowBandWidth / 2:
                break

            for adjID in cellToCellIDs[id].filled(fill_value = -1):
                if adjID != -1:
                    if not evaluatedFlag[adjID]:
                        self.value[adjID] = self._calcTrialValue(adjID, evaluatedFlag)
                        if adjID not in trialIDs:
                            trialIDs.append(adjID)

        self.value = Numeric.array(self.value)

    def _calcTrialValue(self, id, evaluatedFlag):
        adjIDs = self.mesh.getCellToCellIDs()[id]
        adjEvaluatedFlag = MAtake(evaluatedFlag, adjIDs)
        adjValues = MA.masked_array(MAtake(self.value, adjIDs), MA.logical_not(adjEvaluatedFlag))
        indices = MA.argsort(abs(adjValues))
        sign = (self.value[id] > 0) * 2 - 1
        d0 = self.mesh.getCellToCellDistances()[id, indices[0]]
        v0 = self.value[adjIDs[indices[0]]]

        N = Numeric.sum(Numeric.logical_not(adjValues.mask()))

        if N == 0:
            raise Error 
        elif N == 1:
            return v0 + sign * d0
        else:
            d1 = self.mesh.getCellToCellDistances()[id, indices[1]]
            n0 = self.mesh.getCellNormals()[id, indices[0]]
            n1 = self.mesh.getCellNormals()[id, indices[1]]
            v1 = self.value[adjIDs[indices[1]]]

            dotProd = d0 * d1 * Numeric.dot(n0, n1)
            crossProd = d0 * d1 * (n0[0] * n1[1] - n0[1] * n1[0])
            dsq = d0**2 + d1**2 - 2 * dotProd
            
            top = -v0 * (dotProd - d1**2) - v1 * (dotProd - d0**2)
            sqrt = crossProd**2 *(dsq - (v0 - v1)**2)
            sqrt = Numeric.sqrt(max(sqrt, 0))
            
            return (top + sign * sqrt) / dsq
        
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test()         
    



            
            
        
                
