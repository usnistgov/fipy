#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 4/2/04 {4:01:07 PM} { 1:23:41 PM}
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
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

"""

Here we solve the level set equation in two dimensions for a circle. The equation is
given by:

.. raw:: latex

    $$ | \\nabla \\phi | = 1 $$
    $$ \\phi = 0 \;\; \\text{at} \;\;  (x - L / 2)^2 + (y - L / 2)^2 = (L / 4)^2 $$

Do the tests:

   >>> eqn.setInitialEvaluatedCells()
   >>> dY = dy / 2.
   >>> dX = dx / 2.
   >>> mm = min (dX, dY)     
   >>> _testInitialEvaluatedValues = Numeric.array((-1.,  -1.  , -1.  ,  -1.  ,  -1.  ,  -1.  , -1.  ,  -1.  ,  -1.  , -1., -1.,
   ...                                                   -1.,  -1.  , -1.  ,  -1.  ,  -1.  ,  -1.  , -1.  ,  -1.  ,  -1.  , -1., -1.,
   ...                                                   -1.,  -1.  , -1.  ,  -1.  ,  -dY  ,  -dY  , -dY  ,  -1.  ,  -1.  , -1., -1.,
   ...                                                   -1.,  -1.  , -1.  ,  -mm  ,  mm   ,  dY   , mm   ,  -mm  ,  -1.  , -1., -1.,
   ...                                                   -1.,  -1.  , -dX  ,  mm   ,  1.   ,  1    , 1.   ,  mm   ,  -dX  , -1., -1.,
   ...                                                   -1.,  -1.  , -dX  ,  dX   ,  1.   ,  1    , 1.   ,  dX   ,  -dX  , -1., -1.,
   ...                                                   -1.,  -1.  , -dX  ,  mm   ,  1.   ,  1.   , 1.   ,  mm   ,  -dX  , -1., -1.,
   ...                                                   -1.,  -1.  , -1.  ,  -mm  ,  mm   ,  dY   , mm   ,  -mm  ,  -1.  , -1., -1.,
   ...                                                   -1.,  -1.  , -1.  ,  -1.  ,  -dY  ,  -dY  , -dY  ,  -1.  ,  -1.  , -1., -1.,
   ...                                                   -1.,  -1.  , -1.  ,  -1.  ,  -1.  ,  -1.  , -1.  ,  -1.  ,  -1.  , -1., -1.,
   ...                                                   -1.,  -1.  , -1.  ,  -1.  ,  -1.  ,  -1.  , -1.  ,  -1.  ,  -1.  , -1., -1.))
   >>> Numeric.allclose(Numeric.array(eqn.getVar()), _testInitialEvaluatedValues, atol = 1e-10)
   1

   >>> _testInitialTrialCellIDs = Numeric.array((15, 16, 17,
   ...                                               25, 29,
   ...                                               35, 41,
   ...                                               45, 48, 49, 50, 53,
   ...                                               56, 59, 61, 64,
   ...                                               67, 70, 71, 72, 75,
   ...                                                79, 85,
   ...                                               91, 95,
   ...                                               103, 104, 105))
   >>> Numeric.allclose(Numeric.sort(eqn.getInitialTrialCells()), _testInitialTrialCellIDs)
   1

   >>> def evalCell(phix, phiy, dx, dy):
   ...     aa = dy**2 + dx**2
   ...     bb = -2 * ( phix * dy**2 + phiy * dx**2)
   ...     cc = dy**2 * phix**2 + dx**2 * phiy**2 - dx**2 * dy**2
   ...     sqr = Numeric.sqrt(bb**2 - 4. * aa * cc)
   ...     return ((-bb - sqr) / 2. / aa,  (-bb + sqr) / 2. / aa)
   >>> v1 = evalCell(-dY, -mm, dx, dy)[0] 
   >>> v2 = evalCell(-mm, -dX, dx, dy)[0]
   >>> v3 = evalCell(mm,  mm,  dx, dy)[1]
   >>> _testInitialTrialValues = Numeric.array((-1.,  -1.  , -1.  ,  -1.  ,  -1.  ,  -1.  , -1.  ,  -1.  ,  -1.  , -1.  , -1.,
   ...                                                   -1.,  -1.  , -1.  ,  -1.  ,  -3*dY,  -3*dY, -3*dY,  -1.  ,  -1.  , -1.  , -1.,
   ...                                                   -1.,  -1.  , -1.  ,  v1   ,  -dY  ,  -dY  , -dY  ,  v1   ,  -1.  , -1.  , -1.,
   ...                                                   -1.,  -1.  , v2   ,  -mm  ,  mm   ,  dY   , mm   ,  -mm  ,  v2   , -1.  , -1.,
   ...                                                   -1.,  -dX*3, -dX  ,  mm   ,  v3   ,  3*dY , v3   ,  mm   ,  -dX  , -dX*3, -1.,
   ...                                                   -1.,  -dX*3, -dX  ,  dX   ,  3*dX ,  1    , 3*dX ,  dX   ,  -dX  , -dX*3, -1.,
   ...                                                   -1.,  -dX*3, -dX  ,  mm   ,  v3   ,  3*dY , v3   ,  mm   ,  -dX  , -dX*3, -1.,
   ...                                                   -1.,  -1.  , v2   ,  -mm  ,  mm   ,  dY   , mm   ,  -mm  ,  v2   , -1.  , -1.,
   ...                                                   -1.,  -1.  , -1.  ,  v1   ,  -dY  ,  -dY  , -dY  ,  v1   ,  -1.  , -1.  , -1.,
   ...                                                   -1.,  -1.  , -1.  ,  -1.  ,  -3*dY,  -3*dY, -3*dY,  -1.  ,  -1.  , -1.  , -1.,
   ...                                                   -1.,  -1.  , -1.  ,  -1.  ,  -1.  ,  -1.  , -1.  ,  -1.  ,  -1.  , -1.  , -1.))
   >>> Numeric.allclose(Numeric.array(eqn.getVar()), _testInitialTrialValues)
   1
   
"""
__docformat__ = 'restructuredtext'

import Numeric

from fipy.meshes.grid2D import Grid2D
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.variables.cellVariable import CellVariable
from fipy.models.levelSet.distanceFunction.distanceFunctionEquation import DistanceFunctionEquation

dx = 1.
dy = 1.
nx = 11
ny = 11
Lx = nx * dx
Ly = ny * dy

mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)

distanceFunctionVariable = CellVariable(
    name = 'level set variable',
    mesh = mesh,
    value = -1.
    )

distanceFunctionViewer = Grid2DGistViewer(var = distanceFunctionVariable, palette = 'rainbow.gp', minVal = -5., maxVal = 5.)

positiveCells = mesh.getCells(filter = lambda cell: (cell.getCenter()[0] - Lx / 2.)**2 + (cell.getCenter()[1] - Ly / 2.)**2 < (Lx / 4.)**2)
distanceFunctionVariable.setValue(-1.)
distanceFunctionVariable.setValue(1.,positiveCells)

eqn = DistanceFunctionEquation(distanceFunctionVariable)

if __name__ == '__main__':

    distanceFunctionViewer.plot()

    eqn.solve()

    distanceFunctionViewer.plot()

    raw_input('finished')
