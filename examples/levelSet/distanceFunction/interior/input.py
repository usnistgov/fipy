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

Here we solve the level set equation in two dimension for an interior region. The equation is
given by:

.. raw:: latex

    $$ | \\nabla \\phi | = 1 $$
    $$ \\phi = 0 \;\; \\text{at} $$
    $$ x = \\left( d, L - d \\right) \;\; \\text{for} \;\; d \\le y \\le L - d $$
    $$ y = \\left( d, L - d \\right) \;\; \\text{for} \;\; d \\le x \\le L - d $$
    
Do the tests:

   >>> eqn.setInitialEvaluatedCells()
   >>> dX = dx / 2.
   >>> dY = dy / 2.
   >>> mm = min(dX, dY)
   >>> Numeric.allclose(Numeric.array(eqn.getVar()), Numeric.array((  1.  ,   dY  ,   dY  ,  dY  ,  1.  ,
   ...                                                         dX  ,   -mm ,   -dY ,  -mm ,  dX  ,
   ...                                                         dX  ,   -dX ,   -1. ,  -dX ,  dX  ,
   ...                                                         dX  ,   -mm ,   -dY ,  -mm ,  dX  ,
   ...                                                         1.  ,   dY  ,   dY  ,  dY  ,  1.  )), atol = 1e-10)
   1
   >>> Numeric.allclose(eqn.getInitialTrialCells(), Numeric.array((0, 4, 12, 20, 24)), atol = 1e-10)
   1
   >>> def evalCell(phix, phiy, dx, dy):
   ...     aa = dy**2 + dx**2
   ...     bb = -2 * ( phix * dy**2 + phiy * dx**2)
   ...     cc = dy**2 * phix**2 + dx**2 * phiy**2 - dx**2 * dy**2
   ...     sqr = Numeric.sqrt(bb**2 - 4. * aa * cc)
   ...     return ((-bb - sqr) / 2. / aa,  (-bb + sqr) / 2. / aa)
   >>> v1 = evalCell(dY, dX, dx, dy)[1]
   >>> v2 = max(-dY*3, -dX*3)   
   >>> _testInitialTrialValues = Numeric.array((  v1  ,   dY  ,   dY  ,  dY  ,  v1  ,
   ...                                                     dX  ,   -mm ,   -dY ,  -mm ,  dX  ,
   ...                                                     dX  ,   -dX ,   -v1 ,  -dX ,  dX  ,
   ...                                                     dX  ,   -mm ,   -dY ,  -mm ,  dX  ,
   ...                                                     v1  ,   dY  ,   dY  ,  dY  ,  v1  ))
   >>> Numeric.allclose(Numeric.array(eqn.getVar()), _testInitialTrialValues, atol = 1e-10)
   1
   >>> eqn.resetCells()
   >>> eqn.solve()
   >>> _testFinalValues = _testInitialTrialValues
   >>> Numeric.allclose(Numeric.array(eqn.getVar()), _testFinalValues, atol = 1e-10)
   1

"""

import Numeric

from fipy.meshes.grid2D import Grid2D
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.variables.cellVariable import CellVariable
from fipy.models.levelSet.distanceFunction.distanceFunctionEquation import DistanceFunctionEquation

dx = 1.
dy = 1.
nx = 5
ny = 5

Lx = nx * dx
Ly = ny * dy

mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)

distanceFunctionVariable = CellVariable(
    name = 'level set variable',
    mesh = mesh,
    value = -1.
    )

distanceFunctionViewer = Grid2DGistViewer(var = distanceFunctionVariable, palette = 'rainbow.gp', minVal = -5., maxVal = 5.)

positiveCells = mesh.getCells(filter = lambda cell: (cell.getCenter()[0] < dx) or (cell.getCenter()[0] > (Lx - dx)) or (cell.getCenter()[1] < dy) or (cell.getCenter()[1] > (Ly - dy)))

distanceFunctionVariable.setValue(-1.)
distanceFunctionVariable.setValue(1.,positiveCells)

eqn = DistanceFunctionEquation(distanceFunctionVariable)

if __name__ == '__main__':
    distanceFunctionViewer.plot()
    eqn.solve()
    distanceFunctionViewer.plot()
    raw_input('finished')
