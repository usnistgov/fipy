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

Here we solve the level set equation in one dimension. The equation is
given by:

.. raw:: latex

    $$ | \\nabla \\phi | = 1 $$
    $$ \\phi = 0 \;\; \\text{at} \;\; x = L / 2 $$

Do the tests:

   >>> eqn.setInitialEvaluatedCells()
   >>> Numeric.allclose(Numeric.array(eqn.getVar()), Numeric.array((1., 1., 1., 1., dx / 2., -dx / 2., -1., -1., -1., -1.)), atol = 1e-10)
   1
   >>> Numeric.allclose(eqn.getInitialTrialCells(), Numeric.array((3, 6)))
   1
   >>> Numeric.allclose(Numeric.array(eqn.getVar()), Numeric.array((1., 1., 1., 3. * dx / 2., dx / 2., -dx / 2., -3. * dx / 2., -1., -1., -1.)))
   1
   >>> eqn.resetCells()
   >>> eqn.solve()
   >>> Numeric.allclose(Numeric.array(eqn.getVar()), Numeric.array((9. * dx / 2., 7. * dx / 2., 5. * dx / 2., 3. * dx / 2., dx / 2., -dx / 2., -3. * dx / 2., -5. * dx / 2., -7. * dx / 2., -9. * dx / 2.)))
   1

"""
__docformat__ = 'restructuredtext'

import Numeric
   
from fipy.meshes.grid2D import Grid2D
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.variables.cellVariable import CellVariable
from fipy.models.levelSet.distanceFunction.distanceFunctionEquation import DistanceFunctionEquation

dx = 0.5
dy = 2.
nx = 10
ny = 1

L = nx * dx

mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)

distanceFunctionVariable = CellVariable(
    name = 'level set variable',
    mesh = mesh,
    value = -1.
    )

distanceFunctionViewer = Grid2DGistViewer(var = distanceFunctionVariable, palette = 'rainbow.gp', minVal = -5., maxVal = 5.)

positiveCells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] < L / 2.)
distanceFunctionVariable.setValue(-1.)
distanceFunctionVariable.setValue(1.,positiveCells)

eqn = DistanceFunctionEquation(distanceFunctionVariable)

if __name__ == '__main__':
    distanceFunctionViewer.plot()

    eqn.solve()

    distanceFunctionViewer.plot()

    raw_input('finished run()')
