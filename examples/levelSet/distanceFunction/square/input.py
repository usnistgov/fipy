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

Here we solve the level set equation in two dimensions for a square. The equation is
given by:

.. raw:: latex

    $$ | \\nabla \\phi | = 1 $$
    $$ \\phi = 0 \;\; \\text{at} $$
    $$ x = \\left( L / 3, 2 L / 3 \\right) \;\; \\text{for} \;\; L / 3 \\le y \\le 2 L / 3 $$
    $$ y = \\left( L / 3, 2 L / 3 \\right) \;\; \\text{for} \;\; L / 3 \\le x \\le 2 L / 3 $$
    

Do the tests:

   >>> eqn.setInitialEvaluatedCells()
   >>> Numeric.allclose(Numeric.array(eqn.getVar()), Numeric.array((-1., -1., -1., -1., -1.,
   ...                                      -1. ,-1., -dy / 2., -1., -1.,
   ...                                      -1., -dx / 2., min(dx, dy) / 2., -dx / 2., -1.,
   ...                                      -1., -1., -dy / 2., -1., -1.,
   ...                                      -1., -1., -1., -1., -1.)), atol = 1e-10)
   1
   >>> Numeric.allclose(Numeric.sort(eqn.getInitialTrialCells()), Numeric.array((2, 6, 8, 10, 14, 16, 18, 22)), atol = 1e-10)
   1
   >>> def evalCell(phix, phiy, dx, dy):
   ...     aa = dy**2 + dx**2
   ...     bb = -2 * ( phix * dy**2 + phiy * dx**2)
   ...     cc = dy**2 * phix**2 + dx**2 * phiy**2 - dx**2 * dy**2
   ...     sqr = Numeric.sqrt(bb**2 - 4. * aa * cc)
   ...     return ((-bb - sqr) / 2. / aa,  (-bb + sqr) / 2. / aa)
   >>> val = evalCell(-dy / 2., -dx / 2., dx, dy)[0]
   >>> Numeric.allclose(Numeric.array(eqn.getVar()), Numeric.array((-1.          , -1.     , -3. * dy / 2.   , -1.     , -1.,
   ...                                           -1.          , val     , -dy / 2.        , val     , -1.,
   ...                                           -3. * dx / 2., -dx / 2., min(dx, dy) / 2., -dx / 2., -3. * dx / 2.,
   ...                                           -1.          , val     , -dy / 2.        , val     , -1.,
   ...                                           -1.          , -1.     , -3. * dy / 2.   , -1.     , -1.)), atol = 1e-10)
   1
   >>> eqn.resetCells()
   >>> eqn.solve()
   >>> v1 = evalCell(val, -3. * dx / 2., dx, dy)[0]
   >>> v2 = evalCell(-3. * dy / 2., val, dx, dy)[0]
   >>> v3 = evalCell(v2, v1, dx, dy)[0]
   >>> Numeric.allclose(Numeric.array(eqn.getVar()), Numeric.array((v3           , v2      , -3. * dy / 2.   , v2      , v3,
   ...                        v1           , val     , -dy / 2.        , val     , v1           ,
   ...                        -3. * dx / 2., -dx / 2., min(dx, dy) / 2., -dx / 2., -3. * dx / 2.,
   ...                        v1           , val     , -dy / 2.        , val     , v1           ,
   ...                        v3           , v2      , -3. * dy / 2.   , v2      , v3           )), atol = 1e-10)
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

positiveCells = mesh.getCells(filter = lambda cell: Lx / 3. < cell.getCenter()[0] < 2. * Lx / 3. and Ly / 3. < cell.getCenter()[1] < 2. * Ly / 3)
distanceFunctionVariable.setValue(-1.)
distanceFunctionVariable.setValue(1.,positiveCells)

eqn = DistanceFunctionEquation(distanceFunctionVariable)

if __name__ == '__main__':
    distanceFunctionViewer = Grid2DGistViewer(var = distanceFunctionVariable, palette = 'rainbow.gp', minVal = -5., maxVal = 5.)
    distanceFunctionViewer.plot()
    eqn.solve()
    distanceFunctionViewer.plot()
    raw_input('finished')
