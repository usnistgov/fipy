#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 3/4/05 {7:03:02 PM} { 1:23:41 PM}
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

   >>> dX = dx / 2.
   >>> dY = dy / 2.
   >>> mm = dX * dY / Numeric.sqrt(dX**2 + dY**2) 
   >>> def evalCell(phix, phiy, dx, dy):
   ...     aa = dy**2 + dx**2
   ...     bb = -2 * ( phix * dy**2 + phiy * dx**2)
   ...     cc = dy**2 * phix**2 + dx**2 * phiy**2 - dx**2 * dy**2
   ...     sqr = Numeric.sqrt(bb**2 - 4. * aa * cc)
   ...     return ((-bb - sqr) / 2. / aa,  (-bb + sqr) / 2. / aa)
   >>> v1 = evalCell(dY, dX, dx, dy)[1]
   >>> v2 = max(-dY*3, -dX*3)   
   >>> values = Numeric.array((  v1  ,   dY  ,   dY  ,  dY  ,  v1  ,
   ...                           dX  ,   -mm ,   -dY ,  -mm ,  dX  ,
   ...                           dX  ,   -dX ,   -v1 ,  -dX ,  dX  ,
   ...                           dX  ,   -mm ,   -dY ,  -mm ,  dX  ,
   ...                           v1  ,   dY  ,   dY  ,  dY  ,  v1  ))
   >>> Numeric.allclose(var, values, atol = 1e-10)
   1

"""
__docformat__ = 'restructuredtext'

import Numeric

from fipy.meshes.grid2D import Grid2D
from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable

dx = 1.
dy = 1.
nx = 5
ny = 5

Lx = nx * dx
Ly = ny * dy

mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)

initialArray = -Numeric.ones(nx * ny, 'd')

positiveCells = mesh.getCells(filter = lambda cell: (cell.getCenter()[0] < dx) or (cell.getCenter()[0] > (Lx - dx)) or (cell.getCenter()[1] < dy) or (cell.getCenter()[1] > (Ly - dy)))

for cell in positiveCells:
    initialArray[cell.getID()] = 1.

var = DistanceVariable(
    name = 'level set variable',
    mesh = mesh,
    value = initialArray
    )

var.calcDistanceFunction()

if __name__ == '__main__':
    import fipy.viewers
    viewer = fipy.viewers.make(vars = var, limits = {'datamin': -5., 'datamax': 5.})
    viewer.plot()
    raw_input('finished')
