#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 10/14/04 {5:13:24 PM} { 1:23:41 PM}
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

r"""

Here we solve the level set equation in two dimensions for a circle. The 2D
level set equation can be written,

.. raw:: latex

    $$ | \nabla \phi | = 1$$

and the boundary condition for a circle is given by,

.. raw:: latex

    $ \phi = 0 $ at $(x - L / 2)^2 + (y - L / 2)^2 = (L / 4)^2 $.

The solution to this problem will be demonstrated in the following
script. Firstly, setup the parameters.

   >>> dx = 1.
   >>> dy = 1.
   >>> nx = 11
   >>> ny = 11
   >>> Lx = nx * dx
   >>> Ly = ny * dy

Construct the mesh.

   >>> from fipy.meshes.grid2D import Grid2D
   >>> mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)

Construct a `distanceVariable` object.

   >>> from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
   >>> import Numeric
   >>> initialArray = -Numeric.ones(nx * ny, 'd')
   >>> positiveCells = mesh.getCells(filter = lambda cell:
   ...                   (cell.getCenter()[0] - Lx / 2.)**2 +
   ...                   (cell.getCenter()[1] - Ly / 2.)**2 <
   ...                   (Lx / 4.)**2)
   >>> for cell in positiveCells:
   ...     initialArray[cell.getID()] = 1.
   
   >>> var = DistanceVariable(name = 'level set variable',
   ...                        mesh = mesh,
   ...                        value = initialArray)
   >>> var.calcDistanceFunction()
   
   >>> if __name__ == '__main__':
   ...     from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
   ...     viewer = Grid2DGistViewer(var = var, palette = 'rainbow.gp', minVal = -5., maxVal = 5.)
   ...     viewer.plot()
   ...     viewer.plot()

The result can be tested with the following commands.

   >>> dY = dy / 2.
   >>> dX = dx / 2.
   >>> mm = min (dX, dY)
   >>> m1 = dY * dX / Numeric.sqrt(dY**2 + dX**2)
   >>> def evalCell(phix, phiy, dx, dy):
   ...     aa = dy**2 + dx**2
   ...     bb = -2 * ( phix * dy**2 + phiy * dx**2)
   ...     cc = dy**2 * phix**2 + dx**2 * phiy**2 - dx**2 * dy**2
   ...     sqr = Numeric.sqrt(bb**2 - 4. * aa * cc)
   ...     return ((-bb - sqr) / 2. / aa,  (-bb + sqr) / 2. / aa)
   >>> v1 = evalCell(-dY, -m1, dx, dy)[0] 
   >>> v2 = evalCell(-m1, -dX, dx, dy)[0]
   >>> v3 = evalCell(m1,  m1,  dx, dy)[1]
   >>> v4 = evalCell(v3, dY, dx, dy)[1]
   >>> v5 = evalCell(dX, v3, dx, dy)[1]
   >>> import MA
   >>> trialValues = MA.masked_values((
   ...     -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,
   ...     -1000, -1000, -1000, -1000, -3*dY, -3*dY, -3*dY, -1000, -1000, -1000, -1000,
   ...     -1000, -1000, -1000, v1   , -dY  , -dY  , -dY  , v1   , -1000, -1000, -1000,
   ...     -1000, -1000, v2   , -m1  , m1   , dY   , m1   , -m1  , v2   , -1000, -1000,
   ...     -1000, -dX*3, -dX  , m1   ,  v3  , v4   , v3   , m1   , -dX  , -dX*3, -1000,
   ...     -1000, -dX*3, -dX  , dX   , v5   , -1000, v5   , dX   , -dX  , -dX*3, -1000,
   ...     -1000, -dX*3, -dX  , m1   , v3   , v4   , v3   , m1   , -dX  , -dX*3, -1000,
   ...     -1000, -1000, v2   , -m1  , m1   , dY   , m1   , -m1  , v2   , -1000, -1000,
   ...     -1000, -1000, -1000, v1   , -dY  , -dY  , -dY  , v1   , -1000, -1000, -1000,
   ...     -1000, -1000, -1000, -1000, -3*dY, -3*dY, -3*dY, -1000, -1000, -1000, -1000,
   ...     -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000), 
   ...     -1000)
   >>> MA.allclose(var, trialValues)
   1
   
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())
    raw_input("finished")
