#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "interior.py"
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

r"""

Here we solve the level set equation in two dimension for an interior region. The equation is
given by:

.. math::

   \abs{\nabla \phi} $= 1 \\
   \phi &= 0 \;\; \text{at} \begin{cases}
    \text{x = \left( d, L - d \right)$ for $d \le y \le L - d$} \\
    \text{y = \left( d, L - d \right)$ for $d \le x \le L - d$}
    \end{cases}

Do the tests:

>>> var.calcDistanceFunction(order=1) #doctest: +LSM

>>> dX = dx / 2.
>>> dY = dy / 2.
>>> mm = dX * dY / numerix.sqrt(dX**2 + dY**2)
>>> def evalCell(phix, phiy, dx, dy):
...     aa = dy**2 + dx**2
...     bb = -2 * ( phix * dy**2 + phiy * dx**2)
...     cc = dy**2 * phix**2 + dx**2 * phiy**2 - dx**2 * dy**2
...     sqr = numerix.sqrt(bb**2 - 4. * aa * cc)
...     return ((-bb - sqr) / 2. / aa,  (-bb + sqr) / 2. / aa)
>>> v1 = evalCell(dY, dX, dx, dy)[1]
>>> v2 = max(-dY*3, -dX*3)
>>> values = numerix.array((  v1  ,   dY  ,   dY  ,  dY  ,  v1  ,
...                           dX  ,  -mm ,   -dY ,  -mm ,   dX  ,
...                           dX  ,  -dX ,   -v1 ,  -dX ,   dX  ,
...                           dX  ,  -mm ,   -dY ,  -mm ,   dX  ,
...                           v1  ,   dY  ,   dY  ,  dY  ,  v1  ))
>>> print var.allclose(values, atol = 1e-10) #doctest: +LSM
1

"""
__docformat__ = 'restructuredtext'

from fipy import Grid2D, DistanceVariable, Viewer
from fipy.tools import numerix, serialComm

dx = 1.
dy = 1.
nx = 5
ny = 5

Lx = nx * dx
Ly = ny * dy

mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny, communicator=serialComm)

var = DistanceVariable(
    name = 'level set variable',
    mesh = mesh,
    value = -1.,
    hasOld = 1
    )

x, y = mesh.cellCenters
var.setValue(1, where=((x < dx) | (x > (Lx - dx))
                       | (y < dy) | (y > (Ly - dy))))

if __name__ == '__main__':
    var.calcDistanceFunction(order=1)
    viewer = Viewer(vars=var, datamin=-5., datamax=5.)
    viewer.plot()
    raw_input('finished')
