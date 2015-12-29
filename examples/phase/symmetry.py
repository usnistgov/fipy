#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "input.py"
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

This example creates four symmetric quadrilateral regions in a box.
We start with a :class:`~fipy.variables.cellVariable.CellVariable` object that contains the following
values:

.. math::

   \phi(x, y) = x y,
   0 \le x \le L,
   0 \le y \le L

We wish to create 4 symmetric regions such that

.. math::

   \phi(x, y) = \phi(L - x, y) = \phi(L - x, y) = \phi(L - x, L - y),
   0 \le x \le L / 2,
   0 \le y \le L / 2

We create a square domain

>>> from fipy import CellVariable, Grid2D, Viewer
>>> from fipy.tools import numerix

>>> N = 20
>>> L = 1.
>>> dx = L / N
>>> dy = L / N

>>> mesh = Grid2D(
...    dx = dx,
...    dy = dy,
...    nx = N,
...    ny = N)

>>> var = CellVariable(name = "test", mesh = mesh)

First set the values as given in the above equation:

>>> x, y = mesh.cellCenters
>>> var.setValue(x * y)

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var, datamin=0, datamax=L * L / 4.)
...     viewer.plot()

The bottom-left quadrant is mirrored into each of the other three quadrants

>>> q = (x > L / 2.) & (y < L / 2.)
>>> var[q] = var(((L - x)[q],       y[q]))
>>> q = (x < L / 2.) & (y > L / 2.)
>>> var[q] = var((      x[q], (L - y)[q]))
>>> q = (x > L / 2.) & (y > L / 2.)
>>> var[q] = var(((L - x)[q], (L - y)[q]))

>>> if __name__ == '__main__':
...     viewer.plot()

The following code tests the results with a different algorithm:

>>> testResult = numerix.zeros((N // 2, N // 2), 'd')
>>> bottomRight = numerix.zeros((N // 2, N // 2), 'd')
>>> topLeft = numerix.zeros((N // 2, N // 2), 'd')
>>> topRight = numerix.zeros((N // 2, N // 2), 'd')
>>> for j in range(N // 2):
...     for i in range(N // 2):
...         x = dx * (i + 0.5)
...         y = dx * (j + 0.5)
...         testResult[i, j] = x * y
...         bottomRight[i,j] = var(((L - x,), (y,)))
...         topLeft[i,j] = var(((x,), (L - y,)))
...         topRight[i,j] = var(((L - x,), (L - y,)))
>>> numerix.allclose(testResult, bottomRight, atol = 1e-10)
1
>>> numerix.allclose(testResult,topLeft, atol = 1e-10)
1
>>> numerix.allclose(testResult,topRight, atol = 1e-10)
1
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    raw_input("finished")
