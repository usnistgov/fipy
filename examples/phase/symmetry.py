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
>>> from builtins import range
>>> for j in range(N // 2):
...     for i in range(N // 2):
...         x = dx * (i + 0.5)
...         y = dx * (j + 0.5)
...         testResult[i, j] = x * y
...         bottomRight[i, j] = var(((L - x,), (y,)))[0]
...         topLeft[i, j] = var(((x,), (L - y,)))[0]
...         topRight[i, j] = var(((L - x,), (L - y,)))[0]
>>> numerix.allclose(testResult, bottomRight, atol = 1e-10)
1
>>> numerix.allclose(testResult, topLeft, atol = 1e-10)
1
>>> numerix.allclose(testResult, topRight, atol = 1e-10)
1
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    input("finished")
