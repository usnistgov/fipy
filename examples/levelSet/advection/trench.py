r"""

This example creates a trench with the following zero level set:

.. math::

   \phi \left( x, y \right) = 0 & \text{when $y = L_y / 5$ and $x \ge L_x / 2$} \\
   \phi \left( x, y \right) = 0 & \text{when $L_y / 5 \le y \le 3 Ly / 5$ and $x = L_x / 2$$} \\
   \phi \left( x, y \right) = 0 & \text{when $y = 3 Ly / 5$ and $x \le L_x / 2$}

>>> from fipy import CellVariable, Grid2D, DistanceVariable, TransientTerm, FirstOrderAdvectionTerm, AdvectionTerm, Viewer
>>> from fipy.tools import numerix, serialComm

>>> height = 0.5
>>> Lx = 0.4
>>> Ly = 1.
>>> dx = 0.01
>>> velocity = 1.
>>> cfl = 0.1

>>> nx = int(Lx / dx)
>>> ny = int(Ly / dx)
>>> timeStepDuration = cfl * dx / velocity
>>> steps = 200

>>> mesh = Grid2D(dx = dx, dy = dx, nx = nx, ny = ny, communicator=serialComm)

>>> var = DistanceVariable(name = 'level set variable',
...                        mesh = mesh,
...                        value = -1.,
...                        hasOld = 1
...                        )

>>> x, y = mesh.cellCenters
>>> var.setValue(1, where=(y > 0.6 * Ly) | ((y > 0.2 * Ly) & (x > 0.5 * Lx)))

>>> var.calcDistanceFunction() #doctest: +LSM

>>> advEqn = TransientTerm() + FirstOrderAdvectionTerm(velocity)

The trench is then advected with a unit velocity. The following test can be made
for the initial position of the interface:

>>> r1 =  -numerix.sqrt((x - Lx / 2)**2 + (y - Ly / 5)**2)
>>> r2 =  numerix.sqrt((x - Lx / 2)**2 + (y - 3 * Ly / 5)**2)
>>> d = numerix.zeros((len(x), 3), 'd')
>>> d[:, 0] = numerix.where(x >= Lx / 2, y - Ly / 5, r1)
>>> d[:, 1] = numerix.where(x <= Lx / 2, y - 3 * Ly / 5, r2)
>>> d[:, 2] = numerix.where(numerix.logical_and(Ly / 5 <= y, y <= 3 * Ly / 5), x - Lx / 2, d[:, 0])
>>> argmins = numerix.argmin(numerix.absolute(d), axis = 1)
>>> answer = numerix.take(d.ravel(), numerix.arange(len(argmins))*3 + argmins)
>>> print(var.allclose(answer, atol = 1e-1)) #doctest: +LSM
1

Advect the interface and check the position.

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var, datamin=-0.1, datamax=0.1)
... 
...     viewer.plot()

>>> from builtins import range
>>> for step in range(steps):
...     var.updateOld()
...     advEqn.solve(var, dt = timeStepDuration)
...     if __name__ == '__main__':
...         viewer.plot()

>>> distanceMoved = timeStepDuration * steps * velocity
>>> answer = answer - distanceMoved
>>> answer = numerix.where(answer < 0., 0., answer)
>>> var.setValue(numerix.where(var < 0., 0., var))
>>> print(var.allclose(answer, atol = 1e-1)) #doctest: +LSM
1

"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    input('finished')
