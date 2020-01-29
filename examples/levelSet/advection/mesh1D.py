r"""Solve the distance function equation in one dimension and then advect it.

This example first solves the distance function equation in one dimension:

.. math::

   \abs{\nabla \phi} = 1

with :math:`\phi = 0` at :math:`x = L / 5`.

The variable is then advected with,

.. math::

   \frac{ \partial \phi } { \partial t} + \vec{u} \cdot \nabla \phi = 0

The scheme used in the `FirstOrderAdvectionTerm` preserves the `var` as a distance function.

The solution to this problem will be demonstrated in the following
script. Firstly, setup the parameters.

>>> from fipy import CellVariable, Grid1D, DistanceVariable, TransientTerm, FirstOrderAdvectionTerm, AdvectionTerm, Viewer
>>> from fipy.tools import numerix, serialComm

>>> velocity = 1.
>>> dx = 1.
>>> nx = 10
>>> timeStepDuration = 1.
>>> steps = 2
>>> L = nx * dx
>>> interfacePosition = L / 5.

Construct the mesh.

.. index::
   single: Grid1D

>>> mesh = Grid1D(dx=dx, nx=nx, communicator=serialComm)

Construct a `distanceVariable` object.

>>> var = DistanceVariable(name='level set variable',
...                        mesh=mesh,
...                        value=-1.,
...                        hasOld=1)
>>> var.setValue(1., where=mesh.cellCenters[0] > interfacePosition)
>>> var.calcDistanceFunction() #doctest: +LSM

The `advectionEquation` is constructed.

>>> advEqn = TransientTerm() + FirstOrderAdvectionTerm(velocity)

The problem can then be solved by executing a serious of time steps.

>>> from builtins import range
>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var, datamin=-10., datamax=10.)
...     viewer.plot()
...     for step in range(steps):
...         var.updateOld()
...         advEqn.solve(var, dt=timeStepDuration)
...         viewer.plot()

The result can be tested with the following code:

>>> from builtins import range
>>> for step in range(steps):
...     var.updateOld()
...     advEqn.solve(var, dt=timeStepDuration)
>>> x = mesh.cellCenters[0]
>>> distanceTravelled = timeStepDuration * steps * velocity
>>> answer = x - interfacePosition - timeStepDuration * steps * velocity
>>> answer = numerix.where(x < distanceTravelled,
...                        x[0] - interfacePosition, answer)
>>> print(var.allclose(answer)) #doctest: +LSM
1

"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    input("finished")
