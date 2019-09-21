r"""Create a level set variable in one dimension.

The level set
variable calculates its value over the domain to be the distance from
the zero level set. This can be represented succinctly in the
following equation with a boundary condition at the zero level set
such that,

.. math::

   \frac{\partial \phi}{\partial x} = 1

with the boundary condition, :math:`\phi = 0` at :math:`x = L / 2`.

The solution to this problem will be demonstrated in the following
script. Firstly, setup the parameters.

>>> from fipy import CellVariable, Grid1D, DistanceVariable, TransientTerm, FirstOrderAdvectionTerm, AdvectionTerm, Viewer
>>> from fipy.tools import numerix, serialComm

>>> dx = 0.5
>>> nx = 10

Construct the mesh.

.. index::
   single: Grid2D

>>> mesh = Grid1D(dx=dx, nx=nx, communicator=serialComm)

Construct a `distanceVariable` object.

>>> var = DistanceVariable(name='level set variable',
...                        mesh=mesh,
...                        value=-1.,
...                        hasOld=1)
>>> x = mesh.cellCenters[0]
>>> var.setValue(1, where=x > dx * nx / 2)

Once the initial positive and negative regions have been initialized
the `calcDistanceFunction()` method can be used to recalculate `var`
as a distance function from the zero level set.

>>> var.calcDistanceFunction() #doctest: +LSM

The problem can then be solved by executing the :meth:`~fipy.terms.term.Term.solve`
method of the equation.

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var, datamin=-5., datamax=5.)
...     viewer.plot()

The result can be tested with the following commands.

>>> print(numerix.allclose(var, x - dx * nx / 2)) #doctest: +LSM
1

"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    input("finished")

