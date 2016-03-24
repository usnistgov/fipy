#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "circle.py"
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

r"""Solve a circular distance function equation and then advect it.

This example first imposes a circular distance function:

.. math::

   \phi \left( x, y \right) = \left[ \left( x - \frac{ L }{ 2 } \right)^2 + \left( y - \frac{ L }{ 2 } \right)^2 \right]^{1/2} - \frac{L}{4}

The variable is advected with,

.. math::

   \frac{ \partial \phi } { \partial t } + \vec{u} \cdot \nabla \phi = 0

The scheme used in the
:class:`~fipy.terms.firstOrderAdvectionTerm.FirstOrderAdvectionTerm`
preserves the ``var`` as a distance function. The solution to this
problem will be demonstrated in the following script. Firstly, setup
the parameters.

>>> from fipy import CellVariable, Grid2D, DistanceVariable, TransientTerm, FirstOrderAdvectionTerm, AdvectionTerm, Viewer
>>> from fipy.tools import numerix

>>> L = 1.
>>> N = 25
>>> velocity = 1.
>>> cfl = 0.1
>>> velocity = 1.
>>> distanceToTravel = L / 10.
>>> radius = L / 4.
>>> dL = L / N
>>> timeStepDuration = cfl * dL / velocity
>>> steps = int(distanceToTravel / dL / cfl)

Construct the mesh.

>>> mesh = Grid2D(dx=dL, dy=dL, nx=N, ny=N)

Construct a `distanceVariable` object.

>>> var = DistanceVariable(
...     name = 'level set variable',
...     mesh = mesh,
...     value = 1.,
...     hasOld = 1)

Initialise the `distanceVariable` to be a circular distance function.

>>> x, y = mesh.cellCenters
>>> initialArray = numerix.sqrt((x - L / 2.)**2 + (y - L / 2.)**2) - radius
>>> var.setValue(initialArray)

The advection equation is constructed.

>>> advEqn = TransientTerm() + FirstOrderAdvectionTerm(velocity)

The problem can then be solved by executing a serious of time steps.

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var, datamin=-radius, datamax=radius)
...     viewer.plot()
...     for step in range(steps):
...         var.updateOld()
...         advEqn.solve(var, dt=timeStepDuration)
...         viewer.plot()

The result can be tested with the following commands.

>>> for step in range(steps):
...     var.updateOld()
...     advEqn.solve(var, dt=timeStepDuration)
>>> x = numerix.array(mesh.cellCenters[0])
>>> distanceTravelled = timeStepDuration * steps * velocity
>>> answer = initialArray - distanceTravelled
>>> answer = numerix.where(answer < 0., -1001., answer)
>>> solution = numerix.where(answer < 0., -1001., numerix.array(var))
>>> numerix.allclose(answer, solution, atol=4.7e-3)
1

If the advection equation is built with the
:func:`~fipy.terms.advectionTerm.AdvectionTerm`
the result is more accurate,

>>> var.setValue(initialArray)
>>> advEqn = TransientTerm() + AdvectionTerm(velocity)
>>> for step in range(steps):
...     var.updateOld()
...     advEqn.solve(var, dt=timeStepDuration)
>>> solution = numerix.where(answer < 0., -1001., numerix.array(var))
>>> numerix.allclose(answer, solution, atol=1.02e-3)
1

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    input("finished")
