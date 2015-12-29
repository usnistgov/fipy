#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "expandingCircle.py"
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

This example represents an expanding circular interface with an initial
coverage of surfactant. The rate of expansion is dependent on the
coverage of surfactant, The governing equations are given by:

.. math::

   \dot{\theta} &= -\frac{\dot{r}}{r} \theta \\
   \dot{r} &= k \theta

The solution for these set of equations is given by:

.. math::

   r &= \sqrt{2 k r_0 \theta_0 t + r_0^2} \\
   \theta &= \frac{r_0 \theta_0}{\sqrt{2 k r_0 \theta_0 t + r_0^2}}

The following tests can be performed. First test for global
conservation of surfactant:

>>> surfactantBefore = numerix.sum(surfactantVariable * mesh.cellVolumes)
>>> totalTime = 0
>>> steps = 5
>>> for step in range(steps):
...     velocity.setValue(surfactantVariable.interfaceVar * k)
...     distanceVariable.extendVariable(velocity)
...     timeStepDuration = cfl * dx / velocity.max()
...     distanceVariable.updateOld()
...     advectionEquation.solve(distanceVariable, dt = timeStepDuration)
...     surfactantEquation.solve(surfactantVariable, dt=1)
...     totalTime += timeStepDuration #doctest: +LSM
>>> surfactantEquation.solve(surfactantVariable, dt=1)
>>> surfactantAfter = numerix.sum(surfactantVariable * mesh.cellVolumes)
>>> print surfactantBefore.allclose(surfactantAfter)
1

Next test for the correct local value of surfactant:

>>> finalRadius = numerix.sqrt(2 * k * initialRadius * initialSurfactantValue * totalTime + initialRadius**2)
>>> answer = initialSurfactantValue * initialRadius / finalRadius
>>> coverage = surfactantVariable.interfaceVar
>>> error = (coverage / answer - 1)**2 * (coverage > 1e-3)
>>> print numerix.sqrt(numerix.sum(error) / numerix.sum(error > 0)) < 0.04
1

Test for the correct position of the interface:

>>> x, y = mesh.cellCenters
>>> radius = numerix.sqrt((x - L / 2)**2 + (y - L / 2)**2)
>>> solution = radius - distanceVariable
>>> error = (solution / finalRadius - 1)**2 * (coverage > 1e-3)
>>> print numerix.sqrt(numerix.sum(error) / numerix.sum(error > 0)) < 0.02 #doctest: +LSM
1

"""
__docformat__ = 'restructuredtext'

from fipy import CellVariable, SurfactantVariable, Grid2D, DistanceVariable, TransientTerm, ExplicitUpwindConvectionTerm, AdvectionTerm, Viewer
from fipy.tools import numerix

L = 1.
nx = 50
cfl = 0.1
initialRadius = L / 4.
k = 1
dx = L / nx
steps = 20

from fipy.tools import serialComm
mesh = Grid2D(dx = dx, dy = dx, nx = nx, ny = nx, communicator=serialComm)

x, y = mesh.cellCenters
distanceVariable = DistanceVariable(
    name = 'level set variable',
    mesh = mesh,
    value = numerix.sqrt((x - L / 2.)**2 + (y - L / 2.)**2) - initialRadius,
    hasOld = 1)

initialSurfactantValue =  1.

surfactantVariable = SurfactantVariable(
    value = initialSurfactantValue,
    distanceVar = distanceVariable
    )

velocity = CellVariable(
    name = 'velocity',
    mesh = mesh,
    value = 1.,
    )

advectionEquation = TransientTerm() + AdvectionTerm(velocity)

from fipy.variables.surfactantConvectionVariable import SurfactantConvectionVariable
surfactantEquation = TransientTerm() - \
    ExplicitUpwindConvectionTerm(SurfactantConvectionVariable(distanceVariable))

if __name__ == '__main__':

    distanceViewer = Viewer(vars=distanceVariable,
                            datamin=-initialRadius, datamax=initialRadius)
    surfactantViewer = Viewer(vars=surfactantVariable, datamin=0., datamax=100.)
    velocityViewer = Viewer(vars=velocity, datamin=0., datamax=200.)
    distanceViewer.plot()
    surfactantViewer.plot()
    velocityViewer.plot()

    totalTime = 0

    for step in range(steps):
        print 'step',step
        velocity.setValue(surfactantVariable.interfaceVar * k)
        distanceVariable.extendVariable(velocity)
        timeStepDuration = cfl * dx / velocity.max()
        distanceVariable.updateOld()
        advectionEquation.solve(distanceVariable, dt = timeStepDuration)
        surfactantEquation.solve(surfactantVariable, dt=1)

        totalTime += timeStepDuration

        velocityViewer.plot()
        distanceViewer.plot()
        surfactantViewer.plot()

        finalRadius = numerix.sqrt(2 * k * initialRadius * initialSurfactantValue * totalTime + initialRadius**2)
        answer = initialSurfactantValue * initialRadius / finalRadius
        coverage = surfactantVariable.interfaceVar
        error = (coverage / answer - 1)**2 * (coverage > 1e-3)
        print 'error', numerix.sqrt(numerix.sum(error) / numerix.sum(error > 0))




    raw_input('finished')
