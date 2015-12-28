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

r"""
This example first imposes a circular distance function:

.. math::

   \phi \left( x, y \right) = \left[ \left( x - \frac{ L }{ 2 } \right)^2 + \left( y - \frac{ L }{ 2 } \right)^2 \right]^{1/2} - \frac{L}{4}

then the variable is advected with,

.. math::

   \frac{ \partial \phi } { \partial t } + \vec{u} \cdot \nabla \phi = 0

Also a surfactant is present of the interface, governed by the equation:

.. math::

   \frac{d \theta}{d t} = J v \theta

The result can be tested with the following code:


>>> surfactantBefore = numerix.sum(surfactantVariable * mesh.cellVolumes)
>>> for step in range(steps):
...     distanceVariable.updateOld()
...     surfactantEquation.solve(surfactantVariable, dt=1.)
...     advectionEquation.solve(distanceVariable, dt = timeStepDuration)
>>> surfactantEquation.solve(surfactantVariable, dt=1.)
>>> surfactantAfter = numerix.sum(surfactantVariable * mesh.cellVolumes)
>>> print surfactantBefore.allclose(surfactantAfter)
1
>>> areas = (distanceVariable.cellInterfaceAreas < 1e-6) * 1e+10 + distanceVariable.cellInterfaceAreas
>>> answer = initialSurfactantValue * initialRadius / (initialRadius +  distanceToTravel)
>>> coverage = surfactantVariable * mesh.cellVolumes / areas
>>> error = (coverage / answer - 1)**2 * (coverage > 1e-3)
>>> print numerix.sqrt(numerix.sum(error) / numerix.sum(error > 0))
0.00813776069241

"""
__docformat__ = 'restructuredtext'

from fipy import CellVariable, SurfactantVariable, Grid2D, DistanceVariable, TransientTerm, ExplicitUpwindConvectionTerm, AdvectionTerm, Viewer
from fipy.tools import numerix

L = 1.
nx = 50
velocity = 1.
cfl = 0.1
velocity = 1.
distanceToTravel = L / 10.
initialRadius = L / 4.

dx = L / nx
timeStepDuration = cfl * dx / velocity
steps = int(distanceToTravel / dx / cfl)

mesh = Grid2D(dx = dx, dy = dx, nx = nx, ny = nx)

distanceVariable = DistanceVariable(
    name = 'level set variable',
    mesh = mesh,
    value = 1.,
    hasOld = 1
    )

x, y = mesh.cellCenters
cellRadius = numerix.sqrt((x - L / 2.)**2 + (y - L / 2.)**2)
distanceVariable.setValue(cellRadius - initialRadius)

initialSurfactantValue =  1.

surfactantVariable = SurfactantVariable(
    value = initialSurfactantValue,
    distanceVar = distanceVariable
    )

advectionEquation = TransientTerm() + AdvectionTerm(velocity)

from fipy.variables.surfactantConvectionVariable import SurfactantConvectionVariable
surfactantEquation = TransientTerm() - \
    ExplicitUpwindConvectionTerm(SurfactantConvectionVariable(distanceVariable))

if __name__ == '__main__':

    distanceViewer = Viewer(vars=distanceVariable,
                            datamin=-initialRadius, datamax=initialRadius)
    surfactantViewer = Viewer(vars=surfactantVariable, datamin=-1., datamax=100.)
    distanceViewer.plot()
    surfactantViewer.plot()

    print 'total surfactant before:', numerix.sum(surfactantVariable * mesh.cellVolumes)

    for step in range(steps):
        distanceVariable.updateOld()
        surfactantEquation.solve(surfactantVariable, dt=1.)
        advectionEquation.solve(distanceVariable, dt = timeStepDuration)
        distanceViewer.plot()
        surfactantViewer.plot()
    surfactantEquation.solve(surfactantVariable, dt=1.)


    print 'total surfactant after:', numerix.sum(surfactantVariable * mesh.cellVolumes)

    areas = (distanceVariable.cellInterfaceAreas < 1e-6) * 1e+10 + distanceVariable.cellInterfaceAreas
    answer = initialSurfactantValue * initialRadius / (initialRadius +  distanceToTravel)
    coverage = surfactantVariable * mesh.cellVolumes / areas

    error = 0.
    size = 0
    for i in range(len(coverage)):
        if coverage[i] > 1e-3:
            error += (coverage[i] / answer - 1.)**2
            size += 1

    error = numerix.sqrt(error / size)

    print 'error:', error

    raw_input('finished')
