#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "square.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

"""

This example advects a 2 by 2 initially square region outwards.
The example checks for global conservation of surfactant.

Advect the interface and check the position.

   >>> distanceVariable.calcDistanceFunction() #doctest: +LSM
   >>> initialSurfactant = numerix.sum(surfactantVariable)
   >>> for step in range(steps):
   ...     distanceVariable.updateOld()
   ...     surfactantEquation.solve(surfactantVariable, dt=1)
   ...     advectionEquation.solve(distanceVariable, dt = timeStepDuration) #doctest: +LSM
   >>> print numerix.allclose(initialSurfactant, numerix.sum(surfactantVariable)) #doctest: +LSM
   1



"""
__docformat__ = 'restructuredtext'

from fipy import CellVariable, SurfactantVariable, Grid2D, DistanceVariable, TransientTerm, ExplicitUpwindConvectionTerm, AdvectionTerm, Viewer
from fipy.tools import numerix

L = 1.
dx = 0.1
velocity = 1.
cfl = 0.1
distanceToTravel = L / 5.
boxSize = .2

nx = int(L / dx)
ny = int(L / dx)

steps = int(distanceToTravel / dx / cfl)

timeStepDuration = cfl * dx / velocity

mesh = Grid2D(dx = dx, dy = dx, nx = nx, ny = ny)

x0 = (L - boxSize) / 2
x1 = (L + boxSize) / 2

distanceVariable = DistanceVariable(
    mesh = mesh,
    value = 1.,
    hasOld = 1
    )

x, y = mesh.cellCenters
distanceVariable.setValue(-1, where=((x0 < x) & (x < x1)) & ((x0 < y) & (y < x1)))


surfactantVariable = SurfactantVariable(
    distanceVar = distanceVariable,
    value = 1.
    )


from fipy.variables.surfactantConvectionVariable import SurfactantConvectionVariable
surfactantEquation = TransientTerm() - \
    ExplicitUpwindConvectionTerm(SurfactantConvectionVariable(distanceVariable))

advectionEquation = TransientTerm() + AdvectionTerm(velocity)

if __name__ == '__main__':
    distanceViewer = Viewer(vars=distanceVariable,
                            datamin=-.001, datamax=.001)
    surfactantViewer = Viewer(vars=surfactantVariable,
                              datamin=0., datamax=2.)


    distanceVariable.calcDistanceFunction()

    for step in range(steps):
        print numerix.sum(surfactantVariable)
        distanceVariable.updateOld()
        surfactantEquation.solve(surfactantVariable, dt=1)
        advectionEquation.solve(distanceVariable, dt = timeStepDuration)
        distanceViewer.plot()
        surfactantViewer.plot()

    surfactantEquation.solve(surfactantVariable, dt=1)

    distanceViewer.plot()
    surfactantViewer.plot()
    print surfactantVariable
    raw_input('finished')
