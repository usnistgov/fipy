#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 4/2/04 {4:01:07 PM} { 1:23:41 PM}
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

"""

This example follows the advection of a trench as in example
`examples/levelSet/advection/input.py`. In this example there is
a surfactant on the interface. We wish to 


Advect the interface and check the position.

   >>> for step in range(steps):
   ...     it.timestep(dt = timeStepDuration)

   >>> distanceMoved = timeStepDuration * steps * velocity
   >>> answer = answer - distanceMoved
   >>> solution = Numeric.array(distanceVariable)
   >>> answer = Numeric.where(answer < 0., 0., answer)
   >>> solution = Numeric.where(solution < 0., 0., solution)
   >>> Numeric.allclose(answer, solution, atol = 1e-1)
   1
 

   
"""
__docformat__ = 'restructuredtext'

import Numeric
   
from fipy.meshes.grid2D import Grid2D
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.variables.cellVariable import CellVariable
from fipy.models.levelSet.distanceFunction.distanceFunctionEquation import DistanceFunctionEquation
from fipy.models.levelSet.advection.advectionEquation import AdvectionEquation
from fipy.models.levelSet.surfactant.surfactantEquation import SurfactantEquation
from fipy.models.levelSet.advection.higherOrderAdvectionTerm import HigherOrderAdvectionTerm
from fipy.iterators.iterator import Iterator
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.solvers.linearLUSolver import LinearLUSolver
from fipy.boundaryConditions.fixedValue import FixedValue

L = 1.
dx = 0.02
velocity = 1.
cfl = 0.1
distanceToTravel = L / 2.

nx = int(L / dx)
ny = int(L / dx)
steps = int(distanceToTravel / dx / cfl)

timeStepDuration = cfl * dx / velocity

mesh = Grid2D(dx = dx, dy = dx, nx = nx, ny = ny)

distanceVariable = CellVariable(
    name = 'level set variable',
    mesh = mesh,
    value = -1.
    )

surfactantVariable = CellVariable(
    name = 'surfactant variable',
    mesh = mesh,
    value = 1.
    )

positiveCells = mesh.getCells(lambda cell: (cell.getCenter()[1] > L / 2) or (cell.getCenter()[0] > L / 2 and cell.getCenter()[1] > L / 4))

distanceVariable.setValue(1., positiveCells)

distanceEquation = DistanceFunctionEquation(distanceVariable)

surfactantEquation = SurfactantEquation(
    surfactantVariable,
    distanceVariable,
    solver = LinearLUSolver(
        tolerance = 1e-10),
    boundaryConditions = (FixedValue(mesh.getExteriorFaces(), 0),))

advectionEquation = AdvectionEquation(
    distanceVariable,
    advectionCoeff = velocity,
    solver = LinearPCGSolver(
        tolerance = 1.e-15, 
        steps = 1000),
    advectionTerm = HigherOrderAdvectionTerm)

it = Iterator((surfactantEquation, advectionEquation))

def evaluateCoverage(phi, theta):
    alpha = 1.5 * dx
    delta = Numeric.where(phi > -alpha,
                          Numeric.where(phi < alpha,
                                        (1 - 0.5 * Numeric.cos(Numeric.pi * phi / alpha)) / (2 * alpha),
                                        0.),
                          0.)
    return Numeric.sum(delta * phi.getGrad().getMag() * theta * mesh.getCellVolumes())

if __name__ == '__main__':
    distanceViewer = Grid2DGistViewer(var = distanceVariable, palette = 'rainbow.gp', minVal = -.001, maxVal = .001)
    surfactantViewer = Grid2DGistViewer(var = surfactantVariable, palette = 'rainbow.gp', minVal = 0., maxVal = 2.)

    distanceViewer.plot()
    surfactantViewer.plot()

    distanceEquation.solve()

    for step in range(steps):
        print evaluateCoverage(distanceVariable, surfactantVariable)
        it.timestep(dt = timeStepDuration)
        print 'step:',step

        distanceViewer.plot()
        surfactantViewer.plot()

    raw_input('finished')
