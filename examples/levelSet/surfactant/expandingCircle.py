#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "expandingCircle.py"
 #                                    created: 08/10/04 {10:29:10 AM} 
 #                                last update: 08/10/04 {6:09:38 PM} { 1:23:41 PM}
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

This example represents an expanding circular interface with an initial
coverage of surfactant. The rate of expansion is dependent on the
coverage of surfactant, The governing equations are given by:

.. raw:: latex

    $$ \\dot{\\theta} = -\\frac{\\dot{r}}{r} \\theta $$
    $$ \\dot{r} = k \\theta $$

The solution for these set of equations is given by:

.. raw:: latex

    $$ r = \\sqrt{2 k r_0 \\theta_0 t + r_0^2} $$
    $$ \\theta = \\frac{r_0 \\theta_0}{\\sqrt{2 k r_0 \\theta_0 t + r_0^2}} $$
    
The following tests can be performed. First test for global
conservation of surfactant:

   >>> surfactantBefore = Numeric.sum(surfactantVariable * mesh.getCellVolumes())
   >>> totalTime = 0
   >>> for step in range(steps):
   ...     velocity.setValue(surfactantVariable.getInterfaceValue() * k)
   ...     argmax = Numeric.argmax(velocity)
   ...     timeStepDuration = cfl * dx / velocity[argmax]
   ...     it.timestep(dt = timeStepDuration)
   ...     totalTime += timeStepDuration
   >>> surfactantEquation.solve()
   >>> surfactantAfter = Numeric.sum(surfactantVariable * mesh.getCellVolumes())
   >>> Numeric.allclose(surfactantBefore, surfactantAfter)
   1

Next test for the correct local value of surfactant: 

   >>> finalRadius = Numeric.sqrt(2 * k * initialRadius * initialSurfactantValue * totalTime + initialRadius**2)
   >>> answer = initialSurfactantValue * initialRadius / finalRadius
   >>> coverage = surfactantVariable.getInterfaceValue()
   >>> error = 0.
   >>> size = 0
   >>> for i in range(len(coverage)):
   ...     if coverage[i] > 1e-3:
   ...         error += (coverage[i] / answer - 1.)**2
   ...         size += 1
   >>> print Numeric.sqrt(error / size) < 0.04
   1

Test for the correct position of the interface:

   >>> x = mesh.getCellCenters()[:,0]
   >>> y = mesh.getCellCenters()[:,1]
   >>> radius = Numeric.sqrt((x - L / 2)**2 + (y - L / 2)**2)
   >>> solution = radius - distanceVariable
   >>> error = 0.
   >>> size = 0
   >>> for i in range(len(coverage)):
   ...     if coverage[i] > 1e-3:
   ...         error += (solution[i] / finalRadius - 1.)**2
   ...         size += 1
   >>> print Numeric.sqrt(error / size) < 0.02
   1

"""
__docformat__ = 'restructuredtext'

import Numeric
   
from fipy.meshes.grid2D import Grid2D
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.models.levelSet.distanceFunction.extensionEquation import ExtensionEquation
from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
from fipy.models.levelSet.advection.advectionEquation import AdvectionEquation
from fipy.models.levelSet.advection.higherOrderAdvectionTerm import HigherOrderAdvectionTerm
from fipy.models.levelSet.surfactant.surfactantEquation import SurfactantEquation
from fipy.models.levelSet.surfactant.surfactantVariable import SurfactantVariable
from fipy.iterators.iterator import Iterator
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.solvers.linearLUSolver import LinearLUSolver
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.variables.cellVariable import CellVariable

L = 1.
nx = 50
cfl = 0.1
initialRadius = L / 4.
k = 1
dx = L / nx
steps = 20

mesh = Grid2D(dx = dx, dy = dx, nx = nx, ny = nx)

distanceVariable = DistanceVariable(
    name = 'level set variable',
    mesh = mesh,
    value = 1.
    )

cellRadius = Numeric.sqrt((mesh.getCellCenters()[:,0] - L / 2.)**2 + (mesh.getCellCenters()[:,1] - L / 2.)**2)
distanceVariable.setValue(cellRadius - initialRadius)

initialSurfactantValue =  1.

surfactantVariable = SurfactantVariable(
    value = initialSurfactantValue,
    distanceVariable = distanceVariable
    )

velocity = CellVariable(
    name = 'velocity',
    mesh = mesh,
    value = 1.,
    )

advectionEquation = AdvectionEquation(
    distanceVariable,
    advectionCoeff = velocity,
    solver = LinearPCGSolver(
        tolerance = 1.e-15, 
        steps = 1000),
    advectionTerm = HigherOrderAdvectionTerm)

surfactantEquation = SurfactantEquation(
    surfactantVariable,
    distanceVariable,
    solver = LinearLUSolver(
        tolerance = 1e-10),
    boundaryConditions = (FixedValue(mesh.getExteriorFaces(), 0.), ))

extensionEquation = ExtensionEquation(
    distanceVariable,
    velocity)

it = Iterator((extensionEquation, advectionEquation, surfactantEquation))

if __name__ == '__main__':
    
    distanceViewer = Grid2DGistViewer(var = distanceVariable, palette = 'rainbow.gp', minVal = -initialRadius, maxVal = initialRadius)
    surfactantViewer = Grid2DGistViewer(var = surfactantVariable, palette = 'rainbow.gp', minVal = 0., maxVal = 100.)
    velocityViewer = Grid2DGistViewer(var = velocity, palette = 'rainbow.gp', minVal = 0., maxVal = 200.)
    distanceViewer.plot()
    surfactantViewer.plot()
    velocityViewer.plot()

    totalTime = 0

    for step in range(steps):
        print 'step',step
        velocity.setValue(surfactantVariable.getInterfaceValue() * k)

        argmax = Numeric.argmax(velocity)
        timeStepDuration = cfl * dx / velocity[argmax]
        
        it.timestep(dt = timeStepDuration)

        totalTime += timeStepDuration
        
        velocityViewer.plot()
        distanceViewer.plot()
        surfactantViewer.plot(minVal = 0, maxVal = 100)

        finalRadius = Numeric.sqrt(2 * k * initialRadius * initialSurfactantValue * totalTime + initialRadius**2)
        answer = initialSurfactantValue * initialRadius / finalRadius
        coverage = surfactantVariable.getInterfaceValue()
        error = 0.
        size = 0
        for i in range(len(coverage)):
            if coverage[i] > 1e-3:
                error += (coverage[i] / answer - 1.)**2
                size += 1

        print 'error',Numeric.sqrt(error / size)


        

    raw_input('finished')
