#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "expandingCircle.py"
 #                                    created: 08/10/04 {10:29:10 AM} 
 #                                last update: 8/2/05 {5:04:15 PM} { 1:23:41 PM}
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

   >>> surfactantBefore = numerix.sum(surfactantVariable * mesh.getCellVolumes())
   >>> totalTime = 0
   >>> for step in range(steps):
   ...     velocity.setValue(surfactantVariable.getInterfaceVar() * k)
   ...     distanceVariable.extendVariable(velocity)
   ...     timeStepDuration = cfl * dx / numerix.max(velocity)
   ...     distanceVariable.updateOld()
   ...     advectionEquation.solve(distanceVariable, dt = timeStepDuration)
   ...     surfactantEquation.solve(surfactantVariable)
   ...     totalTime += timeStepDuration
   >>> surfactantEquation.solve(surfactantVariable)
   >>> surfactantAfter = numerix.sum(surfactantVariable * mesh.getCellVolumes())
   >>> print surfactantBefore.allclose(surfactantAfter)
   1

Next test for the correct local value of surfactant: 

   >>> finalRadius = numerix.sqrt(2 * k * initialRadius * initialSurfactantValue * totalTime + initialRadius**2)
   >>> answer = initialSurfactantValue * initialRadius / finalRadius
   >>> coverage = surfactantVariable.getInterfaceVar()
   >>> error = 0.
   >>> size = 0
   >>> for i in range(len(coverage)):
   ...     if coverage[i] > 1e-3:
   ...         error += (coverage[i] / answer - 1.)**2
   ...         size += 1
   >>> print numerix.sqrt(error / size) < 0.04
   1

Test for the correct position of the interface:

   >>> x = mesh.getCellCenters()[0]
   >>> y = mesh.getCellCenters()[1]
   >>> radius = numerix.sqrt((x - L / 2)**2 + (y - L / 2)**2)
   >>> solution = radius - distanceVariable
   >>> error = 0.
   >>> size = 0
   >>> for i in range(len(coverage)):
   ...     if coverage[i] > 1e-3:
   ...         error += (solution[i] / finalRadius - 1.)**2
   ...         size += 1
   >>> print numerix.sqrt(error / size) < 0.02
   1

"""
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
   
from fipy.meshes.grid2D import Grid2D
from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
from fipy.models.levelSet.advection.higherOrderAdvectionEquation import buildHigherOrderAdvectionEquation
from fipy.models.levelSet.surfactant.surfactantEquation import SurfactantEquation
from fipy.models.levelSet.surfactant.surfactantVariable import SurfactantVariable
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
    value = numerix.sqrt((mesh.getCellCenters()[0] - L / 2.)**2 + (mesh.getCellCenters()[1] - L / 2.)**2) - initialRadius,
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

advectionEquation = buildHigherOrderAdvectionEquation(
    advectionCoeff = velocity)

surfactantEquation = SurfactantEquation(
    distanceVar = distanceVariable)

if __name__ == '__main__':
    
    import fipy.viewers
    distanceViewer = fipy.viewers.make(vars = distanceVariable, limits = {'datamin': -initialRadius, 'datamax': initialRadius})
    surfactantViewer = fipy.viewers.make(vars = surfactantVariable, limits = {'datamin': 0., 'datamax': 100.})
    velocityViewer = fipy.viewers.make(vars = velocity, limits = {'datamin': 0., 'datamax': 200.})
    distanceViewer.plot()
    surfactantViewer.plot()
    velocityViewer.plot()

    totalTime = 0

    for step in range(steps):
        print 'step',step
        velocity.setValue(surfactantVariable.getInterfaceVar() * k)
        distanceVariable.extendVariable(velocity)
        timeStepDuration = cfl * dx / numerix.max(velocity)
        distanceVariable.updateOld()
        advectionEquation.solve(distanceVariable, dt = timeStepDuration)
        surfactantEquation.solve(surfactantVariable)
        
        totalTime += timeStepDuration
        
        velocityViewer.plot()
        distanceViewer.plot()
        surfactantViewer.plot()

        finalRadius = numerix.sqrt(2 * k * initialRadius * initialSurfactantValue * totalTime + initialRadius**2)
        answer = initialSurfactantValue * initialRadius / finalRadius
        coverage = surfactantVariable.getInterfaceVar()
        error = 0.
        size = 0
        for i in range(len(coverage)):
            if coverage[i] > 1e-3:
                error += (coverage[i] / answer - 1.)**2
                size += 1

        print 'error',numerix.sqrt(error / size)


        

    raw_input('finished')
