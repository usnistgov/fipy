#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 10/6/04 {5:30:55 PM} { 1:23:41 PM}
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
This example first imposes a circular distance function:

.. raw:: latex

    $$ \\phi \\left( x, y \\right) = \\left[ \\left( x - \\frac{ L }{ 2 } \\right)^2 + \\left( y - \\frac{ L }{ 2 } \\right)^2 \\right]^{1/2} - \\frac{L}{4} $$ 

then the variable is advected with,

.. raw:: latex

    $$ \\frac{ \\partial \\phi } { \\partial t } + \\vec{u} \\cdot \\nabla \\phi = 0 $$

The scheme used in the `AdvectionTerm` preserves the `distanceVariable` as a distance function.

This is just like the regular advection circle Grid2D test case except that it is tested the same way as the adaptiveMeshInput test case is tested. This allows for comparison of results between unstructured (adaptive) meshes and structured meshes.

   >>> for step in range(steps):
   ...     it.timestep(dt = timeStepDuration)
   
   >>> distancesFromCenter = Numeric.sqrt((mesh.getCellCenters()[:,0] - L / 2.)**2 + (mesh.getCellCenters()[:,1] - L / 2.)**2)
   >>> targetDistance = (timeStepDuration * steps * velocity) + radius
   >>> bestCells = Numeric.argsort((Numeric.array(distanceVariable) * Numeric.array(distanceVariable)))[:10]
   >>> bestDistances = Numeric.take(distancesFromCenter, bestCells)
   >>> Numeric.allclose(bestDistances, targetDistance, atol = 0.01)
   1

"""
__docformat__ = 'restructuredtext'

import Numeric
   
from fipy.meshes.grid2D import Grid2D
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.variables.cellVariable import CellVariable
from fipy.models.levelSet.distanceFunction.distanceEquation import DistanceEquation
from fipy.models.levelSet.advection.advectionEquation import AdvectionEquation
from fipy.iterators.iterator import Iterator
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.models.levelSet.advection.higherOrderAdvectionTerm import HigherOrderAdvectionTerm
## from fipy.tools.profiler.profiler import Profiler
## from fipy.tools.profiler.profiler import calibrate_profiler

## fudge = calibrate_profiler(10000)
## profile = Profiler('profile', fudge=fudge)

L = 1.
nx = 100
velocity = 1.
cfl = 0.1
velocity = 1.
distanceToTravel = L / 10.
radius = L / 4.

dx = L / nx
timeStepDuration = cfl * dx / velocity
steps = int(distanceToTravel / dx / cfl)

mesh = Grid2D(dx = dx, dy = dx, nx = nx, ny = nx)

distanceVariable = CellVariable(
    name = 'level set variable',
    mesh = mesh,
    value = 1.
    )

initialArray = Numeric.sqrt((mesh.getCellCenters()[:,0] - L / 2.)**2 + (mesh.getCellCenters()[:,1] - L / 2.)**2) - radius

distanceVariable.setValue(initialArray)

advectionEquation = AdvectionEquation(
    distanceVariable,
    advectionCoeff = velocity,
    solver = LinearPCGSolver(
        tolerance = 1.e-15, 
        steps = 1000))

it = Iterator((advectionEquation,))

if __name__ == '__main__':
    distanceViewer = Grid2DGistViewer(var = distanceVariable, palette = 'rainbow.gp', minVal = -radius, maxVal = radius)
    distanceViewer.plot()
    print steps
    for step in range(steps):
        
        it.timestep(dt = timeStepDuration)
        distanceViewer.plot()
## profile.stop()
