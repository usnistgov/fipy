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
This example first imposes a circular distance function:

.. raw:: latex

    $$ \\phi \\left( x, y \\right) = \\left[ \\left( x - \\frac{ L }{ 2 } \\right)^2 + \\left( y - \\frac{ L }{ 2 } \\right)^2 \\right]^{1/2} - \\frac{L}{4} $$ 

then the variable is advected with,

.. raw:: latex

    $$ \\frac{ \\partial \\phi } { \\partial t } + \\vec{u} \\cdot \\nabla \\phi = 0 $$

Also a surfactant is present of the interface, governed by the equation:

.. raw:: latex

    $$ \\frac{d \\theta}{d t} = J v \\theta $$

The result can be tested with the following code:


   >>> for step in range(steps):
   ...     it.timestep(dt = timeStepDuration)
   >>> ids = Numeric.nonzero(Numeric.logical_and(distanceVariable > 0., distanceVariable < dx / 2.))
   >>> surfactantValues = Numeric.take(surfactantVariable, ids)
   >>> cellCenters = Numeric.take(mesh.getCellCenters(), ids)
   >>> finalRadius = Numeric.sqrt((cellCenters[:,0]- L / 2)**2 + (cellCenters[:,1] - L / 2)**2)
   >>> answer =  initialRadius / finalRadius
   >>> Numeric.allclose(answer, Numeric.array(surfactantValues), rtol = 0.03)
   1

"""

import Numeric
   
from fipy.meshes.grid2D import Grid2D
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.variables.cellVariable import CellVariable
from fipy.models.levelSet.distanceFunction.distanceFunctionEquation import DistanceFunctionEquation
from fipy.models.levelSet.advection.advectionEquation import AdvectionEquation
from fipy.models.levelSet.advection.higherOrderAdvectionTerm import HigherOrderAdvectionTerm
from fipy.models.levelSet.surfactant.surfactantEquation import SurfactantEquation
from fipy.iterators.iterator import Iterator
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.solvers.linearLUSolver import LinearLUSolver
from fipy.boundaryConditions.fixedValue import FixedValue


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

distanceVariable = CellVariable(
    name = 'level set variable',
    mesh = mesh,
    value = 1.
    )

cellRadius = Numeric.sqrt((mesh.getCellCenters()[:,0] - L / 2.)**2 + (mesh.getCellCenters()[:,1] - L / 2.)**2)


initialSurfactantValue =  Numeric.where(initialRadius / cellRadius < 2., initialRadius / cellRadius, 2.)

surfactantVariable = CellVariable(
    name = 'surfactant variable',
    mesh = mesh,
    value = initialSurfactantValue
    )

distanceVariable.setValue(cellRadius - initialRadius)

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

it = Iterator((surfactantEquation, advectionEquation))

if __name__ == '__main__':
    
    distanceViewer = Grid2DGistViewer(var = distanceVariable, palette = 'rainbow.gp', minVal = -initialRadius, maxVal = initialRadius)
    surfactantViewer = Grid2DGistViewer(var = surfactantVariable, palette = 'rainbow.gp', minVal = -1., maxVal = 1.)
    distanceViewer.plot()
    surfactantViewer.plot()
    
    for step in range(steps):
        it.timestep(dt = timeStepDuration)
##        array = Numeric.array(distanceVariable)
##        array = array - timeStepDuration * velocity
##        distanceVariable.setValue(array)
        distanceViewer.plot()
        surfactantViewer.plot()
    print 'surfactantVariable',surfactantVariable
    raw_input('finished')
