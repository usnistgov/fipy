#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 4/7/05 {4:50:24 PM} { 1:23:41 PM}
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

The result can be tested with the following code:


   >>> for step in range(steps):
   ...     it.timestep(dt = timeStepDuration)
   
   >>> x = Numeric.array(mesh.getCellCenters())
   >>> distanceTravelled = timeStepDuration * steps * velocity
   >>> answer = initialArray - distanceTravelled
   >>> answer = Numeric.where(answer < 0., -1001., answer)
   >>> solution = Numeric.where(answer < 0., -1001., Numeric.array(distanceVariable))
   >>> Numeric.allclose(answer, solution, atol = 2.5e-3)
   1
   
"""
__docformat__ = 'restructuredtext'

import Numeric
   
from fipy.meshes.tri2D import Tri2D
import fipy.viewers
from fipy.variables.cellVariable import CellVariable
from fipy.models.levelSet.distanceFunction.distanceEquation import DistanceEquation
from fipy.models.levelSet.advection.advectionEquation import AdvectionEquation
from fipy.iterators.iterator import Iterator
from fipy.solvers.linearPCGSolver import LinearPCGSolver


L = 1.
nx = 50
velocity = 1.
cfl = 0.1
velocity = 1.
distanceToTravel = L / 10.
radius = L / 4.

dx = L / nx
timeStepDuration = cfl * dx / velocity
steps = int(distanceToTravel / dx / cfl)

mesh = Tri2D(dx = dx, dy = dx, nx = nx, ny = nx)

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
    distanceViewer = fipy.viewers.make(vars = distanceVariable)
    distanceViewer.plot()

    for step in range(steps):
        
        it.timestep(dt = timeStepDuration)
        distanceViewer.plot()

    raw_input('finished')
