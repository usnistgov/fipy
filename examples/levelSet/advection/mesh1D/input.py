#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 3/8/05 {3:49:57 PM} { 1:23:41 PM}
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

r"""
This example first solves the distance function equation in one dimension:

.. raw:: latex

    $$ |\nabla \phi| = 1 $$

with

.. raw:: latex

    $\phi = 0$ at $x = L / 5$.

The variable is then advected with,

.. raw:: latex

    $$ \frac{ \partial \phi } { \partial t} + \vec{u} \cdot \nabla \phi = 0 $$

The scheme used in the `AdvectionTerm` preserves the `var` as a distance function.

The solution to this problem will be demonstrated in the following
script. Firstly, setup the parameters.

   >>> velocity = 1.
   >>> dx = 1.
   >>> nx = 10
   >>> timeStepDuration = 1.
   >>> steps = 2
   >>> L = nx * dx
   >>> interfacePosition = L / 5.

Construct the mesh.

   >>> from fipy.meshes.grid2D import Grid2D
   >>> mesh = Grid2D(dx = dx, nx = nx)

Construct a `distanceVariable` object. This object is required by the
`distanceEquation`.

   >>> import Numeric
   >>> values = -Numeric.ones(nx, 'd')
   >>> cells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] > interfacePosition)
   >>> for cell in cells:
   ...    values[cell.getID()] = 1

   >>> from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
   >>> var = DistanceVariable(name = 'level set variable',
   ...                        mesh = mesh,
   ...                        value = values)
   >>> var.calcDistanceFunction()
   
The `advectionEquation` is constructed.

   >>> from fipy.models.levelSet.advection.advectionEquation import buildAdvectionEquation
   >>> advEqn = buildAdvectionEquation(advectionCoeff = velocity)

The problem can then be solved by executing a serious of time steps.

   >>> if __name__ == '__main__':
   ...     import fipy.viewers
   ...     viewer = fipy.viewer.make(vars = var, limits = {'datamin': -10., 'datamax': 10.})
   ...     viewer.plot()
   ...     for step in range(steps):
   ...         var.updateOld()
   ...         advEqn.solve(var, dt = timeStepDuration)
   ...         viewer.plot()

The result can be tested with the following code:

   >>> for step in range(steps):
   ...     var.updateOld()
   ...     advEqn.solve(var, dt = timeStepDuration)
   >>> import Numeric
   >>> x = Numeric.array(mesh.getCellCenters()[:,0])
   >>> distanceTravelled = timeStepDuration * steps * velocity
   >>> answer = x - interfacePosition - timeStepDuration * steps * velocity
   >>> answer = Numeric.where(x < distanceTravelled, x[0] - interfacePosition, answer)
   >>> Numeric.allclose(answer, var, atol = 1e-10)
   1
   
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())
    raw_input("finished")
