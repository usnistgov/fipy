#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 7/13/05 {3:31:33 PM} { 1:23:41 PM}
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
This example first imposes a circular distance function:

.. raw:: latex

    $$ \phi \left( x, y \right) = \left[ \left( x - \frac{ L }{ 2 } \right)^2 + \left( y - \frac{ L }{ 2 } \right)^2 \right]^{1/2} - \frac{L}{4} $$ 

The variable is advected with,

.. raw:: latex

    $$ \frac{ \partial \phi } { \partial t } + \vec{u} \cdot \nabla \phi = 0 $$

The scheme used in the `_AdvectionTerm` preserves the `var` as a
distance function.  The solution to this problem will be demonstrated
in the following script. Firstly, setup the parameters.

   >>> L = 1.
   >>> nx = 25
   >>> velocity = 1.
   >>> cfl = 0.1
   >>> velocity = 1.
   >>> distanceToTravel = L / 10.
   >>> radius = L / 4.
   >>> dx = L / nx   
   >>> timeStepDuration = cfl * dx / velocity
   >>> steps = int(distanceToTravel / dx / cfl)

Construct the mesh.

   >>> from fipy.meshes.grid2D import Grid2D
   >>> mesh = Grid2D(dx = dx, dy = dx, nx = nx, ny = nx)

Construct a `distanceVariable` object.

   >>> from fipy.models.levelSet.distanceFunction.distanceVariable \
   ...     import DistanceVariable
   >>> var = DistanceVariable(
   ...     name = 'level set variable',
   ...     mesh = mesh,
   ...     value = 1.)

Initialise the `distanceVariable` to be a circular distance function.

   >>> from fipy.tools import numerix
   >>> initialArray = numerix.sqrt((mesh.getCellCenters()[:,0] - L / 2.)**2 +
   ...                             (mesh.getCellCenters()[:,1] - L / 2.)**2) - \
   ...                             radius
   >>> var.setValue(initialArray)

The `advectionEquation` is constructed.
   
   >>> from fipy.models.levelSet.advection.advectionEquation import \
   ...     buildAdvectionEquation
   >>> advEqn = buildAdvectionEquation(
   ...     advectionCoeff = velocity)

The problem can then be solved by executing a serious of time steps.

   >>> if __name__ == '__main__':
   ...     import fipy.viewers
   ...     viewer = fipy.viewers.make(vars = var, 
   ...         limits = {'datamin': -radius, 'datamax': radius})
   ...     viewer.plot()
   ...     for step in range(steps):
   ...         var.updateOld()
   ...         advEqn.solve(var, dt = timeStepDuration)
   ...         viewer.plot()

The result can be tested with the following commands.

   >>> for step in range(steps):
   ...     var.updateOld()
   ...     advEqn.solve(var, dt = timeStepDuration)
   >>> x = numerix.array(mesh.getCellCenters())
   >>> distanceTravelled = timeStepDuration * steps * velocity
   >>> answer = initialArray - distanceTravelled
   >>> answer = numerix.where(answer < 0., -1001., answer)
   >>> solution = numerix.where(answer < 0., -1001., numerix.array(var))
   >>> numerix.allclose(answer, solution, atol = 4.7e-3)
   1

If the `AdvectionEquation` is built with the `_HigherOrderAdvectionTerm` the result
is more accurate,

   >>> var.setValue(initialArray)
   >>> from fipy.models.levelSet.advection.higherOrderAdvectionEquation \
   ...     import buildHigherOrderAdvectionEquation
   >>> advEqn = buildHigherOrderAdvectionEquation(
   ...     advectionCoeff = velocity)
   >>> for step in range(steps):
   ...     var.updateOld()
   ...     advEqn.solve(var, dt = timeStepDuration)
   >>> solution = numerix.where(answer < 0., -1001., numerix.array(var))
   >>> numerix.allclose(answer, solution, atol = 1.02e-3)
   1

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    raw_input("finished")
