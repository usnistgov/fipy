#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "adaptiveMeshInput.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 10/6/04 {4:48:38 PM} { 1:23:41 PM}
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

    $$ \\phi \\left( x, y \\right) = \\left[ \\left( x - \\frac{ L }{ 2 } \\right)^2 
         + \\left( y - \\frac{ L }{ 2 } \\right)^2 \\right]^{1/2} - \\frac{L}{4} $$ 

then the variable is advected with,

.. raw:: latex

    $$ \\frac{ \\partial \\phi } { \\partial t } + \\vec{u} \\cdot \\nabla \\phi = 0 $$

The scheme used in the `AdvectionTerm` preserves the `distanceVariable` as a distance function.

This example demonstrates the use of the AdaptiveMesh2D and
ReMeshedCellVariable class to dynamically modify the mesh based on the
variable value.  This is done as follows:

First, have the iterator time step

Then, create an \"adaptation variable\" with a mesh equal to the previous
mesh and nodal values equal to the target sizes of the mesh elements at the
corresponding points.  In this example, we would like the mesh to be finer
(lower characteristic lengths) where phi is close to zero so we use a
characterisitc length equal to (3(phi ** 2)) + 0.03

Then, create a new mesh by passing the adaptVar to AdaptiveMesh2D

Then, we need to tell the variable that its mesh is being changed.  To do
this, create a ReMeshedCellVariable, which takes as its arguments an old
variable and a new mesh and returns a vvariable that has the old variable's
values mapped onto the new mesh

Then, re-initialize the equation with the new variable and the iterator
with the new equation

   >>> for step in range(steps):
   ...     it.timestep(dt = timeStepDuration)
   ...     adaptVar = CellVariable(
   ...         name = \"characteristic lengths\",
   ...         value = (3*(Numeric.array(distanceVariable) * Numeric.array(distanceVariable))) + 0.03,
   ...         mesh = mesh)
   ...     newMesh = AdaptiveMesh2D(adaptVar)
   ...     newVariable = ReMeshedCellVariable(distanceVariable, newMesh)
   ...     distanceVariable = newVariable
   ...     mesh = newMesh
   ...     advectionEquation = AdvectionEquation(
   ...         distanceVariable,
   ...         advectionCoeff = velocity,
   ...         solver = LinearPCGSolver(
   ...         tolerance = 1.e-15, 
   ...         steps = 1000))
   ...     it = Iterator((advectionEquation,))
   
We now determine the accuracy of the results.  To do so, we determine
whether the values that are closest to zero really do occur on the circle
that they are supposed to.  For this test, all of the 10 values that are
closest to zero must occur near the designated circle.

   >>> distancesFromCenter = Numeric.sqrt((mesh.getCellCenters()[:,0] - L / 2.)**2 + (mesh.getCellCenters()[:,1] - L / 2.)**2)
   >>> targetDistance = (timeStepDuration * steps * velocity) + radius
   >>> bestCells = Numeric.argsort((Numeric.array(distanceVariable) * Numeric.array(distanceVariable)))[:10]
   >>> bestDistances = Numeric.take(distancesFromCenter, bestCells)
   >>> Numeric.allclose(bestDistances, targetDistance, atol = 0.031)
   1

"""
__docformat__ = 'restructuredtext'

import Numeric
import os   
from fipy.meshes.numMesh.tri2D import Tri2D
from fipy.viewers.pyxviewer import PyxViewer
from fipy.variables.cellVariable import CellVariable, ReMeshedCellVariable
from fipy.models.levelSet.distanceFunction.distanceEquation import DistanceEquation
from fipy.models.levelSet.advection.advectionEquation import AdvectionEquation
from fipy.iterators.iterator import Iterator
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.models.levelSet.advection.higherOrderAdvectionTerm import HigherOrderAdvectionTerm
## from fipy.tools.profiler.profiler import Profiler
## from fipy.tools.profiler.profiler import calibrate_profiler
from fipy.meshes.numMesh.adaptiveMesh import AdaptiveMesh2D
from fipy.meshes.numMesh.gmshExport import exportAsMesh

##fudge = calibrate_profiler(10000)
##profile = Profiler('profile', fudge=fudge)

L = 1.
nx = 25
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
    print steps
    for step in range(steps):
        print step, "steps completed"
        it.timestep(dt = timeStepDuration)
        adaptVar = CellVariable(
            name = "characteristic lengths",
            value = (3 * (Numeric.array(distanceVariable) * Numeric.array(distanceVariable))) + 0.03,
            mesh = mesh)
        newMesh = AdaptiveMesh2D(adaptVar)
        newVariable = ReMeshedCellVariable(distanceVariable, newMesh)
        distanceVariable = newVariable
        mesh = newMesh
        advectionEquation = AdvectionEquation(
            distanceVariable,
            advectionCoeff = velocity,
            solver = LinearPCGSolver(
            tolerance = 1.e-15, 
            steps = 1000))
        it = Iterator((advectionEquation,))
        if(step % 5 == 4):
            viewer = PyxViewer(distanceVariable)
            viewer.plot()

##profile.stop()
