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
from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
from fipy.models.levelSet.distanceFunction.distanceEquation import DistanceEquation
from fipy.models.levelSet.advection.advectionEquation import AdvectionEquation
from fipy.models.levelSet.surfactant.surfactantEquation import SurfactantEquation
from fipy.models.levelSet.advection.higherOrderAdvectionTerm import HigherOrderAdvectionTerm
from fipy.iterators.iterator import Iterator
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.solvers.linearLUSolver import LinearLUSolver
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.models.levelSet.surfactant.surfactantVariable import SurfactantVariable
from fipy.variables.cellVariable import CellVariable
from fipy.models.levelSet.distanceFunction.extensionEquation import ExtensionEquation

Lx = .5
Ly = 1.

dx = 0.005
cfl = 0.5

nx = int(Lx / dx)
ny = int(Ly / dx)
steps = 5000

mesh = Grid2D(dx = dx, dy = dx, nx = nx, ny = ny)

distanceVariable = DistanceVariable(
    name = 'level set variable',
    mesh = mesh,
    value = -1.
    )

positiveCells = mesh.getCells(lambda cell: (cell.getCenter()[1] > 3 * Ly / 4) or
                              (cell.getCenter()[0] > Lx / 2 and cell.getCenter()[1] > Ly / 4) )

distanceVariable.setValue(1., positiveCells)

distanceEquation = DistanceEquation(distanceVariable, terminationValue = 6 * dx)

distanceEquation.solve()

surfactantVariable = SurfactantVariable(
    value = 1,
    distanceVariable = distanceVariable
    )

velocity = CellVariable(
    name = 'velocity',
    mesh = mesh,
    value = 1.,
    )

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

extensionEquation = ExtensionEquation(
    distanceVariable,
    velocity,
    terminationValue = 6 * dx)

it = Iterator((extensionEquation, advectionEquation, surfactantEquation))

faces = mesh.getFacesRight()
faceIDs = [face.getID() for face in faces]
cellIDs = Numeric.take(mesh.getFaceCellIDs()[:,0], faceIDs)

def calcCenterHeight(phi):
    for cellID in cellIDs:
        if phi[cellID] < 0 and phi[cellID + nx] > 0:
            return float(mesh.getCellCenters()[cellID, 1] - phi[cellID])


if __name__ == '__main__':

    distanceViewer = Grid2DGistViewer(var = distanceVariable, palette = 'rainbow.gp', minVal = -.001, maxVal = .001)
    levelSetViewer = Grid2DGistViewer(var = distanceVariable, palette = 'rainbow.gp', minVal = -.5, maxVal = .5)
    surfactantViewer = Grid2DGistViewer(var = surfactantVariable, palette = 'rainbow.gp', minVal = 0., maxVal = 200.)
    velocityViewer = Grid2DGistViewer(var = velocity, palette = 'rainbow.gp', minVal = 0., maxVal = 2.)
    levelSetGradMag = Grid2DGistViewer(var = distanceVariable.getGrad().getMag(), palette = 'rainbow.gp', minVal = .5, maxVal = 1.5)

    distanceViewer.plot()
    surfactantViewer.plot()
    velocityViewer.plot()


##    from fipy.tools.profiler.profiler import calibrate_profiler
##    from fipy.tools.profiler.profiler import Profiler
##    fudge = calibrate_profiler(10000)
##    profile = Profiler('profile.txt', fudge=fudge)
##    from fipy.viewers.gist1DViewer import Gist1DViewer
##    totalTime = 0

##    import fipy.tools.dump as dump

##    stringCenterHeights = dump.read('centerHeights')
    
##    centerHeights = [[0,0]]
##    centerView = Gist1DViewer(vars = (centerHeights, stringCenterHeights), title = 'Center Height')
    totalTime = 0
    argmax = Numeric.argmax(velocity)
    timeStepDuration = cfl * dx / velocity[argmax]
    for step in range(steps):
        print 'step',step

        if step % 10 == 0:
            distanceEquation.solve()

        velocity.setValue(surfactantVariable.getInterfaceValue())

        argmax = Numeric.argmax(velocity)
        timeStepDuration = cfl * dx / velocity[argmax]
        totalTime += timeStepDuration
        it.timestep(dt = timeStepDuration)
        
##        import fipy.tools.memoryLeak as memoryLeak
##        memoryLeak.print_top_100()
##        raw_input("press key to continue")
##    profile.stop()

        distanceViewer.plot()
        surfactantViewer.plot()
        velocityViewer.plot()
        levelSetViewer.plot()
        
##        centerHeights.append([totalTime, calcCenterHeight(distanceVariable) - 0.25])

##        centerView.plot()


    
        
    raw_input('finished')
