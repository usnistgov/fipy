#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "main.py"
 #                                    created: 8/20/04 {10:29:10 AM} 
 #                                last update: 9/3/04 {10:43:20 PM} { 1:23:41 PM}
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

This file build the distance variable and equation.

"""

__docformat__ = 'restructuredtext'

import Numeric

## generate the mesh

xCells = controlParameters["cells below trench"] + int((geometryParameters["depth of trench"] + geometryParameters["boundary layer depth"]) / controlParameters["cell size"])
yCells = int((geometryParameters["trench spacing"] / 2) / controlParameters["cell size"])

mesh = Grid2D(dx = cellSize,
              dy = cellSize,
              nx = xCells,
              ny = yCells)

   
## generate the distance variable and equation

bottomHeight = controlParameters["number of cells below trench"] * cellSize
trenchHeight = bottomDepth + geometryParameters["trench depth"]
trenchWidth = trenchDepth / geometryParameters["aspect ratio"]
sideWidth = (geometryParameters["trench spacing"] - trenchWidth) / 2

distanceVariable = DistanceVariable(
    name = 'level set variable',
    mesh = mesh,
    value = -1.
    )

def electrolyteID(cell):
    x,y = cell.getCenter()
    if y > trenchHeight:
        return 1
    else if y < bottomHeight:
        return 0
    else if x < sideWidth:
        return 0:
    else:
        return 1
    
electrolyteCells = mesh.getCells(electrolyteID)

distanceVariable.setValue(1., electrolyteCells)

terminationValue = controlParameters["number of cells in narrow band"] / 2 * cellSize

distanceEquation = DistanceEquation(distanceVariable, terminationValue = terminationValue)

distanceEquation.solve()

## generate the accelerator variable

acceleratorVariable = SurfactantVariable(
    name = "accelerator",
    value = parameters.experimentalParameters["initial accelerator coverage"],
    distanceVariable = distanceVariable
    )

## generate the cupric variable

cupricVariable = CellVariable(
    name = 'cupric concentration',
    mesh = mesh,
    value = parameters.experimentalParameters["bulk cupric concentration"]
    )
    
## generate the deposition rate variable

depsoitionRateVariable = DepositionRateVariable(
    metalVariable = cupricVariable,
    acceleratorVariable = acceleratorVariable,
    experimentalParameter = parameters.experimentalParameters,
    materialParameters = materialParameters,
    metalParameters = metalParameters
    )

## generate the equations

surfactantEquation = SurfactantEquation(
    surfactantVariable,
    distanceVariable,
    solver = LinearLUSolver(
        tolerance = 1e-10),
    boundaryConditions = (FixedValue(mesh.getExteriorFaces(), 0),))

advectionEquation = HigherOrderAdvectionEquation(
    distanceVariable,
    advectionCoeff = depositionRateVariable

extensionEquation = ExtensionEquation(
    distanceVariable,
    velocity,
    terminationValue = terminationValue)

cupricEquation = MetalIonDiffusionEquation(
    cupricVariable,
    distanceVariable = distanceVariable,
    depositionRate = depositionRate,
    diffusionCoeff = parametes.metalParameters['diffusion coefficient'],
    metalAtomicVolume = parameters.metalParameters['atomic volume'],
    boundaryConditions = (
        FixedValue(parameters.experimentalParameters["bulk cupric concentration"],
                   mesh.getFacesTop(),
                   )
        )
    
it = Iterator((extensionEquation, advectionEquation, surfactantEquation, cupricEquation))

distanceViewer = Grid2DGistViewer(var = distanceVariable,
                                      minVal = -.001,
                                      maxVal = .001)

surfactantViewer = Grid2DGistViewer(var = surfactantVariable,
                                    minVal = 0.,
                                    maxVal = 200.)

cupricViewer = Grid2DGistViewer(var = curicVariable,
                                minVal = 0.,
                                maxVal = parameters.experimentalParameters["bulk cupric concentration"])

viewers = (
    distanceViewer = Grid2DGistViewer(var = distanceVariable,
                                      minVal = -.001,
                                      maxVal = .001),
    surfactantViewer = Grid2DGistViewer(var = surfactantVariable,
                                    minVal = 0.,
                                    maxVal = 200.),
    cupricViewer = Grid2DGistViewer(var = curicVariable,
                                minVal = 0.,
                                maxVal = parameters.experimentalParameters["bulk cupric concentration"])
    )

if __name__ == '__main__':

    for viewer in viewers:
        viewer.plot()

    raw_input("press key to continue")

##    from fipy.tools.profiler.profiler import calibrate_profiler
##    from fipy.tools.profiler.profiler import Profiler
##    fudge = calibrate_profiler(10000)
##    profile = Profiler('profile.txt', fudge=fudge)

    totalTime = 0.
    argmax = Numeric.argmax(depositionRate)
    timeStepDuration = parameters.controlParameters["cfl number"] * cellSize / depeositionRate[argmax]
    levelSetUpdateFrequency = int(depeositionRate[argmax] * timeStepDuration / terminationValue * 4 / 5)

    while totalTime < parameters.controlParameters["simulation duration"]

        if step % levelSetUpdateFrequency == 0:
            distanceEquation.solve()

        argmax = Numeric.argmax(depositionRate)
        timeStepDuration = cfl * dx / depositionRate[argmax]
        totalTime += timeStepDuration

        it.timestep(dt = timeStepDuration)
        
##      profile.stop()

        for viewer in viewers:
            viewer.plot()
        
    raw_input('finished')
