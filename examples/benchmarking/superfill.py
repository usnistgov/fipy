#!/usr/bin/env python

## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "superfill.py"
 #                                     created: 1/19/06 {4:09:41 PM}
 #                                 last update: 2/2/07 {8:47:12 AM}
 # Author: Jonathan Guyer
 # E-mail: <guyer@nist.gov>
 # Author: Daniel Wheeler
 # E-mail: daniel.wheeler@nist.gov
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
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
 # History
 # 
 # modified   by  rev reason
 # ---------- --- --- -----------
 # 1940-01-18 JEG 1.0 original
 # 
 # ########################################################################
 ##

r"""

This example benchmarks the speed and memory usage of solving an
electrochemistry superfill problem. Run:

    $ python setup.py efficiency_test

"""
__docformat__ = 'restructuredtext'

if __name__ == "__main__":
    
    faradaysConstant = 9.6e4
    gasConstant = 8.314
    transferCoefficient = 0.5
    rateConstant0 = 1.76
    rateConstant3 = -245e-6
    catalystDiffusion = 1e-9
    siteDensity = 9.8e-6
    molarVolume = 7.1e-6,
    charge = 2
    metalDiffusionCoefficient = 5.6e-10
    temperature = 298.
    overpotential = -0.3
    bulkMetalConcentration = 250.
    catalystConcentration = 5e-3
    catalystCoverage = 0.
    currentDensity0 = 0.26
    currentDensity1 = 45.
    cflNumber = 0.2
    numberOfCellsInNarrowBand = 10
    cellsBelowTrench = 10
    cellSize = 0.1e-7
    trenchDepth = 0.5e-6
    aspectRatio = 2.
    trenchSpacing = 0.6e-6
    boundaryLayerDepth = 0.3e-6

    from fipy.tools.parser import parse

    numberOfElements = parse('--numberOfElements', action = 'store',
        type = 'int', default = -1)

    from benchmarker import Benchmarker
    bench = Benchmarker()

    numberOfSteps = 5

    bench.start()

    import fipy.tools.numerix as numerix
    if numberOfElements != -1:
        pos = trenchSpacing * cellsBelowTrench / 4 / numberOfElements
        sqr = trenchSpacing * (trenchDepth + boundaryLayerDepth) \
              / (2 * numberOfElements)
        cellSize = pos + numerix.sqrt(pos**2 + sqr)
    else:
        cellSize = 0.1e-7

    yCells = cellsBelowTrench \
             + int((trenchDepth + boundaryLayerDepth) / cellSize)

    xCells = int(trenchSpacing / 2 / cellSize)
    from fipy.meshes.grid2D import Grid2D
    mesh = Grid2D(dx = cellSize,
                  dy = cellSize,
                  nx = xCells,
                  ny = yCells)

    bench.stop('mesh')

    bench.start()

    narrowBandWidth = numberOfCellsInNarrowBand * cellSize
    from fipy.models.levelSet.distanceFunction.distanceVariable import \
        DistanceVariable        

    distanceVar = DistanceVariable(
       name = 'distance variable',
       mesh = mesh,
       value = -1,
       narrowBandWidth = narrowBandWidth,
       hasOld = 1)

    bottomHeight = cellsBelowTrench * cellSize
    trenchHeight = bottomHeight + trenchDepth
    trenchWidth = trenchDepth / aspectRatio
    sideWidth = (trenchSpacing - trenchWidth) / 2
    x, y = mesh.getCellCenters()[...,0], mesh.getCellCenters()[...,1]
    distanceVar.setValue(1, where=(y > trenchHeight) 
                                  | ((y > bottomHeight) 
                                     & (x < xCells * cellSize - sideWidth)))

    distanceVar.calcDistanceFunction(narrowBandWidth = 1e10)
    from fipy.models.levelSet.surfactant.surfactantVariable import \
        SurfactantVariable

    catalystVar = SurfactantVariable(
        name = "catalyst variable",
        value = catalystCoverage,
        distanceVar = distanceVar)

    from fipy.variables.cellVariable import CellVariable
    bulkCatalystVar = CellVariable(
        name = 'bulk catalyst variable',
        mesh = mesh,
        value = catalystConcentration)

    from fipy.variables.cellVariable import CellVariable
    metalVar = CellVariable(
        name = 'metal variable',
        mesh = mesh,
        value = bulkMetalConcentration)

    bench.stop('variables')

    bench.start()

    expoConstant = -transferCoefficient * faradaysConstant \
                   / (gasConstant * temperature)

    tmp = currentDensity1 \
          * catalystVar.getInterfaceVar()

    exchangeCurrentDensity = currentDensity0 + tmp
    expo = numerix.exp(expoConstant * overpotential)
    currentDensity = expo * exchangeCurrentDensity * metalVar \
                     / bulkMetalConcentration

    depositionRateVariable = currentDensity * molarVolume \
                             / (charge * faradaysConstant)

    extensionVelocityVariable = CellVariable(
        name = 'extension velocity',
        mesh = mesh,
        value = depositionRateVariable)   

    from fipy.models.levelSet.surfactant.adsorbingSurfactantEquation \
                import AdsorbingSurfactantEquation

    surfactantEquation = AdsorbingSurfactantEquation(
        surfactantVar = catalystVar,
        distanceVar = distanceVar,
        bulkVar = bulkCatalystVar,
        rateConstant = rateConstant0 \
                       + rateConstant3 * overpotential**3)

    from fipy.models.levelSet.advection.higherOrderAdvectionEquation \
                   import buildHigherOrderAdvectionEquation

    advectionEquation = buildHigherOrderAdvectionEquation(
        advectionCoeff = extensionVelocityVariable)

    from fipy.models.levelSet.electroChem.metalIonDiffusionEquation \
                         import buildMetalIonDiffusionEquation

    metalEquation = buildMetalIonDiffusionEquation(
        ionVar = metalVar,
        distanceVar = distanceVar,
        depositionRate = depositionRateVariable,
        diffusionCoeff = metalDiffusionCoefficient,
        metalIonMolarVolume = molarVolume,
    )

    from fipy.models.levelSet.surfactant.surfactantBulkDiffusionEquation \
                    import buildSurfactantBulkDiffusionEquation

    bulkCatalystEquation = buildSurfactantBulkDiffusionEquation(
        bulkVar = bulkCatalystVar,
        distanceVar = distanceVar,
        surfactantVar = catalystVar,
        diffusionCoeff = catalystDiffusion,
        rateConstant = rateConstant0 * siteDensity
    )

    bench.stop('terms')

    bench.start()

    from fipy.boundaryConditions.fixedValue import FixedValue
    metalEquationBCs = (
            FixedValue(
                mesh.getFacesTop(),
                bulkMetalConcentration
            ),
        )

    catalystBCs = (
            FixedValue(
                mesh.getFacesTop(),
                catalystConcentration
            ),)

    bench.stop('BCs')

    levelSetUpdateFrequency = int(0.8 * narrowBandWidth \
                                  / (cellSize * cflNumber * 2))

    step = 0

    if step % levelSetUpdateFrequency == 0:
        distanceVar.calcDistanceFunction()

    extensionVelocityVariable.setValue(depositionRateVariable())

    distanceVar.updateOld()
    catalystVar.updateOld()
    metalVar.updateOld()
    bulkCatalystVar.updateOld()
    distanceVar.extendVariable(extensionVelocityVariable)
    dt = cflNumber * cellSize / numerix.max(extensionVelocityVariable)
    advectionEquation.solve(distanceVar, dt = dt)
    surfactantEquation.solve(catalystVar, dt = dt)
    metalEquation.solve(metalVar, dt = dt,
                        boundaryConditions = metalEquationBCs)
    bulkCatalystEquation.solve(bulkCatalystVar, dt = dt,
                               boundaryConditions = catalystBCs)

    bench.start()
      
    for step in range(numberOfSteps):

        if step % levelSetUpdateFrequency == 0:
            distanceVar.calcDistanceFunction()

        extensionVelocityVariable.setValue(depositionRateVariable())

        distanceVar.updateOld()
        catalystVar.updateOld()
        metalVar.updateOld()
        bulkCatalystVar.updateOld()
        distanceVar.extendVariable(extensionVelocityVariable)
        dt = cflNumber * cellSize / numerix.max(extensionVelocityVariable)
        advectionEquation.solve(distanceVar, dt = dt)
        surfactantEquation.solve(catalystVar, dt = dt)
        metalEquation.solve(metalVar, dt = dt,
                            boundaryConditions = metalEquationBCs)
        bulkCatalystEquation.solve(bulkCatalystVar, dt = dt,
                                      boundaryConditions = catalystBCs)

    bench.stop('solve')

    print bench.report(numberOfElements=numberOfElements, steps=numberOfSteps)
