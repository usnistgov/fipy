#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "inputLeveler.py"
 #                                    created: 8/26/04 {10:29:10 AM} 
 #                                last update: 9/15/05 {7:03:58 PM} { 1:23:41 PM}
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

r"""
This input file

.. raw:: latex

    \label{inputLeveler} is a demonstration of the use of \FiPy{} for
    modeling copper superfill with leveler and accelerator
    additives. The material properties and experimental parameters
    used are roughly those that have been previously
    published~\cite{NIST:leveler:2005}.

To run this example from the base fipy directory type::
    
    $ examples/levelSet/electroChem/inputLeveler.py

at the command line. The results of the simulation will be displayed
and the word `finished` in the terminal at the end of the
simulation. The simulation will only run for 200 time steps. In order
to alter the number of timesteps, the python function that
encapsulates the system of equations must first be imported (at the
python command line),

    >>> from examples.levelSet.electroChem.inputLeveler import runLeveler

and then the function can be run with a different number of time steps
with the `numberOfSteps` argument as follows,

    >>> runLeveler(numberOfSteps = 10, displayViewers = False, cellSize = 0.25e-7)
    1

Change the `displayViewers` argument to `True` if you wish to see the
results displayed on the screen. This example requires `gmsh` to
construct the mesh.

This

.. raw:: latex

    example is for the case when there are suppressor, accelerator and
    leveler additives in the electrolyte. The suppressor is assumed to
    absorb quickly compared with the other additives. Any unoccupied
    surface sites are immediately covered with suppressor. The
    accelerator additive has more surface affinity than suppressor and
    is thus preferential adsorbed. The accelerator can also remove
    suppressor when the surface reaches full coverage. Similarly, the
    leveler additive has more surface affinity than both the
    suppressor and accelerator. This forms a simple set of assumptions
    for understanding the behavior of these additives.

    The following is a complete description of the equations for the
    model described here. Any equations that have been omitted are the
    same as those given in Example~\ref{inputSimpleTrench}. The
    current density is governed by $$ i = \frac{ c_m }{ c_m^\infty }
    \sum_{ j } \left[ i_j \theta_j \left( \exp{ \frac{\alpha_j F \eta
    }{ R T }} + \exp{ \frac{-\alpha_j F \eta}{ R T }} \right) \right]
    $$ where $j$ represents $s$ for suppressor, $a$ for accelerator
    and $l$ for leveler. This model assumes a linear interpolation
    between the three cases of complete coverage for each
    additive. The governing equations for the surfactants are given
    by, $$ \dot{\theta_{l}} = k_l c_l \left( 1 - \theta_l \right) -
    k_{lc} v \theta_l, $$ $$ \dot{\theta_{a}} = k_a c_a \left( 1 -
    \theta_a - \theta_l \right) - k_l c_l \theta_a - k_{ac}
    \theta_a^{q - 1} $$ and $$ \theta_s = 1 - \theta_a - \theta_l $$
    
"""
__docformat__ = 'restructuredtext'

def runLeveler(kLeveler = 0.018, bulkLevelerConcentration = 0.02, cellSize = 0.1e-7, rateConstant = 0.00026, initialAcceleratorCoverage = 0.0, levelerDiffusionCoefficient = 5e-10, numberOfSteps = 400, displayRate = 10, displayViewers = True):

    
    kLevelerConsumption = 0.0005
    aspectRatio = 1.5
    faradaysConstant = 9.6485e4
    gasConstant = 8.314
    acceleratorDiffusionCoefficient = 4e-10
    siteDensity = 6.35e-6
    atomicVolume = 7.1e-6
    charge = 2
    metalDiffusionCoefficient = 4e-10
    temperature = 298.
    overpotential = -0.25
    bulkMetalConcentration = 250.
    bulkAcceleratorConcentration = 50.0e-3
    initialLevelerCoverage = 0.
    cflNumber = 0.2
    numberOfCellsInNarrowBand = 20
    cellsBelowTrench = 10
    trenchDepth = 0.4e-6
    trenchSpacing = 0.6e-6
    boundaryLayerDepth = 98.7e-6
    i0Suppressor = 0.3
    i0Accelerator = 22.5
    alphaSuppressor = 0.5
    alphaAccelerator = 0.4
    alphaAdsorption = 0.62
    m = 4
    b = 2.65
    A = 0.3
    Ba = -40
    Bb = 60
    Vd = 0.098
    Bd = 0.0008

    etaPrime = faradaysConstant * overpotential / gasConstant / temperature

    from gapFillMesh import TrenchMesh
    from fipy.tools import numerix
    mesh = TrenchMesh(cellSize = cellSize,
                      trenchSpacing = trenchSpacing,
                      trenchDepth = trenchDepth,
                      boundaryLayerDepth = boundaryLayerDepth,
                      aspectRatio = aspectRatio,
                      angle = numerix.pi * 4. / 180.,
                      bowWidth = 0.,
                      overBumpRadius = 0.,
                      overBumpWidth = 0.)

    narrowBandWidth = numberOfCellsInNarrowBand * cellSize
    from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
    distanceVar = DistanceVariable(
        name = 'distance variable',
        mesh = mesh,
        value = -1,
        narrowBandWidth = narrowBandWidth)

    distanceVar.setValue(1, mesh.getElectrolyteCells())
    
    distanceVar.calcDistanceFunction(narrowBandWidth = 1e10)
    from fipy.models.levelSet.surfactant.surfactantVariable import SurfactantVariable
    levelerVar = SurfactantVariable(
        name = "leveler variable",
        value = initialLevelerCoverage,
        distanceVar = distanceVar)

    acceleratorVar = SurfactantVariable(
        name = "accelerator variable",
        value = initialAcceleratorCoverage,
        distanceVar = distanceVar)

    from fipy.variables.cellVariable import CellVariable
    bulkAcceleratorVar = CellVariable(name = 'bulk accelerator variable',
                                      mesh = mesh,
                                      value = bulkAcceleratorConcentration)

    bulkLevelerVar = CellVariable(
        name = 'bulk leveler variable',
        mesh = mesh,
        value = bulkLevelerConcentration)

    metalVar = CellVariable(
        name = 'metal variable',
        mesh = mesh,
        value = bulkMetalConcentration)

    def depositionCoeff(alpha, i0):
        expo = numerix.exp(-alpha * etaPrime)
        return 2 * i0 * (expo - expo * numerix.exp(etaPrime))

    coeffSuppressor = depositionCoeff(alphaSuppressor, i0Suppressor)
    coeffAccelerator = depositionCoeff(alphaAccelerator, i0Accelerator)

    exchangeCurrentDensity = acceleratorVar.getInterfaceVar() * (coeffAccelerator - coeffSuppressor) + coeffSuppressor

    currentDensity = metalVar / bulkMetalConcentration * exchangeCurrentDensity

    depositionRateVariable = currentDensity * atomicVolume / charge / faradaysConstant

    extensionVelocityVariable = CellVariable(
        name = 'extension velocity',
        mesh = mesh,
        value = depositionRateVariable)   

    from fipy.models.levelSet.surfactant.adsorbingSurfactantEquation \
             import AdsorbingSurfactantEquation

    kAccelerator = rateConstant * numerix.exp(-alphaAdsorption * etaPrime)
    kAcceleratorConsumption =  Bd + A / (numerix.exp(Ba * (overpotential + Vd)) + numerix.exp(Bb * (overpotential + Vd)))
    q = m * overpotential + b

    levelerSurfactantEquation = AdsorbingSurfactantEquation(
        levelerVar,
        distanceVar = distanceVar,
        bulkVar = bulkLevelerVar,
        rateConstant = kLeveler,
        consumptionCoeff = kLevelerConsumption * depositionRateVariable)

    accVar1 = acceleratorVar.getInterfaceVar()
    accVar2 = (accVar1 > 0) * accVar1
    accConsumptionCoeff = kAcceleratorConsumption * (accVar2**(q - 1))

    acceleratorSurfactantEquation = AdsorbingSurfactantEquation(
        acceleratorVar,
        distanceVar = distanceVar,
        bulkVar = bulkAcceleratorVar,
        rateConstant = kAccelerator,
        otherVar = levelerVar,
        otherBulkVar = bulkLevelerVar,
        otherRateConstant = kLeveler,
        consumptionCoeff = accConsumptionCoeff)

    from fipy.models.levelSet.advection.higherOrderAdvectionEquation \
         import buildHigherOrderAdvectionEquation

    advectionEquation = buildHigherOrderAdvectionEquation(
        advectionCoeff = extensionVelocityVariable)

    from fipy.boundaryConditions.fixedValue import FixedValue
    from fipy.models.levelSet.electroChem.metalIonDiffusionEquation \
         import buildMetalIonDiffusionEquation

    metalEquation = buildMetalIonDiffusionEquation(
        ionVar = metalVar,
        distanceVar = distanceVar,
        depositionRate = depositionRateVariable,
        diffusionCoeff = metalDiffusionCoefficient,
        metalIonMolarVolume = atomicVolume)

    metalEquationBCs = (
        FixedValue(
        mesh.getTopFaces(),
        bulkMetalConcentration),)

    from fipy.models.levelSet.surfactant.surfactantBulkDiffusionEquation \
         import buildSurfactantBulkDiffusionEquation

    bulkAcceleratorEquation = buildSurfactantBulkDiffusionEquation(
        bulkVar = bulkAcceleratorVar,
        distanceVar = distanceVar,
        surfactantVar = acceleratorVar,
        otherSurfactantVar = levelerVar,
        diffusionCoeff = acceleratorDiffusionCoefficient,
        rateConstant = kAccelerator * siteDensity)

    bulkAcceleratorEquationBCs = (FixedValue(
        mesh.getTopFaces(),
        bulkAcceleratorConcentration),)

    bulkLevelerEquation = buildSurfactantBulkDiffusionEquation(
        bulkVar = bulkLevelerVar,
        distanceVar = distanceVar,
        surfactantVar = levelerVar,
        diffusionCoeff = levelerDiffusionCoefficient,
        rateConstant = kLeveler * siteDensity)

    bulkLevelerEquationBCs =  (FixedValue(
        mesh.getTopFaces(),
        bulkLevelerConcentration),)

    eqnTuple = ( (advectionEquation, distanceVar, ()),
                 (levelerSurfactantEquation, levelerVar, ()),
                 (acceleratorSurfactantEquation, acceleratorVar, ()),
                 (metalEquation, metalVar,  metalEquationBCs),
                 (bulkAcceleratorEquation, bulkAcceleratorVar, bulkAcceleratorEquationBCs),
                 (bulkLevelerEquation, bulkLevelerVar, bulkLevelerEquationBCs))

    levelSetUpdateFrequency = int(0.7 * narrowBandWidth / cellSize / cflNumber / 2)

    totalTime = 0.0

    if displayViewers:
        from fipy.viewers.mayaviViewer.mayaviSurfactantViewer import MayaviSurfactantViewer
        viewers = (
            MayaviSurfactantViewer(distanceVar, acceleratorVar.getInterfaceVar(), zoomFactor = 1e6, limits = { 'datamax' : 0.5, 'datamin' : 0.0 }, smooth = 1, title = 'accelerator coverage'),
            MayaviSurfactantViewer(distanceVar, levelerVar.getInterfaceVar(), zoomFactor = 1e6, limits = { 'datamax' : 0.5, 'datamin' : 0.0 }, smooth = 1, title = 'leveler coverage'))
        
    for step in range(numberOfSteps):

        if displayViewers:
            if step % displayRate == 0:
                for viewer in viewers:
                    viewer.plot()
            
        if step % levelSetUpdateFrequency == 0:
            distanceVar.calcDistanceFunction(deleteIslands = True)

        extensionVelocityVariable.setValue(depositionRateVariable)

        extOnInt = numerix.where(distanceVar > 0,
                                 numerix.where(distanceVar < 2 * cellSize,
                                               extensionVelocityVariable,
                                               0),
                                 0)

        dt = cflNumber * cellSize / max(extOnInt)

        id = max(numerix.nonzero(distanceVar._getInterfaceFlag()))
        distanceVar.extendVariable(extensionVelocityVariable, deleteIslands = True)

        extensionVelocityVariable[mesh.getFineMesh().getNumberOfCells():] = 0.

        for eqn, var, BCs in eqnTuple:
            var.updateOld()

        for eqn, var, BCs in eqnTuple:
            if isinstance(eqn, AdsorbingSurfactantEquation):
                eqn.solve(var, boundaryConditions = BCs, dt = dt)
                res = (0.,0.)
            else:
                res, = eqn.solve(var, boundaryConditions = BCs, dt = dt, returnItems = ('residual',))

        totalTime += dt

    testFile = 'testLeveler.gz'
    import os
    import examples.levelSet.electroChem
    filepath = os.path.join(examples.levelSet.electroChem.__path__[0], testFile)
    from fipy.tools import dump
    data = dump.read(filepath)

    print numerix.allclose(data, levelerVar)

##    viewer[0].plot('accelerator.png')
##    viewer[1].plot('leveler.png')
    
if __name__ == '__main__':
    runLeveler()
    raw_input("finished")    
