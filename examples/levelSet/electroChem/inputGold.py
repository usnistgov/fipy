#!/usr/bin/env python

#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "inputSimpleTrenchSystem.py"
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

    \label{inputGold} is a demonstration of the use of \FiPy{} for
    modeling gold superfill. The material properties and experimental
    parameters used are roughly those that have been previously
    published~\cite{NIST:damascene:2005}.

To run this example from the base fipy directory type::
    
    $ examples/levelSet/electroChem/inputGold.py

at the command line. The results of the simulation will be displayed
and the word `finished` in the terminal at the end of the
simulation. The simulation will only run for 10 time steps. In order
to alter the number of timesteps, the python function that
encapsulates the system of equations must first be imported (at the
python command line),

    >>> from examples.levelSet.electroChem.inputGold import runGold

and then the function can be run with a different number of time steps
with the `numberOfSteps` argument as follows,

    >>> runGold(numberOfSteps = 10, runAsTest = True)
    1

However, do not include the `runAsTest` argument.  This example has a
more realistic default boundary layer depth and thus requires `gmsh`
to construct a more complex mesh.

"""
__docformat__ = 'restructuredtext'

def runGold(faradaysConstant = 9.6e4,
            gasConstant = 8.314,
            transferCoefficient = 0.5,
            consumptionRateConstant = 2.6e+6,
            molarVolume = 10.21e-6,
            charge = 1.0,
            metalDiffusion = 1.7e-9,
            temperature = 298.,
            overpotential = -0.3,
            metalConcentration = 20.0,
            catalystCoverage = 0.15,
            currentDensity0 = 0.26,
            currentDensity1 = 45.,
            cellSize = 0.1e-7,
            trenchDepth = 0.2e-6,
            aspectRatio = 1.47,
            trenchSpacing = 0.5e-6,
            boundaryLayerDepth = 90.0e-6,
            numberOfSteps = 40,
            taperAngle = 6.0,
            A = 3e-2,
            B = 6.5e-1,
            runAsTest = False):
    
    cflNumber = 0.2
    numberOfCellsInNarrowBand = 20
    cellsBelowTrench = 10
    
    from fipy.tools import numerix
    from gapFillMesh import TrenchMesh
    mesh = TrenchMesh(cellSize = cellSize,
                      trenchSpacing = trenchSpacing,
                      trenchDepth = trenchDepth,
                      boundaryLayerDepth = boundaryLayerDepth,
                      aspectRatio = aspectRatio,
                      angle = numerix.pi * taperAngle / 180.,
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
    catalystVar = SurfactantVariable(
        name = "catalyst variable",
        value = catalystCoverage,
        distanceVar = distanceVar)

    from fipy.variables.cellVariable import CellVariable
    metalVar = CellVariable(
        name = 'metal variable',
        mesh = mesh,
        value = metalConcentration)

    exchangeCurrentDensity = 16 * (A + B * catalystVar.getInterfaceVar())
    
    currentDensity = metalVar / metalConcentration * exchangeCurrentDensity

    depositionRateVariable = currentDensity * molarVolume / charge / faradaysConstant

    extensionVelocityVariable = CellVariable(
        name = 'extension velocity',
        mesh = mesh,
        value = depositionRateVariable)   

    from fipy.models.levelSet.surfactant.adsorbingSurfactantEquation \
                import AdsorbingSurfactantEquation

    catalystSurfactantEquation = AdsorbingSurfactantEquation(
        catalystVar,
        distanceVar = distanceVar,
        bulkVar = 0,
        rateConstant = 0,
        consumptionCoeff = consumptionRateConstant * extensionVelocityVariable)

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
        diffusionCoeff = metalDiffusion,
        metalIonMolarVolume = molarVolume)

    metalEquationBCs = (
        FixedValue(
        mesh.getTopFaces(),
        metalConcentration),)

    from fipy.models.levelSet.surfactant.surfactantBulkDiffusionEquation \
                    import buildSurfactantBulkDiffusionEquation

    eqnTuple = ( (advectionEquation, distanceVar, ()),
                 (catalystSurfactantEquation, catalystVar, ()),
                 (metalEquation, metalVar,  metalEquationBCs))
                     
    class PlotVariable(CellVariable):
        def __init__(self, var = None, name = ''):
            CellVariable.__init__(self, mesh = mesh.getFineMesh(), name = name)
            self.var = self._requires(var)

        def _calcValue(self):
            self.value = numerix.array(self.var[:self.mesh.getNumberOfCells()])
    
    levelSetUpdateFrequency = int(0.7 * narrowBandWidth / cellSize / cflNumber / 2)
       
    m1 = 0

    if not runAsTest:
        from fipy.viewers import make
        viewers = (make(PlotVariable(var = distanceVar), limits = {'datamax' : 1e-9, 'datamin' : -1e-9}),
                   make(PlotVariable(var = catalystVar.getInterfaceVar())))

    def doOneStep(step, totalTime):
        if not runAsTest:
            print totalTime,' ',step
        
        if step % levelSetUpdateFrequency == 0:
            distanceVar.calcDistanceFunction(deleteIslands = True)
            
        extensionVelocityVariable.setValue(numerix.array(depositionRateVariable))

        argmax = numerix.argmax(extensionVelocityVariable)
        dt = cflNumber * cellSize / extensionVelocityVariable[argmax]
            
        distanceVar.extendVariable(extensionVelocityVariable, deleteIslands = True)

        for eqn, var, BCs in eqnTuple:
            var.updateOld()

        for eqn, var, BCs in eqnTuple:
            eqn.solve(var, boundaryConditions = BCs, dt = dt)

        if not runAsTest:
            for viewer in viewers:
                viewer.plot()

        return totalTime + dt
        
    totalTime = 0.
    step = 0
    
    metalVar.setValue(metalConcentration)

    while step < numberOfSteps:
        totalTime = doOneStep(step, totalTime)
        step += 1

    if runAsTest:
        
        from fipy.tools import dump
        import os
        import examples.levelSet.electroChem
        data = dump.read(os.path.join(examples.levelSet.electroChem.__path__[0], 'goldData.gz'))
        print catalystVar.allclose(data)
    
if __name__ == '__main__':
    runGold()
