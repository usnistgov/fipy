#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "electroChemicalSystem.py"
 #                                    created: 8/26/04 {10:29:10 AM} 
 #                                last update: 8/26/04 {4:01:07 PM} { 1:23:41 PM}
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

"""

__docformat__ = 'restructuredtext'

import Numeric

from fipy.meshes.grid2D import Grid2D
from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
from fipy.models.levelSet.distanceFunction.distanceEquation import DistanceEquation
from fipy.models.levelSet.surfactant.surfactantVariable import SurfactantVariable
from fipy.variables.cellVariable import CellVariable
from fipy.models.levelSet.electroChem.depositionRateVariable import DepositionRateVariable
from fipy.models.levelSet.surfactant.adsorbingSurfactantEquation import AdsorbingSurfactantEquation
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.models.levelSet.advection.higherOrderAdvectionEquation  import HigherOrderAdvectionEquation
from fipy.models.levelSet.distanceFunction.extensionEquation import ExtensionEquation
from fipy.models.levelSet.electroChem.metalIonDiffusionEquation import MetalIonDiffusionEquation
from fipy.iterators.iterator import Iterator
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from parameters import parameters
from fipy.models.levelSet.surfactant.surfactantBulkDiffusionEquation import SurfactantBulkDiffusionEquation

controlParameters = parameters['control parameters']
geometryParameters = parameters['geometry parameters']
experimentalParameters = parameters['experimental parameters']

yCells = controlParameters["cells below trench"] + int((geometryParameters["trench depth"] + geometryParameters["boundary layer depth"]) / controlParameters["cell size"])
xCells = int((geometryParameters["trench spacing"] / 2) / controlParameters["cell size"])

mesh = Grid2D(dx = controlParameters["cell size"],
              dy = controlParameters["cell size"],
              nx = xCells,
              ny = yCells)

bottomHeight = controlParameters["cells below trench"] * controlParameters['cell size']
trenchHeight = bottomHeight + geometryParameters["trench depth"]
trenchWidth = geometryParameters["trench depth"] / geometryParameters["aspect ratio"]
sideWidth = (geometryParameters["trench spacing"] - trenchWidth) / 2

def electrolyteFunc(cell):
    x,y = cell.getCenter()
    
    if y > trenchHeight:
        return 1
    elif y < bottomHeight:
        return 0
    elif x < sideWidth:
        return 0
    else:
        return 1
            
distanceVar = DistanceVariable(
    name = 'distance variable',
    mesh = mesh,
    value = -1.
    )

electrolyteCells = mesh.getCells(electrolyteFunc)

distanceVar.setValue(1., electrolyteCells)
terminationValue = controlParameters["number of cells in narrow band"] / 2 * controlParameters['cell size']
distanceEquation = DistanceEquation(distanceVar, terminationValue = terminationValue)
distanceEquation.solve(terminationValue = 1e+10)
        
acceleratorVar = SurfactantVariable(
    name = "accelerator variable",
    value = experimentalParameters["initial accelerator coverage"],
    distanceVar = distanceVar
    )

bulkAcceleratorVar = CellVariable(
    name = 'bulk accelerator variable',
    mesh = mesh,
    value = experimentalParameters['bulk accelerator concentration']
    )
        
metalIonVariable = CellVariable(
    name = 'metal ion variable',
    mesh = mesh,
    value = experimentalParameters['bulk metal ion concentration']
    )
    
depositionRateVariable = DepositionRateVariable(
    metalIonVariable,
    name = "deposition rate variable",
    acceleratorVariable = acceleratorVar,
    parameters = parameters
    )

extensionVelocityVariable = CellVariable(
    name = 'extension velocity',
    mesh = mesh,
    value = Numeric.array(depositionRateVariable)
    )

constant = parameters['accelerator properties']['rate constant']['constant']
overPotentialDependence = parameters['accelerator properties']['rate constant']['overpotential dependence']
overPotential = experimentalParameters['overpotential']

rateConstant = constant + overPotentialDependence * overPotential**3

surfactantEquation = AdsorbingSurfactantEquation(
    acceleratorVar,
    distanceVar,
    experimentalParameters['bulk accelerator concentration'],
    rateConstant
    )

advectionEquation = HigherOrderAdvectionEquation(
    distanceVar,
    advectionCoeff = extensionVelocityVariable)
        
extensionEquation = ExtensionEquation(
    distanceVar,
    extensionVelocityVariable,
    terminationValue = terminationValue)
        
metalIonEquation = MetalIonDiffusionEquation(
    metalIonVariable,
    distanceVar = distanceVar,
    depositionRate = depositionRateVariable,
    diffusionCoeff = parameters['metal ion properties']['diffusion coefficient'],
    metalIonAtomicVolume = parameters['metal ion properties']['atomic volume'],
    boundaryConditions = (
        FixedValue(
            mesh.getFacesTop(),
            parameters['experimental parameters']['bulk metal ion concentration'],
            ),
        )
    )

bulkAcceleratorEquation = SurfactantBulkDiffusionEquation(
    bulkAcceleratorVar,
    distanceVar = distanceVar,
    surfactantVar = acceleratorVar,
    diffusionCoeff = parameters['accelerator properties']['diffusion coefficient'],
    rateConstant = rateConstant * parameters['accelerator properties']['site density'],
    boundaryConditions = (
        FixedValue(
            mesh.getFacesTop(),
            parameters['experimental parameters']['bulk accelerator concentration']
            ),
        ) 
    )
    

iterator = Iterator((extensionEquation, advectionEquation, surfactantEquation, metalIonEquation, bulkAcceleratorEquation))

class SurfactantPlotVariable(CellVariable):
    def __init__(self, var = None, name = 'surfactant'):
        CellVariable.__init__(self, mesh = var.getMesh(), name = name)
        self.var = self.requires(var)

    def _calcValue(self):
        self.value = Numeric.array(self.var.getInterfaceValue())

accMod = SurfactantPlotVariable(acceleratorVar)

metalIonPlotVar = (distanceVar < 0) *  experimentalParameters['bulk metal ion concentration'] + metalIonVariable * (distanceVar > 0 )

bulkAcceleratorPlotVar = (distanceVar < 0) *  experimentalParameters['bulk accelerator concentration'] + bulkAcceleratorVar * (distanceVar > 0 )

res = 3
cells = yCells * 2**(res-1)

viewers = (
    Grid2DGistViewer(var = metalIonPlotVar,
                     palette = 'rainbow.gp',
                     maxVal = experimentalParameters['bulk metal ion concentration'],
                     resolution = res,
                     grid = 0,
                     limits = (0, cells, 0, cells),
                     dpi = 100),

    Grid2DGistViewer(var = bulkAcceleratorPlotVar,
                     palette = 'rainbow.gp',
                     maxVal = experimentalParameters['bulk accelerator concentration'],
                     resolution = res,
                     grid = 0,
                     limits = (0, cells, 0, cells),
                     dpi = 100),
    
    Grid2DGistViewer(var = accMod,
                     palette = 'heat.gp',
                     grid = 0,
                     resolution = res,
                     limits = (0, cells, 0, cells),
                     dpi = 75),

    Grid2DGistViewer(var = distanceVar,
                     palette = 'heat.gp',
                     minVal = -1e-8,
                     maxVal = 1e-8,
                     resolution = res,
                     grid = 0,
                     limits = (0, cells, 0, cells),
                     dpi = 100)

            )

            
def plotVariables():

    for viewer in viewers:
        viewer.plot()

def calcTimeStep():
    argmax = Numeric.argmax(extensionVelocityVariable)
    print 'max extension velocity',extensionVelocityVariable[argmax]
    return controlParameters['cfl number'] * controlParameters['cell size'] / extensionVelocityVariable[argmax]

from fipy.viewers.pyxviewer import Grid2DPyxViewer

if __name__ == '__main__':

    plotVariables()

    raw_input("press key to continue")

    totalTime = 0.
    step = 0
        
    levelSetUpdateFrequency = int(0.8 * terminationValue / controlParameters['cell size'] / controlParameters['cfl number'])

    import fipy.tools.dump as dump

    while totalTime < controlParameters['simulation end time']:

        print 'step',step,'time',totalTime
        
        if step % levelSetUpdateFrequency == 0:
            distanceEquation.solve()

        extensionVelocityVariable.setValue(Numeric.array(depositionRateVariable))
        dt = calcTimeStep()
        totalTime += dt
        iterator.timestep(dt = dt)
        
        plotVariables()

        if step % controlParameters['variable dump frequency'] == 0:

            dumpVars = (
                viewers[0].getArray(),
                viewers[1].getArray(),
                viewers[2].getArray(),
                viewers[3].getArray()
                )

            filename = 'dump' + str(step)

            dump.write(dumpVars, filename)

        step += 1
                
    raw_input('finished')
