#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 9/3/04 {10:41:42 PM} { 1:23:41 PM}
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

Also a surfactant is present of the interface, governed by the equation:

.. raw:: latex

    $$ \\frac{d \\theta}{d t} = J v \\theta $$

The result can be tested with the following code:


   >>> surfactantBefore = Numeric.sum(surfactantVariable * mesh.getCellVolumes())
   >>> for step in range(steps):
   ...     surfactantVariable.updateOld()
   ...     distanceVariable.updateOld()
   ...     surfactantEquation.solve(surfactantVariable)
   ...     advectionEquation.solve(distanceVariable, dt = timeStepDuration)
   >>> surfactantEquation.solve(surfactantVariable)
   >>> surfactantAfter = Numeric.sum(surfactantVariable * mesh.getCellVolumes())
   >>> Numeric.allclose(surfactantBefore, surfactantAfter)
   1
   >>> areas = (distanceVariable.getCellInterfaceAreas() < 1e-6) * 1e+10 + distanceVariable.getCellInterfaceAreas()
   >>> answer = initialSurfactantValue * initialRadius / (initialRadius +  distanceToTravel)
   >>> coverage = surfactantVariable * mesh.getCellVolumes() / areas
   >>> error = 0.
   >>> size = 0
   >>> for i in range(len(coverage)):
   ...     if coverage[i] > 1e-3:
   ...         error += (coverage[i] / answer - 1.)**2
   ...         size += 1            
   >>> print Numeric.sqrt(error / size)
   0.00813776069241

"""
__docformat__ = 'restructuredtext'

import Numeric
   
from fipy.meshes.grid2D import Grid2D
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
from fipy.models.levelSet.advection.higherOrderAdvectionEquation import buildHigherOrderAdvectionEquation
from fipy.models.levelSet.surfactant.surfactantEquation import SurfactantEquation
from fipy.models.levelSet.surfactant.surfactantVariable import SurfactantVariable

L = 1.
nx = 50
velocity = 1.
cfl = 0.1
velocity = 1.
distanceToTravel = L / 10.
initialRadius = L / 4.

dx = L / nx
timeStepDuration = cfl * dx / velocity
steps = int(distanceToTravel / dx / cfl)

mesh = Grid2D(dx = dx, dy = dx, nx = nx, ny = nx)

distanceVariable = DistanceVariable(
    name = 'level set variable',
    mesh = mesh,
    value = 1.
    )

cellRadius = Numeric.sqrt((mesh.getCellCenters()[:,0] - L / 2.)**2 + (mesh.getCellCenters()[:,1] - L / 2.)**2)
distanceVariable.setValue(cellRadius - initialRadius)

initialSurfactantValue =  1.

surfactantVariable = SurfactantVariable(
    value = initialSurfactantValue,
    distanceVar = distanceVariable
    )

advectionEquation = buildHigherOrderAdvectionEquation(
    advectionCoeff = velocity)

surfactantEquation = SurfactantEquation(
    distanceVar = distanceVariable)

if __name__ == '__main__':
    
    distanceViewer = Grid2DGistViewer(var = distanceVariable, palette = 'rainbow.gp', minVal = -initialRadius, maxVal = initialRadius)
    surfactantViewer = Grid2DGistViewer(var = surfactantVariable, palette = 'rainbow.gp', minVal = -1., maxVal = 100.)
    distanceViewer.plot()
    surfactantViewer.plot()

    print 'total surfactant before:',Numeric.sum(surfactantVariable * mesh.getCellVolumes())
    
    for step in range(steps):
        surfactantVariable.updateOld()
        distanceVariable.updateOld()
        surfactantEquation.solve(surfactantVariable)
        advectionEquation.solve(distanceVariable, dt = timeStepDuration)
        distanceViewer.plot()
        surfactantViewer.plot()
    surfactantEquation.solve(surfactantVariable)


    print 'total surfactant after:',Numeric.sum(surfactantVariable * mesh.getCellVolumes())

    areas = (distanceVariable.getCellInterfaceAreas() < 1e-6) * 1e+10 + distanceVariable.getCellInterfaceAreas()
    answer = initialSurfactantValue * initialRadius / (initialRadius +  distanceToTravel)
    coverage = surfactantVariable * mesh.getCellVolumes() / areas

    error = 0.
    size = 0
    for i in range(len(coverage)):
        if coverage[i] > 1e-3:
            error += (coverage[i] / answer - 1.)**2
            size += 1
            
    error = Numeric.sqrt(error / size)
    
    print 'error:', error
    
    raw_input('finished')
