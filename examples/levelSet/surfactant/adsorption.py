#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "adsorption.py"
 #                                    created: 9/10/04 {3:23:47 PM}
 #                                last update: 8/2/05 {4:58:59 PM} 
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
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

r"""

This example tests 1D adsorption onto an interface and subsequent
depletion from the bulk. The governing equations are given by,

.. raw:: latex

    $$ c_t = D c_{xx} $$

.. raw:: latex

    $$ D c_x = \Gamma k c (1 - \theta) \;\; \text{at} \;\; x = 0 $$

and

.. raw:: latex

    $$ c = c^{\infty} \;\; \text{at} \;\; x = L $$

and on the interface

.. raw:: latex

    $$ D c_x = -k c (1 - \theta) \;\; \text{at} \;\; x = 0$$

There is a dimensionless number, M, that goverens whether the system
is in an interface limited (M>>1) or diffusion limited (M<<1)
regime. There are analytical solutions for both regimes. The
dimensionless number is given by:

.. raw:: latex

    $$ M = \frac{D}{L^2 k cinf} $$


The test solution provided here is for the case of interface limited
kinetics. The analytical solutions are given by,

.. raw:: latex

    $$ -D \ln \left( 1 - \theta \right) + k L \Gamma_0 \theta = \frac{k D c^{\infty} t}{\Gamma_0} $$

and

.. raw:: latex

    $$ c(x) = \frac{c^{\infty} \left[ k \Gamma_0 (1 - \theta) x / D \right]}{1 + k \Gamma_0 (1 - \theta) L / D$$

Make sure the dimensionless parameter is large enough

   >>> (diffusion / cinf / L / L / rateConstant) > 100
   True
   
Start time steping:

   >>> currentTime = 0.
   >>> for i in range(totalTimeSteps):
   ...     surfactantVar.updateOld()
   ...     bulkVar.updateOld()
   ...     surfEqn.solve(surfactantVar, dt = dt)
   ...     bulkEqn.solve(bulkVar, dt = dt, boundaryConditions = bcs)
   ...     currentTime += dt

Compare the analaytical and numerical results:

   >>> theta = surfactantVar.getInterfaceVar()[1]

   >>> numerix.allclose(currentTimeFunc(theta), currentTime, rtol = 1e-4)
   1
   >>> numerix.allclose(concentrationFunc(theta), bulkVar[1:], rtol = 1e-4)
   1


"""
__docformat__ = 'restructuredtext'

from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
from fipy.variables.cellVariable import CellVariable
from fipy.models.levelSet.surfactant.surfactantVariable import SurfactantVariable
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.models.levelSet.surfactant.surfactantBulkDiffusionEquation import buildSurfactantBulkDiffusionEquation
from fipy.models.levelSet.surfactant.adsorbingSurfactantEquation import AdsorbingSurfactantEquation

# parameter values

diffusion = 101.
rateConstant = .5
L = 1.
siteDensity = 1.
cinf = 1.
nx = 100
totalTimeSteps = 100
dt = 0.001

## build the mesh

from fipy.meshes.grid1D import Grid1D
dx = L / (nx - 1.5)
mesh = Grid1D(nx = nx, dx = dx)

## build the distance variable


value = mesh.getCellCenters()[:,0] - 1.499 * dx
##distanceVar = DistanceVariable(mesh = mesh, value = dx * (numerix.arange(nx) - 0.999))
distanceVar = DistanceVariable(mesh = mesh, value = value, hasOld = 1)

## Build the bulk diffusion equation

bulkVar = CellVariable(mesh = mesh, value = cinf)

surfactantVar = SurfactantVariable(distanceVar = distanceVar)

bulkEqn = buildSurfactantBulkDiffusionEquation(bulkVar,
                                          distanceVar = distanceVar,
                                          surfactantVar = surfactantVar,
                                          diffusionCoeff = diffusion,
                                          rateConstant = rateConstant * siteDensity)

bcs = (FixedValue(mesh.getFacesRight(), cinf),)

## Build the surfactant equation

surfEqn = AdsorbingSurfactantEquation(surfactantVar = surfactantVar,
                                      distanceVar = distanceVar,
                                      bulkVar = bulkVar,
                                      rateConstant = rateConstant)

## Build the analytical solutions,

x = mesh.getCellCenters()[1:,0] - dx

def concentrationFunc(theta):
    tmp = (1 + rateConstant * siteDensity * (1 - theta) * L / diffusion)
    return cinf * (1 + rateConstant * siteDensity * (1 - theta) * x / diffusion) / tmp

import fipy.tools.numerix as numerix

def currentTimeFunc(theta):
    tmp = -diffusion * numerix.log(1 - theta) + rateConstant * siteDensity * L * theta
    return tmp / rateConstant / diffusion/ cinf

## set up the comparison arrays

theta = surfactantVar.getInterfaceVar()[1]
    

if __name__ == "__main__":

    ## set up the viewers

    import fipy.viewers

    ## start time stepping

    currentTime = 0.

    for i in range(totalTimeSteps):

        ## evaluate the analytical and numerical solution and plot

        theta = surfactantVar.getInterfaceVar()[1]
        print "theta:",theta

        ## do a time step
        surfactantVar.updateOld()
        bulkVar.updateOld()
        surfEqn.solve(surfactantVar, dt = dt)
        bulkEqn.solve(bulkVar, dt = dt, boundaryConditions = bcs)
        currentTime += dt

    raw_input("finished")
