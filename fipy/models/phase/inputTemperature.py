#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 12/29/03 {2:45:50 PM} 
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

from meshes.grid2D import Grid2D
from type2PhaseEquation import Type2PhaseEquation
from solvers.linearPCGSolver import LinearPCGSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator
from viewers.grid2DGistViewer import Grid2DGistViewer
from variables.cellVariable import CellVariable
from modularVariable import ModularVariable
from profiler.profiler import Profiler
from profiler.profiler import calibrate_profiler
from temperatureEquation import TemperatureEquation

import Numeric

##timeStepDuration = 0.02
timeStepDuration = 1.

phaseParameters = {
    'tau'                   : 3e-4,
    'time step duration'    : timeStepDuration,
    'epsilon'               : 0.008,
    's'                     : 0.01,
    'alpha'                 : 0.015,
    'c2'                    : 0.02,
    'anisotropy'            : 0.0,
    'symmetry'              : 4.,
    'kappa 1'               : 0.9,
    'kappa 2'               : 20.
    }

temperatureParameters = {
    'time step duration'    : timeStepDuration,
    'temperature diffusion' : 2.25,
    'latent heat'           : 1.,
    'heat capacity'         : 1.
    }


##Length = 0.5
Length = 20.
nx = 20
ny = 20
dx = Length / nx
dy = Length / ny


mesh = Grid2D(dx,dy,nx,ny)
print "built mesh"

phase = CellVariable(
    name = 'PhaseField',
    mesh = mesh,
    value = 0.
    )

theta = ModularVariable(
    name = 'Theta',
    mesh = mesh,
    value = 0.,
    hasOld = 0
    )

temperature = ModularVariable(
    name = 'Theta',
    mesh = mesh,
    value = -0.5,
    hasOld = 0
    )

phaseViewer = Grid2DGistViewer(var = phase)
temperatureViewer = Grid2DGistViewer(var = temperature, minVal = -0.5, maxVal =0.5)

phaseFields = {
    'theta' : theta,
    'temperature' : temperature
    }

temperatureFields = {
    'phase' : phase
    }

def circleCells(cell,L = Length):
    x = cell.getCenter()
    r = L / 4.
    c = (L / 2., L / 2.)
    if (x[0] - c[0])**2 + (x[1] - c[1])**2 < r**2:
        return 1
    else:
        return 0

interiorCells = mesh.getCells(circleCells)

phase.setValue(1.,interiorCells)

print "building phase equation"

phaseEq = Type2PhaseEquation(
    phase,
    solver = LinearPCGSolver(
	tolerance = 1.e-15, 
	steps = 1000
    ),
    boundaryConditions=(
    FixedFlux(mesh.getExteriorFaces(), 0.),
    ),
    parameters = phaseParameters,
    fields = phaseFields
    )

temperatureEq = TemperatureEquation(
    phase,
    solver = LinearPCGSolver(
	tolerance = 1.e-15, 
	steps = 1000
    ),
    boundaryConditions=(
    FixedFlux(mesh.getExteriorFaces(), 0.),
    ),
    parameters = temperatureParameters,
    fields = temperatureFields
    )

it = Iterator((phaseEq, temperatureEq))

# fudge = calibrate_profiler(10000)
# profile = Profiler('profile', fudge=fudge)
print "solving"
phaseViewer.plot()
temperatureViewer.plot()
raw_input()
for i in range(10):
    it.timestep(steps = 10)
    phaseViewer.plot()
    temperatureViewer.plot()
# profile.stop()

raw_input()

