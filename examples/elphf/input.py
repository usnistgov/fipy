#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 12/9/03 {4:51:47 PM} 
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

"""Electrochemical Phase Field input file

    Build a mesh, variable, and diffusion equation with fixed (zero) flux
    boundary conditions at the top and bottom and fixed value boundary
    conditions at the left and right.
    
    Iterates a solution and plots the result with gist.
    
    Iteration is profiled for performance.
"""

from meshes.grid2D import Grid2D
from variables.cellVariable import CellVariable
from viewers.grid2DGistViewer import Grid2DGistViewer
from concentrationEquation import ConcentrationEquation
from solvers.linearCGSSolver import LinearCGSSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator

from profiler.profiler import Profiler
from profiler.profiler import calibrate_profiler

valueLeft=0.3
valueRight=0.6

mesh = Grid2D(
    dx = 1.,
    dy = 1.,
    nx = 4,
    ny = 1)

var1 = CellVariable(
    name = "c1",
    mesh = mesh,
    value = valueLeft,
    viewer = Grid2DGistViewer
    )

var2 = CellVariable(
    name = "c2",
    mesh = mesh,
    value = valueLeft,
    viewer = Grid2DGistViewer
    )
    
parameters = {
    'substitutionals': (var1,var2)
}

eq1 = ConcentrationEquation(
    var = var1,
    diffusivity = 1.,
    solver = LinearCGSSolver(
	tolerance = 1.e-15, 
	steps = 1000
    ),
    boundaryConditions=(
	FixedValue(faces = mesh.getFacesLeft(),value = valueLeft),
	FixedValue(faces = mesh.getFacesRight(),value = valueRight),
	FixedFlux(faces = mesh.getFacesTop(),value = 0.),
	FixedFlux(faces = mesh.getFacesBottom(),value = 0.)
    ),
    parameters = parameters
)

eq2 = ConcentrationEquation(
    var = var2,
    diffusivity = 1.,
    solver = LinearCGSSolver(
	tolerance = 1.e-15, 
	steps = 1000
    ),
    boundaryConditions=(
	FixedValue(faces = mesh.getFacesLeft(),value = valueRight),
	FixedValue(faces = mesh.getFacesRight(),value = valueLeft),
	FixedFlux(faces = mesh.getFacesTop(),value = 0.),
	FixedFlux(faces = mesh.getFacesBottom(),value = 0.)
    ),
    parameters = parameters
)

it = Iterator((eq1,eq2))

for i in range(2):
    # fudge = calibrate_profiler(10000)
    # profile = Profiler('profile', fudge=fudge)
    it.iterate(1,10000.)
    # profile.stop()

    print var1.getValue()
    print var2.getValue()

var1.plot()
var2.plot()



raw_input()

