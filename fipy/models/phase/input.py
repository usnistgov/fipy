#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 11/30/03 {12:39:40 AM} 
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

"""Diffusion equation input file

    Build a mesh, variable, and diffusion equation with fixed (zero) flux
    boundary conditions at the top and bottom and fixed value boundary
    conditions at the left and right.
    
    Iterates a solution and plots the result with gist.
    
    Iteration is profiled for performance.
"""

from meshes.grid2D import Grid2D
from phaseEquation import PhaseEquation
from solvers.linearPCGSolver import LinearPCGSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator
from viewers.grid2DGistViewer import Grid2DGistViewer
from variables.variable import Variable
from profiler.profiler import Profiler
from profiler.profiler import calibrate_profiler

phaseParameters={
    'tau' :        0.1,
    'epsilon' :    0.008,
    's' :          0.01,
    'alpha' :      0.015,
    'c2':          0.0,
    'anisotropy':  0.,
    'symmetry':    4.
    }

valueLeft=1.
valueRight=1.

L = 1.5
nx = 100
dx = L/nx

mesh = Grid2D(dx,1.,nx,1)

phase = Variable(
    name = 'PhaseField',
    mesh = mesh,
    value = 1.,
    viewer = Grid2DGistViewer
    )

theta = Variable(
    name = 'Theta',
    mesh = mesh,
    value = 1.,
    viewer = Grid2DGistViewer,
    hasOld = 0
    )

def func(x):
    if x[0] > L / 2.:
        return 1
    else:
        return 0

rightCells = mesh.getCells(func)

theta.setValue(0.,rightCells)

eq = PhaseEquation(
    phase,
    theta = theta,
    temperature = 1.,
    solver = LinearPCGSolver(
	tolerance = 1.e-15, 
	steps = 1000
    ),
    boundaryConditions=(
    FixedValue(mesh.getFacesLeft(),valueLeft),
    FixedValue(mesh.getFacesRight(),valueRight)),
    parameters = phaseParameters
    )

it = Iterator((eq,))

fudge = calibrate_profiler(10000)
profile = Profiler('profile', fudge=fudge)
it.iterate(10,0.02)
profile.stop()

print phase.getArray()

phase.plot()

raw_input()

