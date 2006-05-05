#!/usr/bin/env python

## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "cahnHilliard.py"
 #                                     created: 1/18/06 {2:36:12 PM}
 #                                 last update: 1/18/06 {4:03:44 PM}
 # Author: Jonathan Guyer
 # E-mail: <guyer@nist.gov>
 # Author: Daniel Wheeler
 # E-mail: daniel.wheeler@nist.gov
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards and
 # Technology by employees of the Federal Government in the course of their
 # official duties.  Pursuant to title 17 Section 105 of the United States
 # Code this document is not subject to copyright protection and is in the
 # public domain.  cahnHilliard.py is an experimental work.  NIST assumes
 # no responsibility whatsoever for its use by other parties, and makes no
 # guarantees, expressed or implied, about its quality, reliability, or any
 # other characteristic.  We would appreciate acknowledgement if the 
 # document is used.
 # 
 # This document can be redistributed and/or modified freely provided that 
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and 
 #  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
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
This example benchmarks the speed and memory usage of solving the
Cahn-Hilliard equation. Run:
    
    $ python setup.py efficiency_test
"""
__docformat__ = 'restructuredtext'

import time

from fipy.tools.parser import parse

from benchmarker import Benchmarker
bench = Benchmarker()

numberOfElements = parse('--numberOfElements', action = 'store', type = 'int', default = 400)


bench.start()

import fipy.tools.numerix as numerix

nx = int(numerix.sqrt(numberOfElements))
ny = int(numerix.sqrt(numberOfElements))

steps = 10

dx = 2.
dy = 2.

L = dx * nx

asq = 1.0
epsilon = 1
diffusionCoeff = 1

from fipy.meshes.grid2D import Grid2D
mesh = Grid2D(dx, dy, nx, ny)

bench.stop('mesh')

bench.start()

from fipy.variables.cellVariable import CellVariable
import RandomArray

var = CellVariable(name = "phase field",
                   mesh = mesh,
                   value = RandomArray.random(nx * ny))

bench.stop('variables')

bench.start()

from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fipy.terms.transientTerm import TransientTerm

faceVar = var.getArithmeticFaceValue()
doubleWellDerivative = asq * ( 1 - 6 * faceVar * (1 - faceVar))

diffTerm2 = ImplicitDiffusionTerm(coeff = (diffusionCoeff * doubleWellDerivative,))
diffTerm4 = ImplicitDiffusionTerm(coeff = (diffusionCoeff, -epsilon**2))
eqch = TransientTerm() - diffTerm2 - diffTerm4

bench.stop('terms')

bench.start()

from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.solvers.linearLUSolver import LinearLUSolver
##solver = LinearLUSolver(tolerance = 1e-15,steps = 1000)
solver = LinearPCGSolver(tolerance = 1e-15,steps = 1000)

bench.stop('solver')

bench.start()

from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.boundaryConditions.nthOrderBoundaryCondition import NthOrderBoundaryCondition
BCs = (FixedFlux(mesh.getFacesRight(), 0),
       FixedFlux(mesh.getFacesLeft(), 0),
       NthOrderBoundaryCondition(mesh.getFacesLeft(), 0, 3),
       NthOrderBoundaryCondition(mesh.getFacesRight(), 0, 3),
       NthOrderBoundaryCondition(mesh.getFacesTop(), 0, 3),
       NthOrderBoundaryCondition(mesh.getFacesBottom(), 0, 3))

bench.stop('BCs')

dexp=-5

dt = numerix.exp(dexp)
dt = min(100, dt)
dexp += 0.01
var.updateOld()
eqch.solve(var, boundaryConditions = BCs, solver = solver, dt = dt)

bench.start()

for step in range(steps):
    dt = numerix.exp(dexp)
    dt = min(100, dt)
    dexp += 0.01
    var.updateOld()
    eqch.solve(var, boundaryConditions = BCs, solver = solver, dt = dt)
            
bench.stop('solve')

print bench.report(numberOfElements=numberOfElements, steps=steps)
