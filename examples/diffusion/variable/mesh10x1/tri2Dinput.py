#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 9/3/04 {10:43:18 PM} 
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

"""

This example is a 1D steady state diffusion test case as in
`./examples/diffusion/variable/mesh2x1/input.py` with then
number of cells set to `nx = 10`.

A simple analytical answer can be used to test the result:

   >>> x = mesh.getCellCenters()[:,0]
   >>> values = Numeric.where(x < 3. * L / 4., 10 * x - 9. * L / 4., x + 18. * L / 4.)
   >>> values = Numeric.where(x < L / 4., x, values)
   >>> Numeric.allclose(values, Numeric.array(var), atol = 1e-8, rtol = 1e-8)
   1

"""

import Numeric

from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.equations.diffusionEquation import DiffusionEquation
from fipy.iterators.iterator import Iterator
from fipy.meshes.numMesh.tri2D import Tri2D
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.variables.cellVariable import CellVariable
from fipy.viewers.pyxviewer import PyxViewer

nx = 10
ny = 1

valueLeft = 0.
fluxRight = 1.
timeStepDuration = 1. 

L = 10.

dx = L / nx
dy = 1.

mesh = Tri2D(dx, dy, nx, ny)
    
var = CellVariable(
    name = "solution variable",
    mesh = mesh,
    value = valueLeft)

x = mesh.getFaceCenters()[:,0]
middleFaces = Numeric.logical_or(x < L / 4.,x >= 3. * L / 4.)
diffCoeff = Numeric.where(middleFaces, 1., 0.1)

eq = DiffusionEquation(
    var,
    transientCoeff = 0. / timeStepDuration, 
    diffusionCoeff = diffCoeff,
    solver = LinearPCGSolver(
    tolerance = 1.e-15, 
    steps = 1000
    ),
    boundaryConditions=(
    FixedValue(mesh.getFacesLeft(),valueLeft),
    FixedFlux(mesh.getFacesRight(),fluxRight),
    FixedFlux(mesh.getFacesTop(),0.),
    FixedFlux(mesh.getFacesBottom(),0.)
    )
    )

it = Iterator((eq,))

it.timestep()

if __name__ == '__main__':
    viewer = PyxViewer(var)
    viewer.plot()
    raw_input('finished')
