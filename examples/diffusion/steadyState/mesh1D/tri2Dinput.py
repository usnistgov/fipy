#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 6/15/04 {10:58:12 AM} 
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
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

"""

To run this example from the base fipy directory type
`./examples/diffusion/steadyState/mesh1D/tri2Dinput.py` at the command
line. A contour plot should appear and the word `finished` in the terminal.

This example is similar to the example found in:
`./examples/diffusion/steadyState/mesh1D/input.py`
However the `mesh` is a `Tri2D` object rather than a `Grid2D` object.

Here the iterator does one time step to implicitly find the steady state
solution.
    
    >>> iterator.timestep()

To test the solution, the analytical result is required. The x
coordinates from the mesh are gathered and the length of the domain,
`Lx`, is calculated.  An array, `analyticalArray`, is calculated to
compare with the numerical result,

    >>> x = mesh.getCellCenters()[:,0]
    >>> Lx = nx * dx
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx

Finally the analytical and numerical results are compared with a
tolerance of `1e-10`. The variable `var` is coerced to a `Numeric.array`
for the comparison.

    >>> import Numeric
    >>> print Numeric.array(var)
    >>> Numeric.allclose(Numeric.array(var), analyticalArray, rtol = 1e-10, atol = 1e-10)
    1

"""

__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.equations.diffusionEquation import DiffusionEquation
from fipy.iterators.iterator import Iterator
from fipy.viewers.pyxviewer import PyxViewer
from fipy.meshes.numMesh.tri2D import Tri2D

nx = 50
ny = 1
dx = 2.
dy = 1.

mesh = Tri2D(dx = dx, dy = dy, nx = nx, ny = ny)

valueLeft = 0
valueRight = 1
var = CellVariable(name = "solution-variable", mesh = mesh, value = valueLeft)

viewer = PyxViewer(var)

boundaryConditions = (FixedValue(mesh.getFacesLeft(),valueLeft), FixedValue(mesh.getFacesRight(),valueRight), FixedFlux(mesh.getFacesTop(),0.), FixedFlux(mesh.getFacesBottom(),0.))

solver = LinearPCGSolver(tolerance = 1.e-15, steps = 1000)

eq = DiffusionEquation(var, transientCoeff = 0., diffusionCoeff = 1., solver = solver, boundaryConditions = boundaryConditions)

iterator = Iterator((eq,))

if __name__ == '__main__':
    iterator.timestep()
    viewer.plot()
    raw_input("finished")
