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

This example takes the user through assembling a simple problem with
FiPy.  It describes a steady 1D diffusion problem with fixed value
boundary conditions such that,

.. raw:: latex

   $$ \partial_t (T \phi) = \partial_x (D \partial_x \phi) $$

with initial conditions,

.. raw:: latex

   $$ \phi = 0 \;\; at \;\; t = 0 \;\; $$

boundary conditions,

.. raw:: latex

   $$ \phi = 0 \;\; at \;\; x = 0 \;\; and \;\; \phi = 1  \;\; at \;\; x = 1 $$

and parameter values,

.. raw:: latex

   $$ T = 0 \;\; and \;\; D = 1 $$

The first step is to create a mesh with 50 elements. The `Grid2D`
object represents a cubic structured grid.

    >>> nx = 10
    >>> ny = 1
    >>> dx = 1.
    >>> dy = 1.
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)

The solution of all equations in FiPy requires a variable They store
values on various parts of the mesh. In this case we need a
`CellVariable` object as the solution is sought on the cell
centers. The boundary conditions are given by `valueLeft = 0` and
`valueRight = 1`. The initial value for the variable is set to `value = valueLeft`.

    >>> valueLeft = 0
    >>> valueRight = 1
    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(name = "solution variable", mesh = mesh, value = valueLeft)

Boundary conditions are given to the equation via a tuple. Boundary
conditions are formed with a value and a set of faces over which they
apply. For example here the exterior faces on the left of the domain
are extracted by `mesh.getFacesLeft()`. These faces and a value
(`valueLeft`) are passed to a `FixedValue` boundary condition. A fixed
flux of zero is set on the top and bottom surfaces to simulate a 1
dimensional problem.

    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> from fipy.boundaryConditions.fixedFlux import FixedFlux
    >>> boundaryConditions = (FixedValue(mesh.getFacesLeft(),valueLeft), FixedValue(mesh.getFacesRight(),valueRight), FixedFlux(mesh.getFacesTop(),0.), FixedFlux(mesh.getFacesBottom(),0.))

A solver is created and passed to the equation. This solver uses an
iterative conjugant gradient method to solve implicitly at each time
step.

    >>> from fipy.solvers.linearPCGSolver import LinearPCGSolver
    >>> solver = LinearPCGSolver(tolerance = 1.e-15, steps = 1000)

An equation knows how to assemble the solution matrix from the
parameters and boundary conditions and how to solve the solution
matrix with the solver. The `transientCoeff` is set to 0 thus
making the problem steady state. The `diffusionCoeff` is set to 1.

    >>> from fipy.equations.diffusionEquation import DiffusionEquation
    >>> eq = DiffusionEquation(var, transientCoeff = 0., diffusionCoeff = 1., solver = solver, boundaryConditions = boundaryConditions)

The `Iterator` object takes a tuple of equations and solves to a
required tolerance for the given equations at each time step.

    >>> from fipy.iterators.iterator import Iterator
    >>> iterator = Iterator((eq,))

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
    >>> Numeric.allclose(Numeric.array(var), analyticalArray, rtol = 1e-10, atol = 1e-10)
    1

"""

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import doctest
    doctest.testmod()

