#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "tri2D.py"
 #
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
 # ###################################################################
 ##

"""

This input file again solves a 1D diffusion problem as in
`./examples/diffusion/steadyState/mesh1D.py`. The difference in this example is
that the solution method is explicit. The equation used is the
`ExplicitDiffusionEquation`. In this case many steps have to be taken to reach
equilibrum. The `timeStepDuration` parameter specifies the size of each time step
and `steps` is the number of time steps.

    >>> dx = 1.
    >>> dy = 1.
    >>> nx = 10
    >>> ny = 1
    >>> valueLeft = 0.
    >>> valueRight = 1.
    >>> timeStepDuration = 0.02
    >>> steps = 10

A loop is required to execute the necessary time steps:

    >>> for step in range(steps):
    ...     eq.solve(var, solver=solver, dt=timeStepDuration)

The result is again tested in the same way:

    >>> Lx = nx * dx
    >>> x = mesh.cellCenters[0]
    >>> print var.allclose(answer, rtol = 1e-8)
    1

"""

from fipy import CellVariable, Tri2D, TransientTerm, ExplicitDiffusionTerm, DefaultSolver, Viewer
from fipy.tools import numerix

dx = 1.
dy = 1.
nx = 10
ny = 1
valueLeft = 0.
valueRight = 1.
timeStepDuration = 0.02

mesh = Tri2D(dx, dy, nx, ny)

var = CellVariable(
    name = "concentration",
    mesh = mesh,
    value = valueLeft)

eq = TransientTerm() == ExplicitDiffusionTerm()

solver = DefaultSolver(tolerance=1e-6, iterations=1000)

var.constrain(valueLeft, mesh.facesLeft)
var.constrain(valueRight, mesh.facesRight)

answer = numerix.array([  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                          0.00000000e+00,  0.00000000e+00,  1.58508452e-07,  6.84325019e-04,
                          7.05111362e-02,  7.81376523e-01,  0.00000000e+00,  0.00000000e+00,
                          0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                          0.00000000e+00,  4.99169535e-05,  1.49682805e-02,  3.82262622e-01,
                          0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                          0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  4.06838361e-06,
                          3.67632029e-03,  1.82227062e-01,  0.00000000e+00,  0.00000000e+00,
                          0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                          0.00000000e+00,  4.99169535e-05,  1.49682805e-02,  3.82262622e-01])

if __name__ == '__main__':
    steps = 1000

    for step in range(steps):
        eq.solve(var, solver = solver, dt = timeStepDuration)
    print var
    viewer = Viewer(vars = var)
    viewer.plot()
    raw_input('finished')
