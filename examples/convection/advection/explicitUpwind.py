"""
This example shows the failure of advecting a square pulse with a first
order explicit upwind scheme.
"""
from __future__ import division
from __future__ import unicode_literals

from builtins import range
from fipy import input
from fipy import CellVariable, Grid1D, TransientTerm, ExplicitUpwindConvectionTerm, LinearLUSolver, Viewer
from fipy.tools import numerix

valueLeft = 0.
valueRight = 0.
L = 10.
nx = 400
dx = L / nx
cfl = 0.1
velocity = -1.
timeStepDuration = cfl * dx / abs(velocity)
steps = 1000

mesh = Grid1D(dx = dx, nx = nx)

startingArray = numerix.zeros(nx, 'd')
startingArray[50:90] = 1.

var = CellVariable(
    name = "advection variable",
    mesh = mesh,
    value = startingArray)

var.constrain(valueLeft, mesh.facesLeft)
var.constrain(valueRight, mesh.facesRight)

eq = TransientTerm() - ExplicitUpwindConvectionTerm(coeff=(velocity,))

if __name__ == '__main__':

    viewer = Viewer(vars=(var,))
    for step in range(steps):
        eq.solve(var,
                 dt = timeStepDuration,
                 solver = LinearLUSolver(tolerance=1.e-15, iterations=2000))
        viewer.plot()
    viewer.plot()
    input('finished')
