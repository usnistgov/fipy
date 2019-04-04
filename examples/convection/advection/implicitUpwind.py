"""
This example shows the failure of advecting a square pulse with a first
order implicit upwind scheme.
"""

from fipy import CellVariable, Grid1D, TransientTerm, PowerLawConvectionTerm, LinearLUSolver, Viewer
from fipy.tools import numerix

valueLeft = 0.
valueRight = 0.
L = 10.
nx = 400
dx = L / nx
cfl = 0.01
velocity = 1.
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

eq = TransientTerm() - PowerLawConvectionTerm(coeff = (velocity,))

if __name__ == '__main__':

    viewer = Viewer(vars=(var,))
    viewer.plot()
    raw_input("press key to continue")
    for step in range(steps):
        eq.solve(var,
                 dt = timeStepDuration,
                 solver = LinearLUSolver(tolerance = 1.e-15))
        viewer.plot()
    viewer.plot()
    raw_input('finished')
