##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

"""
Test case for the `Grid3D`.

   >>> DiffusionTerm().solve(var)
   >>> DiffusionTerm().solve(var2)
   >>> a = numerix.array(var.globalValue)
   >>> b = numerix.array(var2.globalValue)
   >>> c = numerix.ravel(numerix.array((b, b, b)))
   >>> print(numerix.allclose(a, c))
   1

"""
from __future__ import unicode_literals

from fipy import input
from fipy import CellVariable, Grid2D, Grid3D, DiffusionTerm, Viewer
from fipy.tools import numerix

nx = 8 # FIXME: downsized temporarily from 10 due to https://github.com/usnistgov/fipy/issues/622
ny = 5
nz = 3

dx = 1.
dy = 1.
dz = 1.

valueBottomTop = 0.
valueLeftRight = 1.

mesh = Grid3D(dx = dx, dy = dy, dz = dz, nx = nx, ny = ny, nz = nz)

var = CellVariable(name = "solution variable",
                   mesh = mesh,
                   value = valueBottomTop)

var.constrain(valueLeftRight, mesh.facesLeft)
var.constrain(valueLeftRight, mesh.facesRight)
var.constrain(valueBottomTop, mesh.facesTop)
var.constrain(valueBottomTop, mesh.facesBottom)

#do the 2D problem for comparison

nx = 8 # FIXME: downsized temporarily from 10 due to https://github.com/usnistgov/fipy/issues/622
ny = 5

dx = 1.
dy = 1.

mesh2 = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)

var2 = CellVariable(name = "solution variable 2D",
                    mesh = mesh2,
                    value = valueBottomTop)

var2.constrain(valueLeftRight, mesh2.facesLeft)
var2.constrain(valueLeftRight, mesh2.facesRight)
var2.constrain(valueBottomTop, mesh2.facesTop)
var2.constrain(valueBottomTop, mesh2.facesBottom)

eqn = DiffusionTerm()

if __name__ == '__main__':
    eqn.solve(var2)
    viewer = Viewer(var2)
    viewer.plot()
    input("finished")

