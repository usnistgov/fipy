##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

"""
Test case for the `Grid3D`. Diffusion problem with boundary conditions: 0 on front, 10 on back, and 5 on all other sides.

"""
from __future__ import unicode_literals

from fipy import input
from fipy import CellVariable, Grid3D, Viewer

nx = 10
ny = 10
nz = 10

dx = 1.
dy = 1.
dz = 1.

valueFront = 0.
valueBack = 10.
valueSides = 5.

mesh = Grid3D(dx = dx, dy = dy, dz = dz, nx = nx, ny = ny, nz = nz)

var = CellVariable(name = "variable",
                   mesh = mesh,
                   value = valueSides)

##viewer1 = Grid3DPyxViewer(var, zvalue = 1.0)
##viewer3 = Grid3DPyxViewer(var, zvalue = 3.0)
##viewer5 = Grid3DPyxViewer(var, zvalue = 5.0)
##viewer7 = Grid3DPyxViewer(var, zvalue = 7.0)
##viewer9 = Grid3DPyxViewer(var, zvalue = 9.0)

## viewer = Viewer(vars = var)

## viewer.plot()

var.constrain(valueSides, mesh.facesLeft)
var.constrain(valueSides, mesh.facesRight)
var.constrain(valueSides, mesh.facesTop)
var.constrain(valueSides, mesh.facesBottom)
var.constrain(valueFront, mesh.facesFront)
var.constrain(valueBack, mesh.facesBack)

## viewer.plot()

if __name__ == '__main__':
    ##viewer1.plot(resolution = 0.2, xlabel = "X values (Z value = 1)", minval = valueFront, maxval = valueBack)
    ##raw_input("press enter to continue")
    ##viewer3.plot(resolution = 0.2, xlabel = "X values (Z value = 3)", minval = valueFront, maxval = valueBack)
    ##raw_input("press enter to continue")
    ##viewer5.plot(resolution = 0.2, xlabel = "X values (Z value = 5)", minval = valueFront, maxval = valueBack)
    ##raw_input("press enter to continue")
    ##viewer7.plot(resolution = 0.2, xlabel = "X values (Z value = 7)", minval = valueFront, maxval = valueBack)
    ##raw_input("press enter to continue")
    ##viewer9.plot(resolution = 0.2, xlabel = "X values (Z value = 9)", minval = valueFront, maxval = valueBack)
    input("finished")
