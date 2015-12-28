#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "grid3Dinput.py"
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

##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

"""
Test case for the Grid3D.

   >>> DiffusionTerm().solve(var)
   >>> DiffusionTerm().solve(var2)
   >>> a = numerix.array(var.globalValue)
   >>> b = numerix.array(var2.globalValue)
   >>> c = numerix.ravel(numerix.array((b, b, b)))
   >>> print numerix.allclose(a, c)
   1

"""

from fipy import CellVariable, Grid2D, Grid3D, DiffusionTerm, Viewer
from fipy.tools import numerix

nx = 10
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

nx = 10
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
    raw_input("finished")
