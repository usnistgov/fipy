#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "ttri2Dinput.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:02:29 PM} 
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

##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

"""
Test case for the Grid3D.

   >>> a = Numeric.array(var)
   >>> b = Numeric.array(var2)
   >>> c = Numeric.ravel(Numeric.array((b, b, b)))
   >>> print Numeric.allclose(a, c)
   1
   
"""

from fipy.meshes.numMesh.grid3D import Grid3D
from fipy.meshes.grid2D import Grid2D
import Numeric
from fipy.equations.diffusionEquation import DiffusionEquation
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.iterators.iterator import Iterator
from fipy.variables.cellVariable import CellVariable
from fipy.viewers.pyxviewer import Grid3DPyxViewer
from fipy.viewers.pyxviewer import Grid2DPyxViewer

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

viewer = Grid3DPyxViewer(var, zvalue = 1.0)

eq = DiffusionEquation(var,
                       transientCoeff = 0., 
                       diffusionCoeff = 1.,
                       solver = LinearPCGSolver(tolerance = 1.e-15, 
                                                steps = 1000
                                                ),
                       boundaryConditions = (FixedValue(mesh.getFacesLeft(),valueLeftRight),
                                             FixedValue(mesh.getFacesRight(),valueLeftRight),
                                             FixedValue(mesh.getFacesTop(),valueBottomTop),
                                             FixedValue(mesh.getFacesBottom(),valueBottomTop),
                                             FixedFlux(mesh.getFacesFront(), 0.),
                                             FixedFlux(mesh.getFacesBack(), 0.),
                                             )
                       )

it = Iterator((eq,))

it.timestep()

#do the 2D problem for comparison

nx = 10
ny = 5

dx = 1.
dy = 1.

mesh2 = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)

var2 = CellVariable(name = "solution variable 2D",
                   mesh = mesh2,
                   value = valueBottomTop)

viewer2 = Grid2DPyxViewer(var2)

eq2 = DiffusionEquation(var2,
                       transientCoeff = 0., 
                       diffusionCoeff = 1.,
                       solver = LinearPCGSolver(tolerance = 1.e-15, 
                                                steps = 1000
                                                ),
                       boundaryConditions = (FixedValue(mesh2.getFacesLeft(),valueLeftRight),
                                             FixedValue(mesh2.getFacesRight(),valueLeftRight),
                                             FixedValue(mesh2.getFacesTop(),valueBottomTop),
                                             FixedValue(mesh2.getFacesBottom(),valueBottomTop)
                                             )
                       )

it2 = Iterator((eq2,))

it2.timestep()

if __name__ == '__main__':
    viewer.plot(resolution = 0.1)
    raw_input("finished")
