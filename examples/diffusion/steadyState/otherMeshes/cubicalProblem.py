#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "ttri2Dinput.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 9/3/04 {10:37:39 PM} 
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

##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

"""
Test case for the Grid3D. Diffusion problem with boundary conditions: 0 on front, 10 on back, and 5 on all other sides.
   
"""

from fipy.meshes.numMesh.grid3D import Grid3D
import Numeric
from fipy.equations.diffusionEquation import DiffusionEquation
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.iterators.iterator import Iterator
from fipy.variables.cellVariable import CellVariable
from fipy.viewers.pyxviewer import Grid3DPyxViewer

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

viewer1 = Grid3DPyxViewer(var, zvalue = 1.0)
viewer3 = Grid3DPyxViewer(var, zvalue = 3.0)
viewer5 = Grid3DPyxViewer(var, zvalue = 5.0)
viewer7 = Grid3DPyxViewer(var, zvalue = 7.0)
viewer9 = Grid3DPyxViewer(var, zvalue = 9.0)

eq = DiffusionEquation(var,
                       transientCoeff = 0., 
                       diffusionCoeff = 1.,
                       solver = LinearPCGSolver(tolerance = 1.e-15, 
                                                steps = 1000
                                                ),
                       boundaryConditions = (FixedValue(mesh.getFacesLeft(),valueSides),
                                             FixedValue(mesh.getFacesRight(),valueSides),
                                             FixedValue(mesh.getFacesTop(),valueSides),
                                             FixedValue(mesh.getFacesBottom(),valueSides),
                                             FixedValue(mesh.getFacesFront(),valueFront),
                                             FixedValue(mesh.getFacesBack(),valueBack),
                                             )
                       )

it = Iterator((eq,))

it.timestep()

if __name__ == '__main__':
    viewer1.plot(resolution = 0.2, xlabel = "X values (Z value = 1)", minval = valueFront, maxval = valueBack)
    raw_input("press enter to continue")
    viewer3.plot(resolution = 0.2, xlabel = "X values (Z value = 3)", minval = valueFront, maxval = valueBack)
    raw_input("press enter to continue")
    viewer5.plot(resolution = 0.2, xlabel = "X values (Z value = 5)", minval = valueFront, maxval = valueBack)
    raw_input("press enter to continue")
    viewer7.plot(resolution = 0.2, xlabel = "X values (Z value = 7)", minval = valueFront, maxval = valueBack)
    raw_input("press enter to continue")
    viewer9.plot(resolution = 0.2, xlabel = "X values (Z value = 9)", minval = valueFront, maxval = valueBack)
    raw_input("finished")
