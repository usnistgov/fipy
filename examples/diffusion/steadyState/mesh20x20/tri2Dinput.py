#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "ttri2Dinput.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 9/3/04 {10:41:36 PM} 
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

This input file again solves a 1D diffusion problem as in
`./examples/diffusion/steadyState/mesh1D/input.py`. The difference
being that the mesh size is given by 

    >>> nx = 20
    >>> ny = 20

The result is again tested in the same way:

    >>> Lx = nx * dx
    >>> x = mesh.getCellCenters()[:,0]
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    >>> import Numeric
    >>> Numeric.allclose(Numeric.array(var), analyticalArray, rtol = 1e-10, atol = 1e-10)
    1

"""

from fipy.meshes.numMesh.tri2D import Tri2D
from fipy.equations.diffusionEquation import DiffusionEquation
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.iterators.iterator import Iterator
from fipy.variables.cellVariable import CellVariable
from fipy.viewers.pyxviewer import PyxViewer

nx = 20
ny = 20

dx = 1.
dy = 1.

valueLeft = 0.
valueRight = 1.

mesh = Tri2D(dx = dx, dy = dy, nx = nx, ny = ny)

var = CellVariable(name = "solution variable",
                   mesh = mesh,
                   value = valueLeft)

viewer = PyxViewer(var)

eq = DiffusionEquation(var,
                       transientCoeff = 0., 
                       diffusionCoeff = 1.,
                       solver = LinearPCGSolver(tolerance = 1.e-15, 
                                                steps = 1000
                                                ),
                       boundaryConditions = (FixedValue(mesh.getFacesLeft(),valueLeft),
                                             FixedValue(mesh.getFacesRight(),valueRight),
                                             FixedFlux(mesh.getFacesTop(),0.),
                                             FixedFlux(mesh.getFacesBottom(),0.)
                                             )
                       )

it = Iterator((eq,))

it.timestep()

if __name__ == '__main__':
    viewer.plot(resolution = 0.2)
    raw_input("finished")
