#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
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

This input file again solves a 1D diffusion problem as in::
    
    $ examples/diffusion/steadyState/mesh1D/input.py
    
The difference being that the mesh is two dimensional.

The result is again tested in the same way:

    >>> DiffusionTerm().solve(var)
    >>> Lx = nx * dx
    >>> x = mesh.cellCenters[0]
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    >>> print var.allclose(analyticalArray, rtol = 1e-9)
    1

"""

__docformat__ = 'restructuredtext'

from fipy import CellVariable, Grid2D, DiffusionTerm, Viewer

nx = 50
ny = 50

dx = 1.

valueLeft = 0.
valueRight = 1.

mesh = Grid2D(dx = dx, nx = nx, ny = ny)

var = CellVariable(name = "solution variable",
                   mesh = mesh,
                   value = valueLeft)

var.constrain(valueLeft, mesh.facesLeft)
var.constrain(valueRight, mesh.facesRight)
    
if __name__ == '__main__':
    DiffusionTerm().solve(var)
    
    viewer = Viewer(vars=var, datamin=0., datamax=1.)
    viewer.plot()
    raw_input("finished")
