#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "ttri2Dinput.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 12/7/04 {10:26:55 AM} 
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

"""
This input file again solves a steady 1D diffusion problem as in
``examples/diffusion/steadyState/mesh1D/input.py``,
the difference being that the mesh is two dimensional:

    >>> nx = 20
    >>> ny = 20
    >>> dx = 1.
    >>> dy = 1.
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)

We create a `CellVariable` and initialize it to `valueLeft`:
    
    >>> valueLeft = 0.
    >>> valueRight = 1.
    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(name = "solution variable",
    ...                    mesh = mesh,
    ...                    value = valueLeft)

We create a diffusion equation.  This is solved by default with an
iterative conjugate gradient solver.  We apply Dirichlet boundary
conditions to the left and right.  Neumann boundary conditions are
automatically applied to the top and bottom.
    
    >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> ImplicitDiffusionTerm().solve(var = var, 
    ...                               boundaryConditions = (FixedValue(mesh.getFacesLeft(),valueLeft),
    ...                                                     FixedValue(mesh.getFacesRight(),valueRight)))

The result is again tested against the expected linear composition profile:

    >>> Lx = nx * dx
    >>> x = mesh.getCellCenters()[:,0]
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    >>> var.allclose(analyticalArray, rtol = 1e-10, atol = 1e-10)
    1
    
If the problem is run interactively, we can view the result:
    
    >>> if __name__ == '__main__':
    ...     from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
    ...     viewer = Grid2DGistViewer(var)
    ...     viewer.plot()
"""

__docformat__ = 'restructuredtext'

##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())
    raw_input("finished")
