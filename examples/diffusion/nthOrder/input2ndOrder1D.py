#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 3/7/05 {4:23:11 PM} 
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

r"""
In this problem, we demonstrate the use of the `NthOrderDiffusionEquation` class 
in the simple case of steady state 1D diffusion, which was introduced in
``examples/diffusion/steadyState/mesh1D/input.py``,
to solve

.. raw:: latex

   $$ \nabla \cdot (D \nabla \phi) = 0. $$

This examples shows that the `NthOrderDiffusionEquation` is equivalent to the `DiffusionEquation` when 

.. raw:: latex

   $n = 2$.
   
..

We create an appropriate 1D mesh:
    
    >>> nx = 10
    >>> ny = 1
    >>> dx = 1.
    >>> dy = 1.
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)

and initialize the solution variable to `valueLeft`:
    
    >>> valueLeft = 0.
    >>> valueRight = 1.
    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(name = "concentration",
    ...                    mesh = mesh,
    ...                    value = valueLeft)

The order 

.. raw:: latex

   $n$
   
of the `NthOrderDiffusionTerm` is determined by twice the number of
diffusion coefficients it is created with, so a single diffusion
coefficient ``(1.,)`` gives

.. raw:: latex

   $n = 2$.
   
We apply Dirichlet boundary conditions to the left and right and
Neumann boundary conditions to the top and bottom.
   
    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> BCs = (FixedValue(mesh.getFacesLeft(),valueLeft),
    ...        FixedValue(mesh.getFacesRight(),valueRight))

We iterate one timestep to equilibrium:

    >>> from fipy.terms.nthOrderDiffusionTerm import NthOrderDiffusionTerm
    >>> NthOrderDiffusionTerm(coeff = (1.,)).solve(var, boundaryConditions = BCs)
    
The result is tested against the expected linear diffusion profile:
    
    >>> Lx = nx * dx
    >>> x = mesh.getCellCenters()[:,0]
    >>> value = valueLeft + (valueRight - valueLeft) * x / Lx
    >>> var.allclose(value)
    1

If the problem is run interactively, we can view the result:
    
    >>> if __name__ == '__main__':
    ...     import fipy.viewers
    ...     viewer = fipy.viewers.make(vars = var)
    ...     viewer.plot()
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())
    raw_input("finished")
