#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/16/03 {3:23:47 PM}
 #                                last update: 10/27/04 {9:50:47 AM} 
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

Like::
    
    $ examples/diffusion/convection/exponential1D/input.py
    
this example solves a steady-state convection-diffusion equation, but adds a constant source, 

.. raw:: latex

     $S_0 = 1$, such that

     $$ \nabla \cdot \left(D \nabla \phi + \vec{u} \phi \right) + S_0 = 0. $$

Here, the axes are reversed

    >>> nx = 1
    >>> ny = 1000

and

.. raw:: latex

    $ \\vec{u} = (0, 10)$
    
such that

    >>> diffCoeff = 1.
    >>> convCoeff = (0., 10.)
    >>> sourceCoeff = 1.

We define a 1D mesh

    >>> L = 10.
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(L / nx, L / ny, nx, ny)

and impose the boundary conditions

.. raw:: latex

   $$ \phi = \begin{cases}
   0& \text{at $y = 0$,} \\
   1& \text{at $y = L$,}
   \end{cases} $$ 

or

    >>> valueBottom = 0.
    >>> valueTop = 1.
    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> from fipy.boundaryConditions.fixedFlux import FixedFlux
    >>> boundaryConditions = (
    ...     FixedValue(mesh.getFacesTop(), valueTop),
    ...     FixedValue(mesh.getFacesBottom(), valueBottom),
    ...     FixedFlux(mesh.getFacesRight(), 0.),
    ...     FixedFlux(mesh.getFacesLeft(), 0.)
    ...     )

The solution variable is initialized to `valueBottom`:
    
    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(
    ...     name = "concentration",
    ...     mesh = mesh,
    ...     value = valueBottom)

We define the convection-diffusion equation with source

    >>> from fipy.terms.exponentialConvectionTerm import ExponentialConvectionTerm
    >>> from fipy.solvers.linearLUSolver import LinearLUSolver
    >>> from fipy.equations.stdyConvDiffScEquation import SteadyConvectionDiffusionScEquation
    >>> eq = SteadyConvectionDiffusionScEquation(
    ...      var = var,
    ...      diffusionCoeff = diffCoeff,
    ...      convectionCoeff = convCoeff,
    ...      sourceCoeff = sourceCoeff,
    ...      solver = LinearLUSolver(tolerance = 1.e-15),
    ...      convectionScheme = ExponentialConvectionTerm,
    ...      boundaryConditions = boundaryConditions
    ...      )
    
iterate one implicit timestep to equilibrium

    >>> from fipy.iterators.iterator import Iterator
    >>> it = Iterator((eq,))
    >>> it.timestep()

and test the solution against the analytical result:
    
.. raw:: latex

   $$ \phi = -\frac{S_0 y}{u_y} 
   + \left(1 + \frac{S_0 y}{u_y}\right)\frac{1 - \exp(-u_y y / D)}{1 - \exp(-u_y L / D)} $$

or

    >>> axis = 1
    >>> y = mesh.getCellCenters()[:,axis]
    >>> AA = -sourceCoeff * y / convCoeff[axis]
    >>> BB = 1. + sourceCoeff * L / convCoeff[axis]
    >>> import Numeric
    >>> CC = 1. - Numeric.exp(-convCoeff[axis] * y / diffCoeff)
    >>> DD = 1. - Numeric.exp(-convCoeff[axis] * L / diffCoeff)
    >>> analyticalArray = AA + BB * CC / DD
    >>> Numeric.allclose(analyticalArray, var, rtol = 1e-4, atol = 1e-4) 
    1
         
If the problem is run interactively, we can view the result:

    >>> if __name__ == '__main__':
    ...     from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
    ...     viewer = Grid2DGistViewer(var)
    ...     viewer.plot()

"""
__docformat__ = 'restructuredtext'

## from fipy.solvers.linearCGSSolver import LinearCGSSolver
## solver = LinearCGSSolver(tolerance = 1.e-15, steps = 2000),


if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())

    raw_input('finished')
