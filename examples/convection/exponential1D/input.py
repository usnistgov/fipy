#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/16/03 {3:23:47 PM}
 #                                last update: 2/18/05 {3:05:46 PM} 
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

This example solves the steady-state convection-diffusion equation
given by:

.. raw:: latex

   $$ \nabla \cdot \left(D \nabla \phi + \vec{u} \phi \right) = 0 $$

with coefficients

.. raw:: latex

   $D = 1$ and $\vec{u} = (10, 0)$,
   
or

    >>> diffCoeff = 1.
    >>> convCoeff = (10.,0.)
    
We define a 1D mesh

    >>> L = 10.
    >>> nx = 1000
    >>> ny = 1
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(L / nx, L / ny, nx, ny)

and impose the boundary conditions

.. raw:: latex

   $$ \phi = \begin{cases}
   0& \text{at $x = 0$,} \\
   1& \text{at $x = L$,}
   \end{cases} $$ 

or

    >>> valueLeft = 0.
    >>> valueRight = 1.
    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> boundaryConditions = (
    ...     FixedValue(mesh.getFacesLeft(), valueLeft),
    ...     FixedValue(mesh.getFacesRight(), valueRight),
    ...     )

The solution variable is initialized to `valueLeft`:
    
    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(
    ...     name = "concentration",
    ...     mesh = mesh,
    ...     value = valueLeft)

The `SteadyConvectionDiffusionScEquation` object is
used to create the equation.  It needs to be passed a convection term
instantiator as follows:

   >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
   >>> diffTerm = ImplicitDiffusionTerm(coeff = diffCoeff)
   
   >>> from fipy.terms.exponentialConvectionTerm import ExponentialConvectionTerm
   >>> eq = diffTerm + ExponentialConvectionTerm(coeff = convCoeff, diffusionTerm = diffTerm)
   
More details of the benefits and drawbacks of each type of convection
term can be found in the numerical section of the manual. Essentially
the `ExponentialConvectionTerm` and `PowerLawConvectionTerm` will both
handle most types of convection diffusion cases with the
`PowerLawConvectionTerm` being more efficient.

We solve the equation

   >>> from fipy.solvers.linearCGSSolver import LinearCGSSolver
   >>> eq.solve(var = var, 
   ...          solver = LinearCGSSolver(tolerance = 1.e-15, steps = 2000), 
   ...          boundaryConditions = boundaryConditions)
   
and test the solution against the analytical result

.. raw:: latex

   $$ \phi = \frac{1 - \exp(-u_x x / D)}{1 - \exp(-u_x L / D)} $$

or 

    >>> axis = 0
    >>> x = mesh.getCellCenters()[:,axis]
    >>> import Numeric
    >>> CC = 1. - Numeric.exp(-convCoeff[axis] * x / diffCoeff)
    >>> DD = 1. - Numeric.exp(-convCoeff[axis] * L / diffCoeff)
    >>> analyticalArray = CC / DD
    >>> var.allclose(analyticalArray, rtol = 1e-10, atol = 1e-10)
    1
   
If the problem is run interactively, we can view the result:

    >>> if __name__ == '__main__':
    ...     from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
    ...     viewer = Grid2DGistViewer(var)
    ...     viewer.plot()
"""
__docformat__ = 'restructuredtext'
     
if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())
    
    raw_input('finished')
