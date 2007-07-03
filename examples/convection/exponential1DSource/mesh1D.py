#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh1D.py"
 #                                    created: 12/16/03 {3:23:47 PM}
 #                                last update: 3/29/07 {11:50:05 AM} 
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

Like ``examples/diffusion/convection/exponential1D/mesh1D.py``
this example solves a steady-state convection-diffusion equation, but adds a constant source, 

.. raw:: latex

     $S_0 = 1$, such that

     $$ \nabla \cdot \left(D \nabla \phi + \vec{u} \phi \right) + S_0 = 0. $$

..

    >>> diffCoeff = 1.
    >>> convCoeff = (10.,)
    >>> sourceCoeff = 1.

We define a 1D mesh

.. raw:: latex

   \IndexClass{Grid1D}

..

    >>> nx = 1000
    >>> L = 10.
    >>> from fipy.meshes.grid1D import Grid1D
    >>> mesh = Grid1D(dx=L / 1000, nx=nx)

and impose the boundary conditions

.. raw:: latex

   $$ \phi = \begin{cases}
   0& \text{at $x = 0$,} \\
   1& \text{at $x = L$,}
   \end{cases} $$ 
   or
   \IndexClass{FixedValue}
   
..

    >>> valueLeft = 0.
    >>> valueRight = 1.
    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> boundaryConditions = (
    ...     FixedValue(faces=mesh.getFacesRight(), value=valueRight),
    ...     FixedValue(faces=mesh.getFacesLeft(), value=valueLeft),
    ...     )

The solution variable is initialized to `valueLeft`:
    
.. raw:: latex

   \IndexClass{CellVariable}

..

    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(name="variable", mesh=mesh)


We define the convection-diffusion equation with source

.. raw:: latex

   \IndexClass{ImplicitDiffusionTerm}
   \IndexClass{ExponentialConvectionTerm}

..

    >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
    >>> from fipy.terms.exponentialConvectionTerm \
    ...     import ExponentialConvectionTerm
    >>> eq = (ImplicitDiffusionTerm(coeff=diffCoeff)
    ...       + ExponentialConvectionTerm(coeff=convCoeff)
    ...       + sourceCoeff)
    
.. raw:: latex

   \IndexClass{LinearLUSolver}

..
    
    >>> from fipy.solvers.linearLUSolver import LinearLUSolver
    >>> eq.solve(var = var, 
    ...          boundaryConditions = boundaryConditions,
    ...          solver = LinearLUSolver(tolerance = 1.e-15))
    
and test the solution against the analytical result:
    
.. raw:: latex

   $$ \phi = -\frac{S_0 x}{u_x} 
   + \left(1 + \frac{S_0 x}{u_x}\right)\frac{1 - \exp(-u_x x / D)}{1 - \exp(-u_x L / D)} $$
   or
   \IndexModule{numerix}
   \IndexFunction{exp}

..

    >>> axis = 0
    >>> x = mesh.getCellCenters()[axis]
    >>> AA = -sourceCoeff * x / convCoeff[axis]
    >>> BB = 1. + sourceCoeff * L / convCoeff[axis]
    >>> from fipy.tools.numerix import exp
    >>> CC = 1. - exp(-convCoeff[axis] * x / diffCoeff)
    >>> DD = 1. - exp(-convCoeff[axis] * L / diffCoeff)
    >>> analyticalArray = AA + BB * CC / DD
    >>> print var.allclose(analyticalArray, rtol=1e-4, atol=1e-4)
    1
         
If the problem is run interactively, we can view the result:

.. raw:: latex

   \IndexModule{viewers}

..

    >>> if __name__ == '__main__':
    ...     from fipy.viewers import make
    ...     viewer = make(vars=var)
    ...     viewer.plot()

"""
__docformat__ = 'restructuredtext'

## from fipy.solvers.linearCGSSolver import LinearCGSSolver
## solver = LinearCGSSolver(tolerance = 1.e-15, steps = 2000),


if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input('finished')
