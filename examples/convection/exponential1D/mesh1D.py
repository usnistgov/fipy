#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh1D.py"
 #                                    created: 12/16/03 {3:23:47 PM}
 #                                last update: 6/24/08 {7:57:47 AM} 
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

   $D = 1$ and $\vec{u} = (10,)$,
   
or

    >>> diffCoeff = 1.
    >>> convCoeff = (10.,)
    
We define a 1D mesh

.. raw:: latex

   \IndexClass{Grid1D}

..

    >>> from fipy import *

    >>> L = 10.
    >>> nx = 10
    >>> mesh = Grid1D(dx=L / nx, nx=nx)

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
    >>> boundaryConditions = (
    ...     FixedValue(faces=mesh.getFacesLeft(), value=valueLeft),
    ...     FixedValue(faces=mesh.getFacesRight(), value=valueRight),
    ...     )

The solution variable is initialized to `valueLeft`:
    
.. raw:: latex

   \IndexClass{CellVariable}

..

    >>> var = CellVariable(mesh=mesh, name = "variable")

The equation is created with the `ImplicitDiffusionTerm` and
`ExponentialConvectionTerm`. The scheme used by the convection term
needs to calculate a Peclet number and thus the diffusion term
instance must be passed to the convection term.

.. raw:: latex

   \IndexClass{ImplicitDiffusionTerm}
   \IndexClass{ExponentialConvectionTerm}

..

   >>> eq = (ImplicitDiffusionTerm(coeff=diffCoeff)
   ...       + ExponentialConvectionTerm(coeff=convCoeff))
   
More details of the benefits and drawbacks of each type of convection
term can be found in 

.. raw:: latex

   Section~\ref{sec:NumericalSchemes} ``\nameref{sec:NumericalSchemes}''.
   
.. of the manual

Essentially, the `ExponentialConvectionTerm` and `PowerLawConvectionTerm` will
both handle most types of convection-diffusion cases, with the
`PowerLawConvectionTerm` being more efficient.

We solve the equation

   >>> eq.solve(var=var, boundaryConditions=boundaryConditions)
   
and test the solution against the analytical result

.. raw:: latex

   $$ \phi = \frac{1 - \exp(-u_x x / D)}{1 - \exp(-u_x L / D)} $$
   or
   \IndexFunction{exp}

..

    >>> axis = 0
    >>> x = mesh.getCellCenters()[axis]
    >>> CC = 1. - exp(-convCoeff[axis] * x / diffCoeff)
    >>> DD = 1. - exp(-convCoeff[axis] * L / diffCoeff)
    >>> analyticalArray = CC / DD
    >>> print var.allclose(analyticalArray)
    1
   
If the problem is run interactively, we can view the result:

.. raw:: latex

   \IndexModule{viewers}

..

    >>> if __name__ == '__main__':
    ...     viewer = Viewer(vars=var)
    ...     viewer.plot()
"""
__docformat__ = 'restructuredtext'
     
if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    
    raw_input('finished')
