#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/16/03 {3:23:47 PM}
 #                                last update: 10/7/04 {12:54:07 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
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

r"""

This example solves the steady-state convection-diffusion equation
given by:

.. raw:: latex

     $$ \nabla \cdot \left(D \nabla \phi + \vec{u} \phi \right) + S_c = 0 $$

with boundary conditions given by:

.. raw:: latex

   $$ \phi = \begin{cases}
   0& \text{at $x = 0$,} \\
   1& \text{at $x = L$,}
   \end{cases} $$ 

and coefficients

.. raw:: latex

   $D = 1$, $\vec{u} = (10, 0)$, and $S_c = 0$.     

The coefficients are represented by `diffCoeff` and `convCoeff` in
the Python code. The `SteadyConvectionDiffusionScEquation` object is
used to create the equation.  It needs to be passed a convection term
instantiator as follows:

   >>> from fipy.terms.exponentialConvectionTerm import ExponentialConvectionTerm
   >>> from fipy.equations.stdyConvDiffScEquation import SteadyConvectionDiffusionScEquation
   >>> eq = SteadyConvectionDiffusionScEquation(
   ...      var = var,
   ...      diffusionCoeff = diffCoeff,
   ...      convectionCoeff = convCoeff,
   ...      sourceCoeff = sourceCoeff,
   ...      solver = LinearCGSSolver(tolerance = 1.e-15, steps = 2000),
   ...      convectionScheme = ExponentialConvectionTerm,
   ...      boundaryConditions = boundaryConditions
   ...      )

More details of the benefits and drawbacks of each type of convection
term can be found in the numerical section of the manual. Essentially
the `ExponentialConvectionTerm` and `PowerLawConvectionTerm` will both
handle most types of convection diffusion cases with the
`PowerLawConvectionTerm` being more efficient.

The analytical solution test for this problem is given by:

   >>> axis = 0
   >>> x = mesh.getCellCenters()[:,axis]
   >>> AA = -sourceCoeff * x / convCoeff[axis]
   >>> BB = 1. + sourceCoeff * L / convCoeff[axis]
   >>> import Numeric
   >>> CC = 1. - Numeric.exp(-convCoeff[axis] * x / diffCoeff)
   >>> DD = 1. - Numeric.exp(-convCoeff[axis] * L / diffCoeff)
   >>> analyticalArray = AA + BB * CC / DD
   >>> Numeric.allclose(analyticalArray, var, rtol = 1e-10, atol = 1e-10)
   1
"""
__docformat__ = 'restructuredtext'
     
from fipy.meshes.grid2D import Grid2D
from fipy.equations.stdyConvDiffScEquation import SteadyConvectionDiffusionScEquation
from fipy.solvers.linearCGSSolver import LinearCGSSolver
from fipy.iterators.iterator import Iterator
from fipy.variables.cellVariable import CellVariable
from fipy.terms.scSourceTerm import ScSourceTerm
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.terms.exponentialConvectionTerm import ExponentialConvectionTerm
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux

valueLeft = 0.
valueRight = 1.
L = 10.
nx = 1000
ny = 1
diffCoeff = 1.
convCoeff = (10.,0.)
sourceCoeff = 0.

mesh = Grid2D(L / nx, L / ny, nx, ny)

var = CellVariable(
    name = "concentration",
    mesh = mesh,
    value = valueLeft)

boundaryConditions = (
    FixedValue(mesh.getFacesLeft(), valueLeft),
    FixedValue(mesh.getFacesRight(), valueRight),
    FixedFlux(mesh.getFacesTop(), 0.),
    FixedFlux(mesh.getFacesBottom(), 0.)
    )

        
eq = SteadyConvectionDiffusionScEquation(
    var = var,
    diffusionCoeff = diffCoeff,
    convectionCoeff = convCoeff,
    sourceCoeff = sourceCoeff,
    solver = LinearCGSSolver(tolerance = 1.e-15, steps = 2000),
    convectionScheme = ExponentialConvectionTerm,
    boundaryConditions = boundaryConditions
    )

it = Iterator((eq,))

it.timestep()

if __name__ == '__main__':
    viewer = Grid2DGistViewer(var)
##     print var
    viewer.plot()
    raw_input('finished')
