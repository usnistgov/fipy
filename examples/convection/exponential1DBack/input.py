#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/16/03 {3:23:47 PM}
 #                                last update: 10/13/04 {12:09:07 PM} 
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

This example solves the steady-state convection-diffusion equation as described in::
    
    $ examples/diffusion/convection/exponential1D/input.py
    
but with

.. raw:: latex

    $ \\vec{u} = (-10, 0)$.

We test the solution against the analytical result:

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
convCoeff = (-10.,0.)
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
    viewer.plot()
    raw_input('finished')
