#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 10/6/04 {2:23:02 PM} 
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

r"""

This example solves the following equation:

.. raw:: latex

    $$ \frac{\partial^4 \phi}{\partial x^4} = 0 $$

with the following boundary conditions:

.. raw:: latex

    \begin{align*}
    \phi &= \alpha_1 && \text{at $x = 0$} \\
    \frac{\partial \phi}{\partial x} &= \alpha_2 && \text{at $x = L$} \\
    \frac{\partial^2 \phi}{\partial x^2} &= \alpha_3 && \text{at $x = 0$} \\
    \frac{\partial^3 \phi}{\partial x^3} &= \alpha_4 && \text{at $x = L$.}
    \end{align*}

The analytical solution is:

.. raw:: latex

    $$ \phi = \frac{ \alpha_4 }{6} x^3 + \frac{ \alpha_3 }{2} x^2 
    + \left( \alpha_2 - \frac{ \alpha_4 }{2} L^2  - \alpha_3 L \right) x + \alpha_1 $$

Do a timestep to steady state

   >>> it.timestep()

Test the solution:

   >>> x = mesh.getCellCenters()[:,0]
   >>> answer = alpha4 / 6. * x**3 + alpha3 / 2. * x**2 
   >>> answer += (alpha2 - alpha4 / 2. * L**2 - alpha3 * L) * x + alpha1
   >>> Numeric.allclose(answer, var, atol = 1e-10)
   1

"""
__docformat__ = 'restructuredtext'

import Numeric

from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.boundaryConditions.nthOrderBoundaryCondition import NthOrderBoundaryCondition
from fipy.equations.nthOrderDiffusionEquation import NthOrderDiffusionEquation
from fipy.iterators.iterator import Iterator
from fipy.meshes.grid2D import Grid2D
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.variables.cellVariable import CellVariable
from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
from fipy.solvers.linearLUSolver import LinearLUSolver

alpha1 = 2.
alpha2 = 1.
alpha3 = 4.
alpha4 = -3.

L = 1000.

nx = 1000
ny = 1

dx = L / nx
dy = 1.

mesh = Grid2D(dx, dy, nx, ny)
    
var = CellVariable(
    name = "concentration",
    mesh = mesh,
    value = 0.)

eq = NthOrderDiffusionEquation(
    var,
    transientCoeff = 0.0, 
    diffusionCoeff = (-1., 1.),
    solver = LinearLUSolver(tolerance = 1e-11),
    boundaryConditions=(
        FixedValue(mesh.getFacesLeft(), alpha1),
        FixedFlux(mesh.getFacesRight(), alpha2),
        NthOrderBoundaryCondition(mesh.getFacesLeft(), alpha3, 2),
        NthOrderBoundaryCondition(mesh.getFacesRight(), alpha4, 3)))

it = Iterator((eq,))

if __name__ == '__main__':
    viewer = Grid2DGistViewer(var)
    it.timestep()
    viewer.plot()
##     print var
    raw_input('finished')
