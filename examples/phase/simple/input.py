#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 9/10/04 {2:27:45 PM} 
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

To run this example from the base fipy directory type
`./examples/phase/simple/input.py` at the command line.  A gist viewer
object should appear and the word `finished` in the terminal.

This example takes the user through assembling a simple problem with
FiPy.  It describes a steady 1D phase field problem with fixed value
boundary conditions such that,

.. raw:: latex

   $$ \frac{1}{M_\phi}\frac{\partial \phi}{\partial t} 
   = -\frac{\partial f}{\partial \phi} 
   + \kappa_\phi \nabla^2\phi$$.

For solidification problems, the Helmholtz free energy is frequently
given by

.. raw:: latex

   $$f(\phi,T) = \frac{W}{2}g(\phi) + \frac{L(T-T_M)}{T_M}p(\phi)$$.
   
One possible choice for the double-well function is

.. raw:: latex

   $$g(\phi) = \phi^2(1 - \phi)^2$$
   
and for the interpolation function is

.. raw:: latex

   $$p(\phi) = \phi^3(6\phi^2 - 15\phi + 10)$$

Our initial condition is a step function

.. raw:: latex

   $$\phi = 
   \begin{cases}
       1 & \text{for $x \le L/2$} \\
       0 & \text{for $x > L/2$}
   \end{cases}
   \text{at $t = 0$}$$

and we impose the boundary conditions,

.. raw:: latex

   $$ \phi = 
   \begin{cases}
       1 & \text{at $x = 0$} \\
       0 & \text{at $x = L$}
   \end{cases}$$

and parameter values,

.. raw:: latex

   $$ 1/M_\phi = 0 \;\; \\text{and} \;\; D = 1 $$

To test the solution, the analytical result is required. The x
coordinates from the mesh are gathered and the length of the domain,
`Lx`, is calculated.  An array, `analyticalArray`, is calculated to
compare with the numerical result,

    >>> for i in range(20):
    ...     eq.solve(phase)
    >>> x = mesh.getCellCenters()[:,0]
    >>> import Numeric
    >>> analyticalArray = 0.5*(1 - Numeric.tanh((x - L/2)/(2*Numeric.sqrt(kappa/W))))

Finally the analytical and numerical results are compared with a
tolerance of `1e-4`.

    >>> phase.allclose(analyticalArray, rtol = 1e-4, atol = 1e-4)
    1

"""
__docformat__ = 'restructuredtext'

import Numeric

from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fipy.terms.dependentSourceTerm import DependentSourceTerm
from fipy.meshes.grid2D import Grid2D
from fipy.variables.cellVariable import CellVariable
from fipy.viewers.gist1DViewer import Gist1DViewer

nx = 400
ny = 1

valueLeft = 1.
valueRight = 0.
timeStepDuration = 1. 

kappa = 0.0025
W = 1.
Lv = 1.
Tm = 1.
T = Tm
enthalpy = Lv * (T - Tm) / Tm

L = 1.

dx = L / nx
dy = 1.

mesh = Grid2D(dx, dy, nx, ny)
    
phase = CellVariable(
    name = "phase",
    mesh = mesh,
    value = valueLeft)
    
setCells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] > L/2)
phase.setValue(valueLeft)
phase.setValue(valueRight,setCells)

diffusionTerm = ImplicitDiffusionTerm(diffCoeff = kappa)

mPhi = -(30. * phase * (1. - phase) * enthalpy + (1. - 2 * phase) * W)

spTerm = DependentSourceTerm(sourceCoeff = mPhi * (phase - (mPhi < 0.)))

sourceCoeff = (mPhi > 0.) * mPhi * phase

eq = spTerm - sourceCoeff - diffusionTerm
    
if __name__ == '__main__':
    for i in range(20):
        eq.solve(phase)
    viewer = Gist1DViewer(vars = (phase,))    
    viewer.plot()
    raw_input('finished')
