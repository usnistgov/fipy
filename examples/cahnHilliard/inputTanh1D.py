#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "inputTanh1D.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 2/28/05 {4:31:46 PM}
 # Stolen from:
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

This example solves the Cahn-Hilliard equation given by:

.. raw:: latex

    $$ \frac{\partial \phi}{\partial t} = \nabla \cdot D 
	\nabla \left( \frac{\partial f}{\partial \phi} 
	    - \epsilon^2 \nabla^2 \phi \right) $$

where the free energy functional is given by,

.. raw:: latex

    $$ f = \frac{a^2}{2} \phi^2 (1 - \phi)^2. $$
    
We solve the problem on a 1D mesh

    >>> L = 40.
    >>> nx = 1000
    >>> ny = 1
    >>> dx = L / nx
    >>> dy = 1.
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx, dy, nx, ny)

and create the solution variable

    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(
    ...     name = "phase field",
    ...     mesh = mesh,
    ...     value = 1)
    >>> var.setValue(1, cells = mesh.getCells(lambda cell: cell.getCenter()[0] > L / 2))

The boundary conditions for this problem are

.. raw:: latex

   $$
   \left.
       \begin{aligned}
	   \phi &= \frac{1}{2} \\
	   \frac{\partial^3 \phi}{\partial x^3} &= 0
       \end{aligned}
   \right\} \qquad \text{on $x = 0$}
   $$

   and
   
   $$
   \left.
       \begin{aligned}
	   \phi &= 1 \\
	   \frac{\partial^2 \phi}{\partial x^2} &= 0
       \end{aligned}
   \right\} \qquad \text{on $x = L$}
   $$

or

    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> from fipy.boundaryConditions.nthOrderBoundaryCondition \
    ...     import NthOrderBoundaryCondition
    >>> BCs = (
    ...     FixedValue(mesh.getFacesRight(), 1),
    ...     FixedValue(mesh.getFacesLeft(), .5),
    ...     NthOrderBoundaryCondition(mesh.getFacesLeft(), 0, 2),
    ...     NthOrderBoundaryCondition(mesh.getFacesRight(), 0, 3))

Using

    >>> asq = 1.0
    >>> epsilon = 1
    >>> diffusionCoeff = 1

we create the Cahn-Hilliard equation:

    >>> faceVar = var.getArithmeticFaceValue()
    >>> doubleWellDerivative = asq * ( 1 - 6 * faceVar * (1 - faceVar))

    >>> from fipy.terms.nthOrderDiffusionTerm import NthOrderDiffusionTerm
    >>> from fipy.terms.transientTerm import TransientTerm
    >>> diffTerm2 = NthOrderDiffusionTerm(coeff = (diffusionCoeff * doubleWellDerivative,))
    >>> diffTerm4 = NthOrderDiffusionTerm(coeff = (diffusionCoeff, epsilon**2))
    >>> eqch = TransientTerm() == diffTerm2 - diffTerm4

    >>> from fipy.solvers.linearLUSolver import LinearLUSolver
    >>> solver = LinearLUSolver(tolerance = 1e-15, steps = 100)

The solution to this 1D problem over an infinite domain is given by,

.. raw:: latex

   $$ \phi(x) = \frac{1}{1 + \exp{\left(-\frac{a}{\epsilon} x \right)}} $$

or

    >>> import Numeric
    >>> a = Numeric.sqrt(asq)
    >>> answer = 1 / (1 + 
    ...     Numeric.exp(-a * (mesh.getCellCenters()[:,0]) / epsilon))

If we are running interactively, we create a viewer to see the results

    >>> if __name__ == '__main__':
    ...     from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
    ...     viewer = Grid2DGistViewer(var, minVal=0., maxVal=1.0, palette = 'rainbow.gp')
    ...     viewer.plot()

We iterate the solution to equilibrium and, if we are running interactively, 
we update the display and output data about the progression of the solution

    >>> dexp=-5
    >>> for step in range(100):
    ...     dt = Numeric.exp(dexp)
    ...     dt = min(10, dt)
    ...     dexp += 0.5
    ...     eqch.solve(var, boundaryConditions = BCs, solver = solver, dt = dt)
    ...     if __name__ == '__main__':
    ...         diff = abs(answer - Numeric.array(var))
    ...         maxarg = Numeric.argmax(diff)
    ...         print 'maximum error:',diff[maxarg]
    ...         print 'element id:',maxarg
    ...         print 'value at element ',maxarg,' is ',var[maxarg]
    ...         print 'solution value',answer[maxarg]
    ... 
    ...         viewer.plot()

We compare the analytical solution with the numerical result,

   >>> Numeric.allclose(var, answer, atol = 1e-4)
   1

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())
    
    raw_input('finished')


