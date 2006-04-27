#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "inputTanh1D.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 8/10/05 {3:15:17 PM}
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
    >>> dx = L / nx
    >>> from fipy.meshes.grid1D import Grid1D
    >>> mesh = Grid1D(dx=dx, nx=nx)

and create the solution variable

    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(
    ...     name="phase field",
    ...     mesh=mesh,
    ...     value=1)

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
    ...     FixedValue(faces=mesh.getFacesRight(), value=1),
    ...     FixedValue(faces=mesh.getFacesLeft(), value=.5),
    ...     NthOrderBoundaryCondition(faces=mesh.getFacesLeft(), value=0, order=2),
    ...     NthOrderBoundaryCondition(faces=mesh.getFacesRight(), value=0, order=3))

Using

    >>> asq = 1.0
    >>> epsilon = 1
    >>> diffusionCoeff = 1

we create the Cahn-Hilliard equation:

    >>> faceVar = var.getArithmeticFaceValue()
    >>> doubleWellDerivative = asq * ( 1 - 6 * faceVar * (1 - faceVar))

    >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
    >>> from fipy.terms.transientTerm import TransientTerm
    >>> diffTerm2 = ImplicitDiffusionTerm(
    ...     coeff = (diffusionCoeff * doubleWellDerivative,))
    >>> diffTerm4 = ImplicitDiffusionTerm(coeff=(diffusionCoeff, epsilon**2))
    >>> eqch = TransientTerm() == diffTerm2 - diffTerm4

    >>> from fipy.solvers.linearLUSolver import LinearLUSolver
    >>> solver = LinearLUSolver(tolerance=1e-15, steps=100)

The solution to this 1D problem over an infinite domain is given by,

.. raw:: latex

   $$ \phi(x) = \frac{1}{1 + \exp{\left(-\frac{a}{\epsilon} x \right)}} $$

or

    >>> import fipy.tools.numerix as numerix
    >>> a = numerix.sqrt(asq)
    >>> answer = 1 / (1 + 
    ...     numerix.exp(-a * (mesh.getCellCenters()[:,0]) / epsilon))

If we are running interactively, we create a viewer to see the results

    >>> if __name__ == '__main__':
    ...     from fipy.viewers import make
    ...     viewer = make(vars=var,
    ...                   limits={'datamin': 0., 'datamax': 1.0})
    ...     viewer.plot()

We iterate the solution to equilibrium and, if we are running interactively, 
we update the display and output data about the progression of the solution

    >>> dexp=-5
    >>> for step in range(100):
    ...     dt = numerix.exp(dexp)
    ...     dt = min(10, dt)
    ...     dexp += 0.5
    ...     eqch.solve(var=var, boundaryConditions=BCs, solver=solver, dt=dt)
    ...     if __name__ == '__main__':
    ...         diff = abs(answer - numerix.array(var))
    ...         maxarg = numerix.argmax(diff)
    ...         print 'maximum error:',diff[maxarg]
    ...         print 'element id:',maxarg
    ...         print 'value at element ',maxarg,' is ',var[maxarg]
    ...         print 'solution value',answer[maxarg]
    ... 
    ...         viewer.plot()

We compare the analytical solution with the numerical result,

   >>> print var.allclose(answer, atol=1e-4)
   1

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    
    raw_input('finished')


