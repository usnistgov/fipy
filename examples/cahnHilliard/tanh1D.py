#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "inputTanh1D.py"
 #
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
 # ###################################################################
 ##

r"""

This example solves the Cahn-Hilliard equation given by,

.. math::

   \frac{\partial \phi}{\partial t} = \nabla \cdot D
    \nabla \left( \frac{\partial f}{\partial \phi}
        - \epsilon^2 \nabla^2 \phi \right)

where the free energy functional is given by,

.. math::

   f = \frac{a^2}{2} \phi^2 (1 - \phi)^2

The Cahn-Hilliard equation can be rewritten in the following form,

.. math::

   \frac{\partial \phi}{\partial t} = \nabla \cdot D
    \left( \frac{\partial^2 f}{\partial \phi^2} \nabla \phi
        - \epsilon^2 \nabla^3 \phi \right)

The above form of the equation makes the non-linearity part of the
diffusion coefficient for the first term on the RHS. This is the correct
way to express the equation to :term:`FiPy`.

We solve the problem on a 1D mesh

.. index:: Grid2D

>>> from fipy import CellVariable, Grid1D, NthOrderBoundaryCondition, DiffusionTerm, TransientTerm, LinearLUSolver, DefaultSolver, Viewer
>>> from fipy.tools import numerix

>>> L = 40.
>>> nx = 1000
>>> dx = L / nx
>>> mesh = Grid1D(dx=dx, nx=nx)

and create the solution variable

.. index:: CellVariable

>>> var = CellVariable(
...     name="phase field",
...     mesh=mesh,
...     value=1.)

The boundary conditions for this problem are

.. math::

   \left.
       \begin{aligned}
	   \phi &= \frac{1}{2} \\
	   \frac{\partial^2 \phi}{\partial x^2} &= 0
       \end{aligned}
   \right\} \qquad \text{on $x = 0$}

and

.. math::

   \left.
       \begin{aligned}
	   \phi &= 1 \\
	   \frac{\partial^2 \phi}{\partial x^2} &= 0
       \end{aligned}
   \right\} \qquad \text{on $x = L$}

or

.. index:: NthOrderBoundaryCondition

>>> BCs = (
...     NthOrderBoundaryCondition(faces=mesh.facesLeft, value=0, order=2),
...     NthOrderBoundaryCondition(faces=mesh.facesRight, value=0, order=2))

>>> var.constrain(1, mesh.facesRight)
>>> var.constrain(.5, mesh.facesLeft)

Using

>>> asq = 1.0
>>> epsilon = 1
>>> diffusionCoeff = 1

we create the Cahn-Hilliard equation:

>>> faceVar = var.arithmeticFaceValue
>>> freeEnergyDoubleDerivative = asq * ( 1 - 6 * faceVar * (1 - faceVar))

>>> diffTerm2 = DiffusionTerm(
...     coeff=diffusionCoeff * freeEnergyDoubleDerivative)
>>> diffTerm4 = DiffusionTerm(coeff=(diffusionCoeff, epsilon**2))
>>> eqch = TransientTerm() == diffTerm2 - diffTerm4

.. index:: LinearLUSolver, DefaultSolver

>>> import fipy.solvers.solver
>>> if fipy.solvers.solver == 'pysparse':
...     solver = LinearLUSolver(tolerance=1e-15, iterations=100)
... else:
...     solver = DefaultSolver()

The solution to this 1D problem over an infinite domain is given by,

.. math::

   \phi(x) = \frac{1}{1 + \exp{\left(-\frac{a}{\epsilon} x \right)}}

or

.. index:: sqrt, exp

>>> a = numerix.sqrt(asq)
>>> answer = 1 / (1 + numerix.exp(-a * (mesh.cellCenters[0]) / epsilon))

If we are running interactively, we create a viewer to see the results

.. index::
   module: fipy.viewers

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var, datamin=0., datamax=1.0)
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

    input('finished')
