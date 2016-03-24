#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "coupled.py"
 #
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
 # ###################################################################
 ##

r"""Solve the biharmonic equation as a coupled pair of diffusion equations.

:term:`FiPy` has only first order time derivatives so equations such
as the biharmonic wave equation written as

.. math::

   \frac{\partial^4 v}{\partial x^4} + \frac{\partial^2 v}{\partial t^2} &= 0

cannot be represented as a single equation. We need to decompose the
biharmonic equation into two equations that are first order in time in
the following way,

.. math::

   \frac{\partial^2 v_0}{\partial x^2} + \frac{\partial v_1}{\partial t} &= 0 \\
   \frac{\partial^2 v_1}{\partial x^2} - \frac{\partial v_0}{\partial t} &= 0

Historically, :term:`FiPy` required systems of coupled equations to be
solved successively, "sweeping" the equations to convergence. As a
practical example, we use the following system

.. math::

   \frac{\partial v_0}{\partial t} &= 0.01 \nabla^2 v_0 - \nabla^2 v_1 \\
   \frac{\partial v_1}{\partial t} &= \nabla^2 v_0 + 0.01 \nabla^2 v_1

subject to the boundary conditions

.. math::
   :nowrap:

   \begin{align*}
   v_0|_{x=0} &= 0 & v_0|_{x=1} &= 1 \\
   v_1|_{x=0} &= 1 & v_1|_{x=1} &= 0
   \end{align*}

This system closely resembles the pure biharmonic equation, but has an
additional diffusion contribution to improve numerical stability.  The
example system is solved with the following block of code using
explicit coupling for the cross-coupled terms.

>>> from fipy import Grid1D, CellVariable, TransientTerm, DiffusionTerm, Viewer

>>> m = Grid1D(nx=100, Lx=1.)

>>> v0 = CellVariable(mesh=m, hasOld=True, value=0.5)
>>> v1 = CellVariable(mesh=m, hasOld=True, value=0.5)

>>> v0.constrain(0, m.facesLeft)
>>> v0.constrain(1, m.facesRight)

>>> v1.constrain(1, m.facesLeft)
>>> v1.constrain(0, m.facesRight)

>>> eq0 = TransientTerm() == DiffusionTerm(coeff=0.01) - v1.faceGrad.divergence
>>> eq1 = TransientTerm() == v0.faceGrad.divergence + DiffusionTerm(coeff=0.01)

>>> vi = Viewer((v0, v1))

>>> for t in range(100):
...     v0.updateOld()
...     v1.updateOld()
...     res0 = res1 = 1e100
...     while max(res0, res1) > 0.1:
...         res0 = eq0.sweep(var=v0, dt=1e-5)
...         res1 = eq1.sweep(var=v1, dt=1e-5)
...     if t % 10 == 0:
...         vi.plot()

The uncoupled method still works, but it can be advantageous to solve
the two equations simultaneously. In this case, by coupling the
equations, we can eliminate the explicit sources and dramatically
increase the time steps:

>>> v0.value = 0.5
>>> v1.value = 0.5

>>> eqn0 = TransientTerm(var=v0) == DiffusionTerm(0.01, var=v0) - DiffusionTerm(1, var=v1)
>>> eqn1 = TransientTerm(var=v1) == DiffusionTerm(1, var=v0) + DiffusionTerm(0.01, var=v1)

>>> eqn = eqn0 & eqn1

>>> for t in range(1):
...     v0.updateOld()
...     v1.updateOld()
...     eqn.solve(dt=1.e-3)
...     vi.plot()

It is also possible to pose the same equations in vector form:

>>> v = CellVariable(mesh=m, hasOld=True, value=[[0.5], [0.5]], elementshape=(2,))

>>> v.constrain([[0], [1]], m.facesLeft)
>>> v.constrain([[1], [0]], m.facesRight)

>>> eqn = TransientTerm([[1, 0],
...                      [0, 1]]) == DiffusionTerm([[[0.01, -1],
...                                                  [1, 0.01]]])

>>> vi = Viewer((v[0], v[1]))

>>> for t in range(1):
...     v.updateOld()
...     eqn.solve(var=v, dt=1.e-3)
...     vi.plot()

Whether you pose your problem in coupled or vector form should be dictated by
the underlying physics. If :math:`v_0` and :math:`v_1` represent the
concentrations of two conserved species, then it is natural to write two
seperate governing equations and to couple them. If they represent two
components of a vector field, then the vector formulation is obviously more
natural. FiPy will solve the same matrix system either way.
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')
