#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "input.py"
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

r"""Solve a fourth-order diffusion problem.

This example uses the :class:`~fipy.terms.diffusionTerm.DiffusionTerm` class to solve the equation

.. math::

   \frac{\partial^4 \phi}{\partial x^4} = 0

on a 1D mesh of length

>>> L = 1000.

We create an appropriate mesh

.. index:: Grid1D

>>> from fipy import CellVariable, Grid1D, NthOrderBoundaryCondition, DiffusionTerm, Viewer, GeneralSolver

>>> nx = 500
>>> dx = L / nx
>>> mesh = Grid1D(dx=dx, nx=nx)

and initialize the solution variable to 0

.. index:: CellVariable

>>> var = CellVariable(mesh=mesh, name='solution variable')

For this problem, we impose the boundary conditions:

.. math::

   \phi &= \alpha_1 \quad \text{at $x = 0$} \\
   \frac{\partial \phi}{\partial x} &= \alpha_2 \quad \text{at $x = L$} \\
   \frac{\partial^2 \phi}{\partial x^2} &= \alpha_3 \quad \text{at $x = 0$} \\
   \frac{\partial^3 \phi}{\partial x^3} &= \alpha_4 \quad \text{at $x = L$.}

or

>>> alpha1 = 2.
>>> alpha2 = 1.
>>> alpha3 = 4.
>>> alpha4 = -3.

.. index:: NthOrderBoundaryCondition

>>> BCs = (NthOrderBoundaryCondition(faces=mesh.facesLeft, value=alpha3, order=2),
...        NthOrderBoundaryCondition(faces=mesh.facesRight, value=alpha4, order=3))
>>> var.faceGrad.constrain([alpha2], mesh.facesRight)
>>> var.constrain(alpha1, mesh.facesLeft)

We initialize the steady-state equation

>>> eq = DiffusionTerm(coeff=(1, 1)) == 0

and use the :class:`~fipy.solvers.pysparse.linearLUSolver.LinearLUSolver` for stability.

We perform one implicit timestep to achieve steady state

>>> eq.solve(var=var,
...          boundaryConditions=BCs,
...          solver=GeneralSolver())

The analytical solution is:

.. math::

   \phi = \frac{ \alpha_4 }{6} x^3 + \frac{ \alpha_3 }{2} x^2
   + \left( \alpha_2 - \frac{ \alpha_4 }{2} L^2  - \alpha_3 L \right) x + \alpha_1

or

>>> analytical = CellVariable(mesh=mesh, name='analytical value')
>>> x = mesh.cellCenters[0]
>>> analytical.setValue(alpha4 / 6. * x**3 + alpha3 / 2. * x**2 + \
...                     (alpha2 - alpha4 / 2. * L**2 - alpha3 * L) * x + alpha1)

>>> print var.allclose(analytical, rtol=1e-4)
1

If the problem is run interactively, we can view the result:

.. index::
   module: fipy.viewers

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=(var, analytical))
...     viewer.plot()

.. image:: input4thOrder1D.*
   :width: 90%
   :align: center
   :alt: solution to biharmonic equation

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    input('finished')
