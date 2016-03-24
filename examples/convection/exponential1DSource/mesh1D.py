#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "mesh1D.py"
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

r"""Solve the steady-state convection-diffusion equation with a constant source.

Like :mod:`examples.convection.exponential1D.mesh1D`
this example solves a steady-state convection-diffusion equation, but adds a constant source,
:math:`S_0 = 1`, such that

.. math::

   \nabla \cdot \left(D \nabla \phi + \vec{u} \phi \right) + S_0 = 0.

>>> diffCoeff = 1.
>>> convCoeff = (10.,)
>>> sourceCoeff = 1.

We define a 1D mesh

.. index:: Grid1D

>>> from fipy import CellVariable, Grid1D, DiffusionTerm, ExponentialConvectionTerm, DefaultAsymmetricSolver, Viewer
>>> from fipy.tools import numerix

>>> nx = 1000
>>> L = 10.
>>> mesh = Grid1D(dx=L / 1000, nx=nx)

>>> valueLeft = 0.
>>> valueRight = 1.

The solution variable is initialized to ``valueLeft``:

.. index:: CellVariable

>>> var = CellVariable(name="variable", mesh=mesh)

and impose the boundary conditions

.. math::

   \phi = \begin{cases}
   0& \text{at $x = 0$,} \\
   1& \text{at $x = L$,}
   \end{cases}

with

>>> var.constrain(valueLeft, mesh.facesLeft)
>>> var.constrain(valueRight, mesh.facesRight)

We define the convection-diffusion equation with source

>>> eq = (DiffusionTerm(coeff=diffCoeff)
...       + ExponentialConvectionTerm(coeff=convCoeff)
...       + sourceCoeff)

.. .. index:: DefaultAsymmetricSolver

>>> eq.solve(var=var,
...          solver=DefaultAsymmetricSolver(tolerance=1.e-15, iterations=10000))

and test the solution against the analytical result:

.. math::

   \phi = -\frac{S_0 x}{u_x}
   + \left(1 + \frac{S_0 x}{u_x}\right)\frac{1 - \exp(-u_x x / D)}{1 - \exp(-u_x L / D)}

or

.. index:: exp

>>> axis = 0
>>> x = mesh.cellCenters[axis]
>>> AA = -sourceCoeff * x / convCoeff[axis]
>>> BB = 1. + sourceCoeff * L / convCoeff[axis]
>>> CC = 1. - numerix.exp(-convCoeff[axis] * x / diffCoeff)
>>> DD = 1. - numerix.exp(-convCoeff[axis] * L / diffCoeff)
>>> analyticalArray = AA + BB * CC / DD
>>> print var.allclose(analyticalArray, rtol=1e-4, atol=1e-4)
1

If the problem is run interactively, we can view the result:

.. index::
   module: fipy.viewers

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var)
...     viewer.plot()
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')
