#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "circle.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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

r"""Solve an advection-diffusion equation with a Robin boundary condition.

This example demonstrates how to apply a Robin boundary condition to
an advection-diffusion equation. The equation we wish to solve is
given by,

.. math::

   0 &=\frac{\partial^2 C}{\partial x^2}
    - P \frac{\partial C}{\partial x} -D C \qquad 0 < x <1\\
   x=0: P &= -\frac{\partial C}{\partial x} + P C \\
   x=1: \frac{\partial C}{\partial x} & = 0

The analytical solution for this equation is given by,

.. math::

   C \left( x \right) =
   \frac{ 2 P \exp{\left(\frac{P x}{2}\right)}
          \left[ \left(P + A \right) \exp{\left(\frac{A}{2} \left(1 - x\right)\right)} -
                 \left(P - A \right) \exp{\left(-\frac{A}{2} \left(1 - x\right)\right)} \right]}
        { \left(P + A \right)^2 \exp{\left(\frac{A}{2}\right)} -
          \left(P - A \right)^2 \exp{\left(-\frac{A}{2}\right)}}

where

.. math::

   A = \sqrt{P^2 + 4D}

..

>>> from fipy import CellVariable, Grid1D, DiffusionTerm, PowerLawConvectionTerm, ImplicitSourceTerm, Viewer
>>> from fipy.tools import numerix
>>> nx = 100
>>> dx = 1.0 / nx

>>> mesh = Grid1D(nx=nx, dx=dx)
>>> C = CellVariable(mesh=mesh)

>>> D = 2.0
>>> P = 3.0

>>> C.faceGrad.constrain([-P + P * C.faceValue], mesh.facesLeft)
>>> C.faceGrad.constrain([0], mesh.facesRight)

>>> eq = PowerLawConvectionTerm((P,)) == \
...      DiffusionTerm() - ImplicitSourceTerm(D)

>>> A = numerix.sqrt(P**2 + 4 * D)

>>> x = mesh.cellCenters[0]
>>> CAnalytical = CellVariable(mesh=mesh)
>>> CAnalytical.setValue(2 * P * numerix.exp(P * x / 2) * ((P + A) * numerix.exp(A / 2 * (1 - x))
...             - (P - A) * numerix.exp(-A / 2 *(1 - x)))/
...             ((P + A)**2*numerix.exp(A / 2)- (P - A)**2 * numerix.exp(-A / 2)))

>>> if __name__ == '__main__':
...     C.name = 'C'
...     viewer = Viewer(vars=(C, CAnalytical))

>>> if __name__ == '__main__':
...     restol = 1e-5
...     anstol = 1e-3
... else:
...     restol = 0.5
...     anstol = 0.15

>>> res = 1e+10

>>> while res > restol:
...     res = eq.sweep(var=C)
...     if __name__ == '__main__':
...         viewer.plot()

>>> print C.allclose(CAnalytical, rtol=anstol, atol=anstol)
True

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')
