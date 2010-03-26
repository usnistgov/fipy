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

r"""

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
          \left[ \left(P + A \right) \exp{\left(\frac{A}{2} \left(x - 1\right)\right)} -
                 \left(P - A \right) \exp{\left(-\frac{A}{2} \left(x - 1\right)\right)} \right]}
        { \left(P + A \right)^2 \exp{\left(\frac{A}{2}\right)} -
          \left(P - A \right)^2 \exp{\left(-\frac{A}{2}\right)}}

where

.. math::
    
   A = \sqrt{P + 4D^2}
   
..

>>> from fipy import *
>>> nx = 100
>>> dx = 1.0 / nx

>>> mesh = Grid1D(nx=nx, dx=dx)
>>> C = CellVariable(mesh=mesh)

>>> D = 2.0
>>> P = 3.0

From the main equation, the flux into the domain at :math:`x=0` is given
by

.. math::
    
   \frac{\partial C}{\partial x}
    - P C

Using the boundary condition at :math:`x=0` this flux should be equal to
:math:`-P`. Setting the :math:`x=1` boundary condition to be a fixed value
equal to :math:`C \left( 1 \right)` fixes the edge derivative on both the
convection and diffusion terms to be zero.

>>> BCs = (FixedFlux(faces=mesh.getFacesLeft(), value=-P),
...        FixedValue(faces=mesh.getFacesRight(), value=C.getFaceValue()))

>>> eq = PowerLawConvectionTerm((P,)) == \
...      DiffusionTerm() - ImplicitSourceTerm(D)

>>> A = numerix.sqrt(P**2 + 4 * D)

>>> x = mesh.getCellCenters()[0]
>>> CAnalytical = CellVariable(mesh=mesh)
>>> CAnalytical.setValue(2 * P * exp(P * x / 2) * ((P + A) * exp(A / 2 * (1 - x))
...             - (P - A) * exp(-A / 2 *(1 - x)))/
...             ((P + A)**2*exp(A / 2)- (P - A)**2 * exp(-A / 2)))

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=(C, CAnalytical))

>>> res = 1e+10
>>> while res > 1e-5:
...     res = eq.sweep(var=C, boundaryConditions=BCs)
...     if __name__ == '__main__':
...         viewer.plot()

>>> print C.allclose(CAnalytical, rtol=1.e-3, atol=1.e-3)
True

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input('finished')
