#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "source.py"
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

r"""Solve a convection problem with a source.

This example solves the equation

.. math::

   \frac{\partial \phi}{\partial x} + \alpha \phi = 0
   
with :math:`\phi \left( 0 \right) = 1` at :math:`x = 0`. The boundary
condition at :math:`x = L` is an outflow boundary condition requiring
the use of an artificial constraint to be set on the right hand side
faces. Exterior faces without constraints are considered to have zero
outflow. An :class:`~fipy.terms.implicitSourceTerm.ImplicitSourceTerm`
object will be used to represent this term. The derivative of
:math:`\phi` can be represented by a
:class:`~fipy.terms.ConvectionTerm` with a constant unitary velocity
field from left to right. The following is an example code that
includes a test against the analytical result.

>>> from fipy import *

>>> L = 10.
>>> nx = 5000
>>> dx =  L / nx
>>> mesh = Grid1D(dx=dx, nx=nx)
>>> phi0 = 1.0
>>> alpha = 1.0
>>> phi = CellVariable(name=r"$\phi$", mesh=mesh, value=phi0)
>>> solution = CellVariable(name=r"solution", mesh=mesh, value=phi0 * numerix.exp(-alpha * mesh.cellCenters[0]))

>>> if __name__ == "__main__":
...     viewer = Viewer(vars=(phi, solution))
...     viewer.plot()
...     raw_input("press key to continue")

>>> phi.constrain(phi0, mesh.facesLeft)
>>> ## fake outflow condition
>>> phi.faceGrad.constrain([0], mesh.facesRight)

>>> eq = PowerLawConvectionTerm((1,)) + ImplicitSourceTerm(alpha)
>>> eq.solve(phi)
>>> print numerix.allclose(phi, phi0 * numerix.exp(-alpha * mesh.cellCenters[0]), atol=1e-3)
True

>>> if __name__ == "__main__":
...     viewer = Viewer(vars=(phi, solution))
...     viewer.plot()
...     raw_input("finished")    


"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
