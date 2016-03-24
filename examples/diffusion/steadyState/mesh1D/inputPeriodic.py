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

r"""

One can then solve the same problem as in
`examples/diffusion/steadyState/mesh1D/input.py` but with a periodic
mesh and no boundary conditions. The periodic mesh is used to simulate
periodic boundary conditions.

>>> from fipy import PeriodicGrid1D, CellVariable, TransientTerm, DiffusionTerm, Viewer

>>> nx = 50
>>> dx = 1.
>>> mesh = PeriodicGrid1D(nx = nx, dx = dx)

The variable is initially a line varying form `valueLeft` to `valueRight`.

>>> valueLeft = 0
>>> valueRight = 1
>>> x = mesh.cellCenters[0]

>>> Lx = nx * dx
>>> initialArray = valueLeft + (valueRight - valueLeft) * x / Lx
>>> var = CellVariable(name = "solution variable", mesh = mesh,
...                                                value = initialArray)

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var, datamin=0., datamax=1.)
...     viewer.plot()
...     raw_input("press key to continue")


A `TransientTerm` is used to provide some fixed point, otherwise the
solver has no fixed value and can become unstable.

>>> eq = TransientTerm(coeff=1e-8) - DiffusionTerm()
>>> eq.solve(var=var, dt=1.)

>>> if __name__ == '__main__':
...     viewer.plot()

The result of the calculation will be the average value over the domain.

>>> print var.allclose((valueLeft + valueRight) / 2., rtol = 1e-5)
1

"""

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    input("finished")
