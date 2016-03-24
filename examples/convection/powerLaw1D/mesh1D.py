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

"""

This example solves the steady-state convection-diffusion equation as
described in :mod:`examples.diffusion.convection.exponential1D.mesh1D` but
uses the :class:`~fipy.terms.powerLawConvectionTerm.PowerLawConvectionTerm` rather than the
:class:`~fipy.terms.exponentialConvectionTerm.ExponentialConvectionTerm`.

>>> from fipy import CellVariable, Grid1D, DiffusionTerm, PowerLawConvectionTerm, Viewer
>>> from fipy.tools import numerix

>>> L = 10.
>>> nx = 1000
>>> mesh = Grid1D(dx = L / nx, nx = nx)

>>> valueLeft = 0.
>>> valueRight = 1.

>>> var = CellVariable(name = "concentration",
...                    mesh = mesh,
...                    value = valueLeft)

>>> var.constrain(valueLeft, mesh.facesLeft)
>>> var.constrain(valueRight, mesh.facesRight)

>>> diffCoeff = 1.
>>> convCoeff = (10.,)

>>> eq = (DiffusionTerm(coeff=diffCoeff)
...       + PowerLawConvectionTerm(coeff=convCoeff))

>>> eq.solve(var = var)

We test the solution against the analytical result:

>>> axis = 0
>>> x = mesh.cellCenters[axis]
>>> CC = 1. - numerix.exp(-convCoeff[axis] * x / diffCoeff)
>>> DD = 1. - numerix.exp(-convCoeff[axis] * L / diffCoeff)
>>> analyticalArray = CC / DD
>>> print var.allclose(analyticalArray, rtol = 1e-2, atol = 1e-2)
1

If the problem is run interactively, we can view the result:

>>> if __name__ == '__main__':
...     viewer = Viewer(vars = var)
...     viewer.plot()
"""

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')
