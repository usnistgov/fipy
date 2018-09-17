#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "peclet.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

r"""

This example tests diffusion-convection for increasing Peclet numbers.
This test case has been introduced because :class:`~fipy.solvers.pysparse.linearCGSSolver.LinearCGSSolver` was not
working with Peclet numbers over 1. :class:`~fipy.solvers.pysparse.linearLUSolver.LinearLUSolver` is now the default
for :class:`~fipy.terms.convectionTerm.ConvectionTerm`. For ``nx = 1000`` the :class:`~fipy.solvers.pysparse.linearGMRESSolver.LinearGMRESSolver` does not work.

>>> from fipy import CellVariable, Grid1D, TransientTerm, DiffusionTerm, PowerLawConvectionTerm, Viewer
>>> from fipy.tools import numerix

>>> L = 1.
>>> nx = 1000
>>> dx =  L / nx
>>> mesh = Grid1D(dx=dx , nx=nx)

>>> valueLeft = 0.
>>> valueRight = 1.

>>> var = CellVariable(name = "solution variable", mesh=mesh, value=valueLeft)

>>> var.constrain(valueLeft, mesh.facesLeft)
>>> var.constrain(valueRight, mesh.facesRight)

>>> if __name__ == '__main__':
...     viewer = Viewer(vars = var)

>>> convCoeff = 1.0
>>> peclet = 1e-3
>>> allcloseList = []
>>> while peclet < 1e4:
...     var[:] = valueLeft
...     diffCoeff = convCoeff * dx / peclet
...     eq = (TransientTerm(1e-4)
...           == DiffusionTerm(coeff=diffCoeff)
...           + PowerLawConvectionTerm(coeff=(convCoeff,)))
...     eq.solve(var=var, dt=1.)
...     x = mesh.cellCenters[0]
...     arg0 = -convCoeff * x / diffCoeff
...     arg0 = numerix.where(arg0 < -200, -200, arg0)
...     arg1 = -convCoeff * L / diffCoeff
...     arg1 = (arg1 >= -200) * (arg1 + 200) - 200
...     CC = 1. - numerix.exp(arg0)
...     DD = 1. - numerix.exp(arg1)
...     analyticalArray = CC / DD
...     allcloseList.append(var.allclose(CC / DD, rtol = 1e-2, atol = 1e-2).value)
...     if __name__ == '__main__':
...         viewer.plot()
...         raw_input("Peclet number: " + str(peclet) + ", press key")
...     peclet *= 10

>>> print allcloseList
[True, True, True, True, True, True, True]

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
