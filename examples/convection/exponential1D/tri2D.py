#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "tri2D.py"
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

This example solves the steady-state convection-diffusion equation as described in
:mod:`examples.diffusion.convection.exponential1D.input` but uses a
:class:`~fipy.meshes.tri2D.Tri2D` mesh.

Here the axes are reversed (``nx = 1``, ``ny = 1000``) and

.. math::

   \vec{u} = (0, 10)

>>> from fipy import CellVariable, Tri2D, DiffusionTerm, ExponentialConvectionTerm, DefaultAsymmetricSolver, Viewer
>>> from fipy.tools import numerix

>>> L = 10.
>>> nx = 1
>>> ny = 1000
>>> mesh = Tri2D(dx = L / ny, dy = L / ny, nx = nx, ny = ny)

>>> valueBottom = 0.
>>> valueTop = 1.

>>> var = CellVariable(name = "concentration",
...                    mesh = mesh,
...                    value = valueBottom)

>>> var.constrain(valueBottom, mesh.facesBottom)
>>> var.constrain(valueTop, mesh.facesTop)

>>> diffCoeff = 1.
>>> convCoeff = numerix.array(((0.,), (10.,)))

>>> eq = (DiffusionTerm(coeff=diffCoeff)
...       + ExponentialConvectionTerm(coeff=convCoeff))

>>> eq.solve(var = var,
...          solver=DefaultAsymmetricSolver(iterations=10000))

The analytical solution test for this problem is given by:

>>> axis = 1
>>> y = mesh.cellCenters[axis]
>>> CC = 1. - numerix.exp(-convCoeff[axis] * y / diffCoeff)
>>> DD = 1. - numerix.exp(-convCoeff[axis] * L / diffCoeff)
>>> analyticalArray = CC / DD

>>> print var.allclose(analyticalArray, rtol = 1e-6, atol = 1e-6)
1

>>> if __name__ == '__main__':
...     viewer = Viewer(vars = var)
...     viewer.plot()
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input('finished')
